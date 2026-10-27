"""Bounded OpenAI-compatible model calls with strict choice caching."""

from __future__ import annotations

import errno
import hashlib
import hmac
import http.client
import json
import logging
import math
import os
import re
import socket
import ssl
import threading
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable, Iterator, Mapping
from contextlib import contextmanager, suppress
from pathlib import Path

from mendelprime.domain import MendelPrimeError, validate_identifier
from mendelprime.runtime.http import tls_context
from mendelprime.runtime.storage import atomic_replace_bytes, private_directory, read_regular_bounded

LOGGER = logging.getLogger(__name__)
MAX_RESPONSE_BYTES = 8 * 1024 * 1024
LLM_DEFAULTS = {
    "model": "gpt-6.1-sol",
    "api_key_env": "MENDELPRIME_API_KEY",
    "wire_api": "chat_completions",
    "transport_retries": 1,
}
_MAX_CACHE_BYTES = 64 * 1024
_CACHE_VERSION = 1
_KEY_PATTERN = re.compile(r"[0-9a-f]{64}")
_STAGE_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.:-]*")
_ENVIRONMENT_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
_REPORTED_MODEL = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:/-]{0,255}")
_RETRYABLE_HTTP_STATUSES = {408, 409, 429, 500, 502, 503, 504}
_RETRYABLE_ERRNOS = {
    value
    for name in (
        "ECONNABORTED",
        "ECONNREFUSED",
        "ECONNRESET",
        "EHOSTDOWN",
        "EHOSTUNREACH",
        "ENETDOWN",
        "ENETRESET",
        "ENETUNREACH",
        "EPIPE",
        "ETIMEDOUT",
    )
    if (value := getattr(errno, name, None)) is not None
}


def _contract(field: str, code: str) -> MendelPrimeError:
    return MendelPrimeError(code, category="contract", field=field)


def _positive_seconds(value: object, field: str) -> float:
    if type(value) not in {int, float} or not math.isfinite(value) or value <= 0:
        raise _contract(field, "positive_finite_number_required")
    return float(value)


def _has_control(value: str) -> bool:
    return any(unicodedata.category(character) == "Cc" for character in value)


def _plain_nonempty(value: object) -> bool:
    return (
        isinstance(value, str) and bool(value.strip()) and value == value.strip() and not _has_control(value)
    )


def canonical_base_url(value: object) -> str:
    invalid_text = (
        not isinstance(value, str)
        or not value
        or _has_control(value)
        or any(character.isspace() for character in value)
        or "\\" in value
        or "?" in value
        or "#" in value
    )
    if invalid_text:
        raise _contract("llm.base_url", "absolute_api_root_required")
    try:
        parsed = urllib.parse.urlsplit(value)
        port = parsed.port
    except ValueError:
        raise _contract("llm.base_url", "absolute_api_root_required") from None
    scheme = parsed.scheme.lower()
    invalid_url = (
        scheme not in {"http", "https"}
        or not parsed.netloc
        or parsed.hostname is None
        or parsed.username is not None
        or parsed.password is not None
        or parsed.query
        or parsed.fragment
    )
    if invalid_url:
        raise _contract("llm.base_url", "absolute_api_root_required")
    try:
        host = parsed.hostname.encode("idna").decode("ascii").lower()
    except UnicodeError:
        raise _contract("llm.base_url", "absolute_api_root_required") from None
    if not host or _has_control(host):
        raise _contract("llm.base_url", "absolute_api_root_required")
    default_port = 80 if scheme == "http" else 443
    rendered_host = f"[{host}]" if ":" in host else host
    authority = rendered_host if port in {None, default_port} else f"{rendered_host}:{port}"
    path = parsed.path.rstrip("/")
    if path.endswith(("/chat/completions", "/responses")):
        raise _contract("llm.base_url", "api_root_not_endpoint_required")
    return urllib.parse.urlunsplit((scheme, authority, path, "", ""))


def _validated_settings(
    settings: object, *, wire_code="wire_api_required", config_order=False
) -> dict[str, object]:
    if not isinstance(settings, Mapping):
        raise _contract("llm", "mapping_required")
    allowed = {"base_url", *LLM_DEFAULTS}
    if set(settings) - allowed:
        raise _contract("llm", "unknown_field")
    if "base_url" not in settings:
        raise _contract("llm.base_url", "required")
    values = LLM_DEFAULTS | dict(settings)
    model, key_env = values["model"], values["api_key_env"]
    wire_api, retries = values["wire_api"], values["transport_retries"]
    invalid_wire = wire_api not in {"chat_completions", "responses"}
    if config_order and invalid_wire:
        raise _contract("llm.wire_api", wire_code)
    base_url = canonical_base_url(settings["base_url"]) if config_order else None
    if not _plain_nonempty(model):
        raise _contract("llm.model", "nonempty_text_required")
    if not isinstance(key_env, str) or not _ENVIRONMENT_NAME.fullmatch(key_env):
        raise _contract("llm.api_key_env", "environment_name_required")
    if not config_order and invalid_wire:
        raise _contract("llm.wire_api", wire_code)
    if type(retries) is not int or not 0 <= retries <= 3:
        raise _contract("llm.transport_retries", "integer_0_to_3_required")
    return {
        "base_url": base_url if config_order else canonical_base_url(settings["base_url"]),
        "model": model,
        "api_key_env": key_env,
        "wire_api": wire_api,
        "transport_retries": retries,
    }


def _validated_request(request: object) -> dict[str, object]:
    if not isinstance(request, dict):
        raise _contract("request", "mapping_required")
    required = {"case_key", "stage_id", "system_prompt", "user_prompt", "max_tokens", "timeout_seconds"}
    if not required <= set(request):
        raise _contract("request", "missing_field")
    if set(request) - required - {"cache_choice", "max_http_attempts"}:
        raise _contract("request", "unknown_field")
    case_key = request["case_key"]
    if type(case_key) is not tuple or len(case_key) != 2:
        raise _contract("case_key", "case_key_required")
    validate_identifier(case_key[0], "cohort")
    validate_identifier(case_key[1], "case_id")
    stage = request["stage_id"]
    if not isinstance(stage, str) or not _STAGE_PATTERN.fullmatch(stage):
        raise _contract("stage_id", "stage_identifier_required")
    if not isinstance(request["system_prompt"], str):
        raise _contract("system_prompt", "text_required")
    user_prompt = request["user_prompt"]
    if not isinstance(user_prompt, str) or not user_prompt.strip():
        raise _contract("user_prompt", "nonempty_text_required")
    if type(request["max_tokens"]) is not int or request["max_tokens"] < 1:
        raise _contract("max_tokens", "positive_integer_required")
    timeout = request["timeout_seconds"]
    _positive_seconds(timeout, "timeout_seconds")
    cache_choice = request.get("cache_choice", False)
    if type(cache_choice) is not bool:
        raise _contract("cache_choice", "boolean_required")
    if "max_http_attempts" in request:
        cap = request["max_http_attempts"]
        if type(cap) is not int or not 1 <= cap <= 4:
            raise _contract("max_http_attempts", "integer_1_to_4_required")
    return {
        "case_key": case_key,
        "stage_id": stage,
        "system_prompt": request["system_prompt"],
        "user_prompt": user_prompt,
        "max_tokens": request["max_tokens"],
        "timeout_seconds": timeout,
        "cache_choice": cache_choice,
        **({"max_http_attempts": request["max_http_attempts"]} if "max_http_attempts" in request else {}),
    }


def _canonical_json(value: object) -> bytes:
    return json.dumps(
        value, ensure_ascii=True, allow_nan=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def _endpoint(settings: Mapping[str, object]) -> str:
    suffix = "responses" if settings["wire_api"] == "responses" else "chat/completions"
    return f"{settings['base_url']}/{suffix}"


def _wire_body(settings: Mapping[str, object], request: Mapping[str, object]) -> bytes:
    if settings["wire_api"] == "responses":
        payload = {
            "model": settings["model"],
            "instructions": request["system_prompt"],
            "input": request["user_prompt"],
            "max_output_tokens": request["max_tokens"],
            "store": False,
        }
    else:
        payload = {
            "model": settings["model"],
            "messages": [
                {"role": "system", "content": request["system_prompt"]},
                {"role": "user", "content": request["user_prompt"]},
            ],
            "temperature": 0,
            "max_tokens": request["max_tokens"],
            "store": False,
        }
    return json.dumps(payload, ensure_ascii=True, separators=(",", ":")).encode("ascii")


def _provider_identity(settings: Mapping[str, object], api_key_bytes: bytes) -> str:
    fields = {
        "credential_scope": hashlib.sha256(api_key_bytes).hexdigest(),
        "endpoint": _endpoint(settings),
        "model": settings["model"],
        "version": 1,
        "wire_api": settings["wire_api"],
    }
    return hashlib.sha256(_canonical_json(fields)).hexdigest()


def _request_key(provider_identity: str, request: Mapping[str, object]) -> str:
    cohort, case_id = request["case_key"]
    value = {
        "case": {"case_id": case_id, "cohort": cohort},
        "max_tokens": request["max_tokens"],
        "provider_identity": provider_identity,
        "stage_id": request["stage_id"],
        "system_prompt": request["system_prompt"],
        "timeout_seconds": float(request["timeout_seconds"]),
        "user_prompt": request["user_prompt"],
        "version": _CACHE_VERSION,
        **({"max_http_attempts": request["max_http_attempts"]} if "max_http_attempts" in request else {}),
    }
    return hashlib.sha256(_canonical_json(value)).hexdigest()


def _strict_choice(content: object) -> bool:
    return isinstance(content, str) and content.strip().upper() in {"A", "B"}


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError("duplicate cache field")
        value[key] = item
    return value


def _cache_load(root: Path, key: str) -> tuple[dict[str, object] | None, str]:
    try:
        private_directory(root, create=False)
        raw, _before, _after = read_regular_bounded(
            root / f"{key}.json", _MAX_CACHE_BYTES, follow_symlinks=False
        )
    except FileNotFoundError:
        return None, "miss"
    except OverflowError:
        return None, "rejected"
    except OSError as error:
        rejected = error.errno == errno.ELOOP or ("regular" in str(error).lower() and root.is_dir())
        return None, "rejected" if rejected else "unavailable"
    try:
        record = json.loads(raw.decode("utf-8"), object_pairs_hook=_unique_object)
        fields = {"version", "key", "content", "reported_model", "payload_sha256"}
        valid_record = (
            isinstance(record, dict)
            and set(record) == fields
            and type(record["version"]) is int
            and record["version"] == _CACHE_VERSION
            and record["key"] == key
            and _strict_choice(record["content"])
        )
        if not valid_record:
            return None, "rejected"
        checksum = record["payload_sha256"]
        payload = {name: record[name] for name in ("version", "key", "content", "reported_model")}
        expected = hashlib.sha256(_canonical_json(payload)).hexdigest()
        if not isinstance(checksum, str) or not _KEY_PATTERN.fullmatch(checksum):
            return None, "rejected"
        if not hmac.compare_digest(checksum, expected):
            return None, "rejected"
        model = record["reported_model"]
        if model is not None and (not isinstance(model, str) or not model.strip()):
            return None, "rejected"
    except (RecursionError, UnicodeDecodeError, ValueError, TypeError):
        return None, "rejected"
    return {"content": record["content"], "reported_model": model, "http_attempts": 0}, "hit"


def _cache_save(root: Path, key: str, reply: Mapping[str, object]) -> bool:
    payload = {
        "version": _CACHE_VERSION,
        "key": key,
        "content": reply["content"],
        "reported_model": reply["reported_model"],
    }
    record = payload | {"payload_sha256": hashlib.sha256(_canonical_json(payload)).hexdigest()}
    raw = _canonical_json(record)
    if len(raw) > _MAX_CACHE_BYTES:
        raise _contract("cache.reply", "cache_record_too_large")
    try:
        private_directory(root)
        atomic_replace_bytes(root / f"{key}.json", raw, prefix=f".{key}.")
        return True
    except OSError:
        return False


def _response_error(code: str, attempts: int) -> MendelPrimeError:
    return MendelPrimeError(code, category="response", attempts=attempts)


def _require_response(condition: bool, attempts: int) -> None:
    if not condition:
        raise _response_error("invalid_response", attempts)


def _text_parts(parts: object, kind: str, attempts: int) -> str:
    _require_response(isinstance(parts, list), attempts)
    text = []
    for part in parts:
        _require_response(isinstance(part, dict), attempts)
        if part.get("type") == kind:
            value = part.get("text")
            _require_response(isinstance(value, str), attempts)
            text.append(value)
    return "".join(text)


def _response_content(payload: object, wire_api: str, attempts: int) -> str:
    _require_response(isinstance(payload, dict), attempts)
    if wire_api == "responses":
        if "status" in payload and payload["status"] != "completed":
            raise _response_error("incomplete_response", attempts)
        if "output_text" in payload:
            content = payload["output_text"]
            _require_response(isinstance(content, str), attempts)
        else:
            output = payload.get("output")
            _require_response(isinstance(output, list), attempts)
            text: list[str] = []
            for item in output:
                _require_response(isinstance(item, dict), attempts)
                if item.get("type") != "message":
                    continue
                text.append(_text_parts(item.get("content"), "output_text", attempts))
            content = "".join(text)
    else:
        choices = payload.get("choices")
        valid_choices = isinstance(choices, list) and choices and isinstance(choices[0], dict)
        _require_response(bool(valid_choices), attempts)
        choice = choices[0]
        if "finish_reason" in choice and choice["finish_reason"] != "stop":
            raise _response_error("incomplete_response", attempts)
        message = choice.get("message")
        _require_response(isinstance(message, dict) and "content" in message, attempts)
        if isinstance(message["content"], str):
            content = message["content"]
        elif isinstance(message["content"], list):
            content = _text_parts(message["content"], "text", attempts)
        else:
            raise _response_error("invalid_response", attempts)
    if not content.strip():
        raise _response_error("empty_content", attempts)
    return content


def _parse_reply(raw: bytes, wire_api: str, attempts: int) -> dict[str, object]:
    if len(raw) > MAX_RESPONSE_BYTES:
        raise _response_error("response_too_large", attempts)
    if not raw:
        raise _response_error("empty_response", attempts)
    try:
        text = raw.decode("utf-8")
        payload = json.loads(text)
        pairs = json.loads(text, object_pairs_hook=lambda value: value)
    except (UnicodeDecodeError, json.JSONDecodeError, RecursionError):
        raise _response_error("invalid_response_json", attempts) from None
    content = _response_content(payload, wire_api, attempts)
    model = payload.get("model") if isinstance(payload, dict) else None
    duplicate_model = (
        isinstance(pairs, list)
        and sum(1 for item in pairs if isinstance(item, tuple) and item[0] == "model") > 1
    )
    if duplicate_model or not isinstance(model, str) or not _REPORTED_MODEL.fullmatch(model):
        model = None
    return {"content": content, "reported_model": model, "http_attempts": attempts}


def _retryable(error: Exception) -> bool:
    if isinstance(error, urllib.error.HTTPError):
        return error.code in _RETRYABLE_HTTP_STATUSES
    reason = error.reason if isinstance(error, urllib.error.URLError) else error
    if isinstance(reason, (ssl.CertificateError, ssl.SSLError)):
        return False
    if isinstance(reason, (TimeoutError, socket.timeout)):
        return True
    if isinstance(
        reason,
        (ConnectionError, socket.gaierror, http.client.IncompleteRead, http.client.RemoteDisconnected),
    ):
        return True
    return isinstance(reason, OSError) and reason.errno in _RETRYABLE_ERRNOS


def _transport_error(error: Exception, attempts: int) -> MendelPrimeError:
    retryable = _retryable(error)
    if isinstance(error, urllib.error.HTTPError):
        code = "http_redirect" if 300 <= error.code <= 399 else "http_error"
        return MendelPrimeError(
            code, category="model", attempts=attempts, http_status=error.code, retryable=retryable
        )
    reason = error.reason if isinstance(error, urllib.error.URLError) else error
    if isinstance(reason, (ssl.CertificateError, ssl.SSLError)):
        code = "tls_error"
    elif isinstance(reason, (TimeoutError, socket.timeout)) or (
        isinstance(reason, OSError) and reason.errno == getattr(errno, "ETIMEDOUT", None)
    ):
        code = "transport_timeout"
    elif retryable:
        code = "connection_error"
    else:
        code = "transport_error"
    return MendelPrimeError(code, category="model", attempts=attempts, retryable=retryable)


class RequestBudget:
    """Bound concurrent model calls within one process and one application run."""

    def __init__(self, maximum: int) -> None:
        if type(maximum) is not int or maximum < 1:
            raise _contract("maximum", "positive_integer_required")
        self._semaphore = threading.BoundedSemaphore(maximum)

    @contextmanager
    def slot(self, timeout_seconds: float) -> Iterator[float]:
        timeout = _positive_seconds(timeout_seconds, "timeout_seconds")
        started = time.monotonic()
        if not self._semaphore.acquire(timeout=timeout):
            raise MendelPrimeError("budget_timeout", category="model", retryable=True)
        try:
            yield time.monotonic() - started
        finally:
            self._semaphore.release()


class LLMClient:
    """Validate, cache, budget, send, and audit one provider configuration."""

    def __init__(
        self,
        settings: Mapping[str, object],
        budget: RequestBudget,
        queue_timeout_seconds: float = 1800,
        *,
        cache_root: Path | None = None,
        observer: Callable[[dict[str, object]], None] | None = None,
        responses: Callable[[dict[str, object], dict[str, object]], None] | None = None,
        environ: Mapping[str, str] | None = None,
    ) -> None:
        if not isinstance(budget, RequestBudget):
            raise _contract("budget", "request_budget_required")
        if cache_root is not None and not isinstance(cache_root, Path):
            raise _contract("cache_root", "path_or_none_required")
        if observer is not None and not callable(observer):
            raise _contract("observer", "callable_or_none_required")
        if responses is not None and not callable(responses):
            raise _contract("responses", "callable_or_none_required")
        values = os.environ if environ is None else environ
        if not isinstance(values, Mapping):
            raise _contract("environ", "mapping_required")
        self._settings = _validated_settings(settings)
        api_key = values.get(self._settings["api_key_env"])
        try:
            key_bytes = api_key.encode("utf-8") if isinstance(api_key, str) else None
        except UnicodeEncodeError:
            key_bytes = None
        invalid_key = (
            not isinstance(api_key, str) or not api_key.strip() or _has_control(api_key) or key_bytes is None
        )
        if invalid_key:
            raise _contract("llm.api_key_env", "valid_api_key_required")
        self._api_key = api_key
        self._identity = _provider_identity(self._settings, key_bytes)
        self._budget = budget
        self._queue_timeout = _positive_seconds(queue_timeout_seconds, "queue_timeout_seconds")
        self._cache_root = cache_root
        self._observer = observer
        self._responses = responses
        redirect = urllib.request.HTTPRedirectHandler()
        redirect.max_repeats = 0
        redirect.max_redirections = 0
        context = tls_context()
        handlers = [redirect] if context is None else [redirect, urllib.request.HTTPSHandler(context=context)]
        self._opener = urllib.request.build_opener(*handlers)

    def complete(self, request: dict[str, object]) -> dict[str, object]:
        value = _validated_request(request)
        started = time.monotonic()
        request_id = _request_key(self._identity, value)
        input_characters = len(value["system_prompt"]) + len(value["user_prompt"])
        audit = value, request_id, input_characters, started
        cache_status = "disabled"
        if value["cache_choice"] and self._cache_root is not None:
            try:
                cached, cache_status = _cache_load(self._cache_root, request_id)
            except Exception:
                LOGGER.warning("choice cache read failed")
                cached, cache_status = None, "unavailable"
            if cached is not None and cache_status == "hit":
                self._emit(*audit, "success", "hit", cached)
                self._observe_response(value, cached)
                return cached
        queue_started = time.monotonic()
        queue_seconds = 0.0
        provider_seconds = 0.0
        acquired = False
        failure = None
        try:
            with self._budget.slot(self._queue_timeout) as queue_seconds:
                acquired = True
                provider_started = time.monotonic()
                try:
                    reply = self._send(value)
                finally:
                    provider_seconds = time.monotonic() - provider_started
        except Exception as caught:
            failure = (
                caught
                if isinstance(caught, MendelPrimeError)
                else MendelPrimeError("transport_error", category="model")
            )
        if failure is not None:
            if not acquired:
                queue_seconds = time.monotonic() - queue_started
            self._emit(*audit, "error", cache_status, None, queue_seconds, provider_seconds, failure)
            raise failure from None
        if value["cache_choice"] and self._cache_root is not None:
            cache_status = self._store_choice(request_id, reply, cache_status)
        self._emit(*audit, "success", cache_status, reply, queue_seconds, provider_seconds)
        self._observe_response(value, reply)
        return reply

    def _store_choice(self, request_id: str, reply: Mapping[str, object], cache_status: str) -> str:
        if not _strict_choice(reply["content"]):
            return "rejected"
        try:
            if _cache_save(self._cache_root, request_id, reply):
                return cache_status
        except MendelPrimeError:
            return "rejected"
        except Exception:
            pass
        LOGGER.warning("choice cache write failed")
        return "write_failed"

    def _send(self, request: Mapping[str, object]) -> dict[str, object]:
        body = _wire_body(self._settings, request)
        maximum_attempts = min(self._settings["transport_retries"] + 1, request.get("max_http_attempts", 4))
        for attempt in range(1, maximum_attempts + 1):
            failure = None
            try:
                outgoing = urllib.request.Request(
                    _endpoint(self._settings),
                    data=body,
                    headers={"Authorization": f"Bearer {self._api_key}", "Content-Type": "application/json"},
                    method="POST",
                )
                # The stdlib checks zero redirect limits only after history exists.
                outgoing.redirect_dict = {}
                with self._opener.open(outgoing, timeout=request["timeout_seconds"]) as response:
                    raw = response.read(MAX_RESPONSE_BYTES + 1)
                return _parse_reply(raw, self._settings["wire_api"], attempt)
            except MendelPrimeError as error:
                failure = error
            except Exception as error:
                retryable = _retryable(error)
                if isinstance(error, urllib.error.HTTPError):
                    with suppress(Exception):
                        error.close()
                if retryable and attempt < maximum_attempts:
                    time.sleep(min(2 ** (attempt - 1), 8))
                    continue
                failure = _transport_error(error, attempt)
            raise failure from None
        raise MendelPrimeError("transport_error", category="model", attempts=maximum_attempts) from None

    def _emit(
        self,
        request: Mapping[str, object],
        request_id: str,
        input_characters: int,
        started: float,
        status: str,
        cache_status: str,
        reply: Mapping[str, object] | None,
        queue_seconds: float = 0.0,
        provider_seconds: float = 0.0,
        error: MendelPrimeError | None = None,
    ) -> None:
        if self._observer is None:
            return
        event = {
            "case_key": request["case_key"],
            "stage_id": request["stage_id"],
            "request_id": request_id,
            "input_characters": input_characters,
            "status": status,
            "cache_status": cache_status,
            "http_attempts": error.attempts if error else reply["http_attempts"],
            "queue_seconds": queue_seconds,
            "provider_seconds": provider_seconds,
            "elapsed_seconds": time.monotonic() - started,
            "error_code": error.code if error else None,
            "http_status": error.http_status if error else None,
        }
        try:
            self._observer(event)
        except Exception:
            LOGGER.warning("model progress callback failed")

    def _observe_response(self, request: Mapping[str, object], reply: Mapping[str, object]) -> None:
        if self._responses is None:
            return
        try:
            self._responses(dict(request), dict(reply))
        except Exception:
            LOGGER.warning("model response observer failed")


__all__ = ["LLMClient", "LLM_DEFAULTS", "RequestBudget", "canonical_base_url"]
