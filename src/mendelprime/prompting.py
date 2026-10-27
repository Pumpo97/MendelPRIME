"""Frozen prompt assets and structured calls with caller-owned recovery policy."""

import json
import re
from functools import lru_cache
from importlib.resources import files

from .domain import MendelPrimeError


def text(name):
    return _asset(name, ".txt")


@lru_cache(maxsize=None)
def _asset(name, extension):
    for section in ("ranking", "explanation"):
        path = files("mendelprime").joinpath("prompts", section, name + extension)
        if path.is_file():
            return path.read_text(encoding="utf-8")
    raise MendelPrimeError(
        "prompt_not_found" if extension == ".txt" else "template_not_found", field="prompt"
    )


def template(name, key=None):
    value = json.loads(_asset(name, ".json"))
    return value if key is None else value[key]


def parse_object(value):
    cleaned = re.sub(r"^```(?:json)?\s*|\s*```$", "", value.strip(), flags=re.IGNORECASE)
    for candidate in (cleaned, cleaned[cleaned.find("{") : cleaned.rfind("}") + 1]):
        try:
            result = json.loads(candidate)
        except (json.JSONDecodeError, RecursionError):
            continue
        if isinstance(result, dict):
            return result
    raise MendelPrimeError("invalid_json_object", category="response")


def call_json(model, request, decode, attempts=1, repair=None):
    code = "not_attempted"
    for attempt in range(attempts):
        current = dict(request)
        if attempt and repair is not None:
            current["user_prompt"] = repair(code)
        try:
            return decode(parse_object(model.complete(current)["content"]))
        except MendelPrimeError as error:
            code = error.code if error.category in {"model", "response"} else "stage_execution_failed"
        except Exception:
            code = "stage_execution_failed"
    raise MendelPrimeError(code, category="response", attempts=attempts)


def call_with_repair(
    model,
    case_key,
    stage,
    system,
    payload,
    max_tokens,
    decode,
    *,
    single_retry=False,
    max_input_characters=None,
):
    """One format repair; reports share that retry with transport failures."""
    request = {
        "case_key": case_key,
        "stage_id": stage,
        "system_prompt": system,
        "user_prompt": json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")),
        "max_tokens": max_tokens,
        "timeout_seconds": 300.0,
        "cache_choice": False,
        **({"max_http_attempts": 1} if single_retry else {}),
    }
    for attempt in range(2):
        if (
            max_input_characters is not None
            and len(request["system_prompt"]) + len(request["user_prompt"]) > max_input_characters
        ):
            raise MendelPrimeError("report_packet_budget_exceeded", category="contract")
        try:
            reply = model.complete(request)
        except MendelPrimeError as error:
            if single_retry and attempt == 0 and error.retryable:
                continue
            raise
        try:
            return decode(parse_object(reply["content"]))
        except MendelPrimeError as error:
            if attempt:
                raise
            code = error.code if error.category in {"model", "response"} else "assessment_record_invalid"
        except (KeyError, TypeError, ValueError) as error:
            if attempt:
                raise MendelPrimeError("assessment_record_invalid", category="response") from error
            code = "assessment_record_invalid"
        repair = {
            **payload,
            "format_repair": {
                "error": code,
                **({"draft": reply["content"]} if max_input_characters is None else {}),
                "instruction": (
                    "Correct only the identified schema, source, or evidence violation. Return "
                    "the complete JSON again; do not invent sources or force a desired class."
                ),
            },
        }
        request = {
            **request,
            "stage_id": stage + "_repair",
            "user_prompt": json.dumps(repair, ensure_ascii=False, sort_keys=True, separators=(",", ":")),
        }
