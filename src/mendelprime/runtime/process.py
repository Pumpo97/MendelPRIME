"""Private external-process execution with process-group cleanup."""

import os
import signal
import subprocess
import time
from pathlib import Path

_TIMEOUT_RETURN_CODE = 124
_TERMINATION_GRACE_SECONDS = 0.5


def batch_outcome(target, status, error_code="", annotation=None):
    return {"target": target, "annotation": annotation, "status": status, "error_code": error_code}


def fail_batch(outcomes, eligible, error_code):
    for index, target, *_ in eligible:
        outcomes[index] = batch_outcome(target, "failed", error_code)


def complete_batch(targets, outcomes, missing_code):
    return tuple(
        outcome or batch_outcome(target, "failed", missing_code)
        for target, outcome in zip(targets, outcomes, strict=True)
    )


def _terminate_process_tree(process):
    if os.name == "nt":
        try:
            subprocess.run(
                ["taskkill", "/PID", str(process.pid), "/T", "/F"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
                shell=False,
            )
        except OSError:
            process.kill()
        try:
            process.wait(timeout=_TERMINATION_GRACE_SECONDS)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
        return

    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        process.wait()
        return
    deadline = time.monotonic() + _TERMINATION_GRACE_SECONDS
    while time.monotonic() < deadline:
        process.poll()
        try:
            os.killpg(process.pid, 0)
        except ProcessLookupError:
            process.wait()
            return
        time.sleep(min(0.02, max(0.0, deadline - time.monotonic())))
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def run_process(*, command, root: Path, timeout_seconds, artifacts=()):
    """Run one argv in a private directory and return a sanitized receipt."""
    started = time.monotonic()
    returncode = None
    error_code = ""
    stdout_path, stderr_path = root / "stdout.log", root / "stderr.log"
    creationflags = subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0
    try:
        with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
            process = subprocess.Popen(
                command,
                cwd=root,
                stdout=stdout,
                stderr=stderr,
                shell=False,
                start_new_session=os.name != "nt",
                creationflags=creationflags,
            )
            try:
                returncode = process.wait(timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                _terminate_process_tree(process)
                returncode = _TIMEOUT_RETURN_CODE
                error_code = "tool_timeout"
            except BaseException:
                _terminate_process_tree(process)
                raise
    except OSError:
        error_code = "tool_start_failed"

    if not error_code and returncode != 0:
        error_code = "tool_nonzero_exit"
    if not error_code:
        if any(not artifact.is_file() for artifact in artifacts):
            error_code = "tool_artifact_missing"
        elif any(artifact.stat().st_size == 0 for artifact in artifacts):
            error_code = "tool_artifact_empty"
    return {
        "returncode": returncode,
        "error_code": error_code,
        "seconds": round(time.monotonic() - started, 6),
    }
