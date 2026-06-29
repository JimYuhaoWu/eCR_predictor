"""
Shared local-CLI-or-API backend dispatch for fusion-stage tools.

Every external tool (NetMHCpan, NetCTLpan, AGGRESCAN3D, CamSol, …) is configured
under fusion.tools.<tool> in config.yaml with a `backend` key:

  backend: local       run a local command-line binary
  backend: api         submit to a remote HTTP API (submit → poll → fetch)
  backend: disabled    skip this tool (gate degrades gracefully)

This mirrors the AF3 stage's backend selection. Gate modules build the
tool-specific command / payload and parse the tool-specific output; this module
only handles the local-vs-api plumbing common to all of them.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Callable


class ToolDisabled(RuntimeError):
    """Raised when a tool's backend is set to 'disabled'."""


class ToolNotAvailableError(RuntimeError):
    """Raised when a local binary is missing or an API is unreachable."""


def resolve_backend(tool_cfg: dict[str, Any]) -> str:
    """Return 'local' | 'api' | 'disabled' for a tool config block."""
    backend = str(tool_cfg.get("backend", "disabled")).lower()
    if backend not in ("local", "api", "disabled"):
        raise ValueError(
            f"Invalid tool backend {backend!r}. Choose local | api | disabled."
        )
    return backend


def api_key_from_env(api_cfg: dict[str, Any]) -> str:
    """
    Resolve an API key from the env var named in api_cfg['api_key_env'],
    falling back to api_cfg['api_key']. Returns '' if neither is set.
    """
    env_name = api_cfg.get("api_key_env", "")
    if env_name:
        key = os.environ.get(env_name, "")
        if key:
            return key
    return api_cfg.get("api_key", "")


# ---------------------------------------------------------------------------
# Local CLI invocation
# ---------------------------------------------------------------------------

def run_cli(
    command: str,
    args: list[str],
    input_files: dict[str, str] | None = None,
    stdin_text: str | None = None,
    timeout: int | None = None,
) -> str:
    """
    Run a local command and return its stdout.

    Parameters
    ----------
    command : binary name or path (resolved on PATH)
    args : argument list; the literal token '{<key>}' is replaced by the temp
           path of input_files[<key>] so commands can reference written inputs.
    input_files : {key: text} written to temp files, exposed as '{key}' in args
    stdin_text : optional text piped to the process stdin
    timeout : seconds (None = no limit)

    Raises ToolNotAvailableError if the binary is missing.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        subst: dict[str, str] = {}
        for key, text in (input_files or {}).items():
            p = tmpdir / key
            p.write_text(text, encoding="utf-8")
            subst[key] = str(p)

        resolved_args = [a.format(**subst) if "{" in a else a for a in args]
        cmd = [command, *resolved_args]

        try:
            result = subprocess.run(
                cmd,
                input=stdin_text,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except FileNotFoundError:
            raise ToolNotAvailableError(
                f"Local binary not found: {command!r}. "
                f"Install it or set this tool's backend to 'api' / 'disabled'."
            )

        if result.returncode != 0:
            print(
                f"  WARNING: {command} exited {result.returncode}:\n"
                f"{result.stderr[-500:]}",
                file=sys.stderr,
            )
        return result.stdout


# ---------------------------------------------------------------------------
# Remote API invocation (submit → poll → fetch)
# ---------------------------------------------------------------------------

def submit_and_poll(
    api_cfg: dict[str, Any],
    payload: dict[str, Any],
    parse_job_id: Callable[[dict], str | None],
    parse_status: Callable[[dict], str],
    parse_result: Callable[[dict], Any],
    done_states: tuple[str, ...] = ("completed", "success"),
    fail_states: tuple[str, ...] = ("failed", "error", "cancelled"),
) -> Any:
    """
    Generic submit/poll loop for a remote tool API, mirroring the Chai-1 client.

    The caller supplies small callbacks to extract the job id, status, and final
    result from each JSON response, keeping this loop tool-agnostic.

    Returns parse_result(final_response), or None on failure/timeout.
    """
    import requests

    url = api_cfg.get("url", "")
    if not url:
        raise ToolNotAvailableError("API backend selected but fusion tool 'url' is empty.")
    key = api_key_from_env(api_cfg)
    headers = {"Content-Type": "application/json"}
    if key:
        headers["Authorization"] = f"Bearer {key}"

    poll_interval = int(api_cfg.get("poll_interval", 10))
    timeout = int(api_cfg.get("timeout", 600))

    resp = requests.post(url, headers=headers, json=payload, timeout=30)
    if resp.status_code not in (200, 201, 202):
        print(f"  WARNING: API submit failed ({resp.status_code}): {resp.text[:300]}",
              file=sys.stderr)
        return None

    data = resp.json()
    job_id = parse_job_id(data)
    if job_id is None:
        # Some APIs return the result synchronously on submit.
        status = parse_status(data)
        if status in done_states:
            return parse_result(data)
        return None

    elapsed = 0
    while True:
        r = requests.get(f"{url}/{job_id}", headers=headers, timeout=15)
        if r.status_code != 200:
            print(f"  WARNING: API poll error ({r.status_code}): {r.text[:300]}",
                  file=sys.stderr)
            return None
        data = r.json()
        status = parse_status(data)
        print(f"  API job {job_id}: {status} ({elapsed}s)", file=sys.stderr)
        if status in done_states:
            return parse_result(data)
        if status in fail_states:
            print(f"  API job {job_id} ended with status: {status}", file=sys.stderr)
            return None
        if elapsed >= timeout:
            print(f"  API job {job_id} timed out after {elapsed}s.", file=sys.stderr)
            return None
        time.sleep(poll_interval)
        elapsed += poll_interval
