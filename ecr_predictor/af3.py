"""
AlphaFold 3 structure prediction — three backends.

  local  — run run_alphafold3.sh directly on the current machine
  hpcc   — upload input JSON to HPCC via SFTP, sbatch via SSH, poll, download CIF
  online — AlphaFold Server REST API (scaffold only, not implemented)

Select backend in config.yaml:  af3.backend: local | hpcc | online

SSH auth (hpcc backend): password-based via paramiko.
  Password is read from the ECR_HPCC_PASSWORD environment variable first,
  then falls back to af3.hpcc.ssh_password in config.yaml.
  Set it with: export ECR_HPCC_PASSWORD='yourpassword'

JSON format follows the AF3 input spec:
  https://github.com/google-deepmind/alphafold3#input-format
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Generator

import pandas as pd

AF3_OUTPUT_DIR = Path("af3_outputs")


# ---------------------------------------------------------------------------
# JSON input construction
# ---------------------------------------------------------------------------

def _build_af3_json(job_name: str, protein_sequence: str, dna_sequence: str) -> dict:
    """
    Build an AF3 JSON input dict for a DBD–DNA complex.

    Chain A = protein (DBD), chain B = DNA (single-stranded query sequence).
    modelSeeds uses [1,2,3,4,5] for 5 models, matching AF3 defaults.
    """
    return {
        "name": job_name,
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": protein_sequence,
                }
            },
            {
                "dna": {
                    "id": ["B"],
                    "sequence": dna_sequence,
                }
            },
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1,
    }


def _write_job_jsons(
    hits: pd.DataFrame,
    dna_sequence: str,
    job_dir: Path,
    top_n: int,
) -> list[tuple[str, Path]]:
    """
    Write AF3 JSON input files for the top-N hits.

    Returns list of (job_name, json_path) for hits that have a protein sequence.
    """
    job_dir.mkdir(parents=True, exist_ok=True)

    # Rank: FIMO-validated first, then fimo_pvalue, then motif_score
    ranked = hits.copy()
    if "fimo_validated" in ranked.columns:
        ranked["_fv"] = ranked["fimo_validated"].map(
            lambda x: 0 if str(x).lower() == "true" else 1
        )
        ranked["_fp"] = pd.to_numeric(ranked.get("fimo_pvalue", pd.NA), errors="coerce")
        ranked["_ms"] = pd.to_numeric(ranked.get("motif_score", pd.NA), errors="coerce")
        ranked = ranked.sort_values(["_fv", "_fp", "_ms"], ascending=[True, True, False])
        ranked = ranked.drop(columns=["_fv", "_fp", "_ms"])

    jobs = []
    for _, row in ranked.head(top_n).iterrows():
        gene = str(row["gene_name"])
        protein_seq = row.get("sequence_aa", "")
        if not protein_seq or pd.isna(protein_seq):
            print(f"  SKIP {gene}: no sequence_aa (re-run cli.py with --include-sequence).", file=sys.stderr)
            continue

        job_name = gene.replace(" ", "_")
        job_input = _build_af3_json(job_name, str(protein_seq), dna_sequence)
        json_path = job_dir / f"{job_name}.json"
        json_path.write_text(json.dumps(job_input, indent=2), encoding="utf-8")
        jobs.append((job_name, json_path))
        print(f"  Wrote {json_path.name}", file=sys.stderr)

    return jobs


# ---------------------------------------------------------------------------
# Local backend
# ---------------------------------------------------------------------------

def _run_local(
    jobs: list[tuple[str, Path]],
    output_dir: Path,
    af3_cfg: dict[str, Any],
) -> dict[str, Path | None]:
    """
    Run AF3 locally via run_alphafold3.sh for each job.

    If af3.local.module_load is set, the command is wrapped in:
      bash -c "source /etc/profile.d/modules.sh && module load <X> && run_alphafold3.sh ..."
    This handles environments (e.g. HPCC interactive sessions) where
    'module' is a shell function not available to plain subprocess calls.

    Returns {job_name: cif_path | None}
    """
    local_cfg = af3_cfg.get("local", {})
    af3_script = local_cfg.get("af3_script", "run_alphafold3.sh")
    module_load = local_cfg.get("module_load", "").strip()
    results: dict[str, Path | None] = {}

    for job_name, json_path in jobs:
        job_out = output_dir / job_name
        job_out.mkdir(parents=True, exist_ok=True)

        af3_call = (
            f'{af3_script} "{json_path.parent}" "{job_out}" "{json_path.name}"'
        )
        if module_load:
            # Source the module system then load the module before calling AF3.
            # /etc/profile.d/modules.sh is the standard location on most HPC systems;
            # adjust if your cluster uses a different path.
            bash_cmd = (
                f'source /etc/profile.d/modules.sh 2>/dev/null || true && '
                f'module load {module_load} && '
                f'{af3_call}'
            )
            cmd = ["bash", "-c", bash_cmd]
        else:
            cmd = ["bash", "-c", af3_call]

        print(f"  [local] Running AF3 for {job_name}...", file=sys.stderr)
        result = subprocess.run(cmd, capture_output=False, text=True)
        if result.returncode != 0:
            print(f"  WARNING: AF3 failed for {job_name} (exit {result.returncode})", file=sys.stderr)
            results[job_name] = None
            continue

        cif = _find_best_cif(job_out / job_name)
        results[job_name] = cif

    return results


# ---------------------------------------------------------------------------
# HPCC backend (paramiko SSH + SFTP + Slurm)
# ---------------------------------------------------------------------------

_SLURM_TEMPLATE = """\
#!/bin/bash
#SBATCH -N {nodes}
#SBATCH -c {cpus}
#SBATCH -J af3_{job_name}
#SBATCH --gres=gpu:{gpus}
#SBATCH -p {partition}
#SBATCH -q {qos}
#SBATCH -o {remote_workdir}/af3_{job_name}_%J.log
#SBATCH --mem={mem}

module load {module}

input_dir={remote_workdir}
output_dir={remote_workdir}

run_alphafold3.sh "$input_dir" "$output_dir" "{json_filename}"
"""


def _get_password(hpcc: dict[str, Any]) -> str:
    """
    Resolve the SSH password.
    Priority: ECR_HPCC_PASSWORD env var → config ssh_password field.
    """
    pw = os.environ.get("ECR_HPCC_PASSWORD", "")
    if pw:
        return pw
    pw = hpcc.get("ssh_password", "")
    if pw:
        return pw
    raise ValueError(
        "HPCC password not found. Set the ECR_HPCC_PASSWORD environment variable:\n"
        "  export ECR_HPCC_PASSWORD='yourpassword'\n"
        "Or set af3.hpcc.ssh_password in config.yaml."
    )


@contextmanager
def _ssh_client(hpcc: dict[str, Any]) -> Generator:
    """
    Open a paramiko SSHClient connection to the HPCC and yield it.
    Closes the connection on exit.
    """
    try:
        import paramiko
    except ImportError:
        raise ImportError("paramiko is required: pip install paramiko")

    host = hpcc["host"]
    user = hpcc["user"]
    port = int(hpcc.get("port", 22))
    password = _get_password(hpcc)

    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(host, port=port, username=user, password=password, timeout=15)
    try:
        yield client
    finally:
        client.close()


def _run_cmd(client, cmd: str) -> tuple[str, str, int]:
    """Run a command over an open paramiko SSHClient. Returns (stdout, stderr, exit_code)."""
    _, stdout, stderr = client.exec_command(cmd)
    exit_code = stdout.channel.recv_exit_status()
    return stdout.read().decode(), stderr.read().decode(), exit_code


def _sftp_put(client, local: Path, remote: str) -> None:
    """Upload a single file via SFTP."""
    with client.open_sftp() as sftp:
        sftp.put(str(local), remote)


def _sftp_get_dir(client, remote: str, local: Path) -> None:
    """
    Recursively download a remote directory via SFTP.
    Creates local directory structure as needed.
    """
    local.mkdir(parents=True, exist_ok=True)
    with client.open_sftp() as sftp:
        _sftp_get_recursive(sftp, remote, local)


def _sftp_get_recursive(sftp, remote: str, local: Path) -> None:
    import stat as stat_module
    for entry in sftp.listdir_attr(remote):
        remote_path = f"{remote}/{entry.filename}"
        local_path = local / entry.filename
        if stat_module.S_ISDIR(entry.st_mode):
            local_path.mkdir(exist_ok=True)
            _sftp_get_recursive(sftp, remote_path, local_path)
        else:
            sftp.get(remote_path, str(local_path))


def _submit_slurm_job(
    job_name: str,
    json_path: Path,
    hpcc: dict[str, Any],
    client,
) -> str | None:
    """
    Upload JSON + slurm script to HPCC, submit with sbatch.
    Returns the Slurm job ID string, or None on failure.
    """
    remote_workdir = hpcc["remote_workdir"]

    # Ensure remote workdir exists
    _run_cmd(client, f"mkdir -p {remote_workdir}")

    # Upload JSON
    _sftp_put(client, json_path, f"{remote_workdir}/{json_path.name}")

    # Write and upload slurm script
    slurm_content = _SLURM_TEMPLATE.format(
        nodes=hpcc.get("slurm_nodes", 1),
        cpus=hpcc.get("slurm_cpus", 8),
        gpus=hpcc.get("slurm_gpus", 1),
        partition=hpcc.get("slurm_partition", "a40-tmp"),
        qos=hpcc.get("slurm_qos", "gpu"),
        mem=hpcc.get("slurm_mem", "80G"),
        module=hpcc.get("slurm_module", "alphafold/3_a40-tmp"),
        remote_workdir=remote_workdir,
        job_name=job_name,
        json_filename=json_path.name,
    )
    slurm_local = json_path.parent / f"{job_name}.slurm"
    slurm_local.write_text(slurm_content, encoding="utf-8")
    _sftp_put(client, slurm_local, f"{remote_workdir}/{slurm_local.name}")

    # Submit
    stdout, stderr, rc = _run_cmd(client, f"cd {remote_workdir} && sbatch {slurm_local.name}")
    if rc != 0:
        print(f"  WARNING: sbatch failed for {job_name}:\n{stderr}", file=sys.stderr)
        return None

    # Parse job ID from "Submitted batch job 12345"
    for token in stdout.split():
        if token.isdigit():
            print(f"  Submitted {job_name} → Slurm job {token}", file=sys.stderr)
            return token

    print(f"  WARNING: could not parse job ID from: {stdout!r}", file=sys.stderr)
    return None


def _poll_slurm_job(
    job_id: str,
    job_name: str,
    hpcc: dict[str, Any],
    client,
) -> bool:
    """
    Poll squeue until the job is no longer listed, then check sacct.
    Returns True if COMPLETED, False on timeout or failure.
    """
    poll_interval = int(hpcc.get("poll_interval", 60))
    timeout = int(hpcc.get("timeout", 7200))
    elapsed = 0

    print(f"  Polling job {job_id} every {poll_interval}s (timeout {timeout}s)...", file=sys.stderr)
    while True:
        stdout, _, _ = _run_cmd(client, f"squeue -j {job_id} -h -o '%T' 2>/dev/null")
        state = stdout.strip()

        if not state:
            sacct_out, _, _ = _run_cmd(client, f"sacct -j {job_id} -n -o State -X 2>/dev/null")
            final = sacct_out.strip().split("\n")[0].strip() if sacct_out.strip() else "UNKNOWN"
            if "COMPLETED" in final:
                print(f"  Job {job_id} ({job_name}) COMPLETED.", file=sys.stderr)
                return True
            else:
                print(f"  Job {job_id} ({job_name}) ended with state: {final}", file=sys.stderr)
                return False

        print(f"  Job {job_id} state: {state} ({elapsed}s elapsed)", file=sys.stderr)

        if timeout and elapsed >= timeout:
            print(f"  Timeout waiting for job {job_id}.", file=sys.stderr)
            return False

        time.sleep(poll_interval)
        elapsed += poll_interval


# ---------------------------------------------------------------------------
# Run log — persists job state across interruptions
# ---------------------------------------------------------------------------

def _log_path(job_dir: Path) -> Path:
    return job_dir / "run_log.json"


def _load_log(job_dir: Path) -> dict[str, Any]:
    path = _log_path(job_dir)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return {}


def _save_log(job_dir: Path, log: dict[str, Any]) -> None:
    _log_path(job_dir).write_text(
        json.dumps(log, indent=2), encoding="utf-8"
    )


def _update_log(job_dir: Path, job_name: str, **fields) -> dict[str, Any]:
    """Update a single job's entry in the log and persist it."""
    log = _load_log(job_dir)
    entry = log.setdefault(job_name, {})
    entry.update(fields)
    _save_log(job_dir, log)
    return log


# ---------------------------------------------------------------------------
# HPCC orchestration with resume support
# ---------------------------------------------------------------------------

def _run_hpcc(
    jobs: list[tuple[str, Path]],
    output_dir: Path,
    af3_cfg: dict[str, Any],
) -> dict[str, Path | None]:
    """
    Submit jobs to HPCC via paramiko, poll for completion, retrieve outputs.

    Persists state to af3_outputs/jobs/run_log.json after every transition.
    On re-run:
      - already completed + CIF downloaded → returned immediately (no SSH)
      - job_id known but not yet complete   → re-polls (no re-submit)
      - no entry or failed/cancelled        → re-submits from scratch

    Returns {job_name: best_cif_path | None}
    """
    hpcc = af3_cfg.get("hpcc", {})
    if not hpcc.get("host") or not hpcc.get("user") or not hpcc.get("remote_workdir"):
        raise ValueError(
            "af3.hpcc.host, .user, and .remote_workdir must be set in config.yaml"
        )

    remote_workdir = hpcc["remote_workdir"]
    job_dir = jobs[0][1].parent  # af3_outputs/jobs/
    results: dict[str, Path | None] = {}

    # --- Fast path: return any already-downloaded CIFs without opening SSH ---
    pending_jobs = []
    for job_name, json_path in jobs:
        log = _load_log(job_dir)
        entry = log.get(job_name, {})
        if entry.get("status") == "completed" and entry.get("local_cif"):
            cif = Path(entry["local_cif"])
            if cif.exists():
                print(f"  {job_name}: already downloaded ({cif}), skipping.", file=sys.stderr)
                results[job_name] = cif
                continue
        pending_jobs.append((job_name, json_path))

    if not pending_jobs:
        return results

    # --- Open one SSH connection for all remaining work ---
    with _ssh_client(hpcc) as client:

        # Submit or re-use existing job IDs
        job_ids: dict[str, str] = {}
        for job_name, json_path in pending_jobs:
            log = _load_log(job_dir)
            entry = log.get(job_name, {})
            existing_jid = entry.get("job_id")
            existing_status = entry.get("status", "")

            if existing_jid and existing_status in ("submitted", "running"):
                # Connection was interrupted mid-poll — re-attach to the job
                print(f"  {job_name}: resuming poll for existing job {existing_jid}.", file=sys.stderr)
                job_ids[job_name] = existing_jid
            elif existing_status in ("failed", "cancelled", ""):
                # Submit fresh
                jid = _submit_slurm_job(job_name, json_path, hpcc, client)
                if jid:
                    job_ids[job_name] = jid
                    _update_log(job_dir, job_name,
                                job_id=jid, status="submitted",
                                remote_workdir=remote_workdir,
                                submitted_at=_now())
                else:
                    _update_log(job_dir, job_name, status="failed")
                    results[job_name] = None
            # status == "completed" without a local CIF: fall through to poll
            elif existing_jid:
                job_ids[job_name] = existing_jid

        # Poll each job to completion
        for job_name, jid in job_ids.items():
            _update_log(job_dir, job_name, status="running")
            ok = _poll_slurm_job(jid, job_name, hpcc, client)
            final_status = "completed" if ok else "failed"
            _update_log(job_dir, job_name, status=final_status)

            if not ok:
                results[job_name] = None
                continue

            # Download output
            job_name_lower = job_name.lower()
            local_out = output_dir / job_name_lower
            remote_out = f"{remote_workdir}/{job_name_lower}"
            try:
                _sftp_get_dir(client, remote_out, local_out)
            except Exception as e:
                print(f"  WARNING: SFTP download failed for {job_name}: {e}", file=sys.stderr)
                _update_log(job_dir, job_name, status="download_failed")
                results[job_name] = None
                continue

            cif = _find_best_cif(local_out)
            _update_log(job_dir, job_name,
                        local_cif=str(cif) if cif else None)
            results[job_name] = cif

    return results


def _now() -> str:
    """Return current UTC time as ISO string."""
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------------------
# Online backend — Chai-1 API (https://chaidiscovery.com)
# ---------------------------------------------------------------------------
#
# Chai-1 is an AF3-class model that supports protein–DNA complexes.
# Sign up at https://chaidiscovery.com to get an API key, then set:
#
#   export ECR_CHAI_API_KEY='your_key'
#   # or set af3.online.api_key in config.yaml
#
# API reference: https://chaidiscovery.com/docs/api
# ---------------------------------------------------------------------------

_CHAI_API_BASE = "https://api.chaidiscovery.com/v1"
_CHAI_POLL_INTERVAL = 30   # seconds between status polls
_CHAI_TIMEOUT = 3600       # max wait per job (seconds)


def _chai_api_key(online_cfg: dict[str, Any]) -> str:
    key = os.environ.get("ECR_CHAI_API_KEY", "")
    if key:
        return key
    key = online_cfg.get("api_key", "")
    if key:
        return key
    raise ValueError(
        "Chai-1 API key not found.\n"
        "Sign up at https://chaidiscovery.com, then:\n"
        "  export ECR_CHAI_API_KEY='your_key'\n"
        "Or set af3.online.api_key in config.yaml."
    )


def _chai_fasta(job_name: str, protein_sequence: str, dna_sequence: str) -> str:
    """
    Build a FASTA string for Chai-1: protein chain A + DNA chain B.
    Chai-1 uses sequence type tags in the FASTA header.
    """
    return (
        f">protein|name={job_name}_A\n{protein_sequence}\n"
        f">dna|name={job_name}_B\n{dna_sequence}\n"
    )


def _chai_submit(fasta: str, api_key: str) -> str | None:
    """Submit a prediction to the Chai-1 API. Returns job_id or None."""
    import requests
    resp = requests.post(
        f"{_CHAI_API_BASE}/predictions",
        headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
        json={"fasta": fasta, "use_msa_server": True},
        timeout=30,
    )
    if resp.status_code not in (200, 201, 202):
        print(f"  WARNING: Chai-1 submit failed ({resp.status_code}): {resp.text}", file=sys.stderr)
        return None
    return resp.json().get("id") or resp.json().get("job_id")


def _chai_poll(job_id: str, api_key: str) -> str | None:
    """
    Poll until the Chai-1 job finishes.
    Returns the CIF download URL, or None on failure/timeout.
    """
    import requests
    elapsed = 0
    while True:
        resp = requests.get(
            f"{_CHAI_API_BASE}/predictions/{job_id}",
            headers={"Authorization": f"Bearer {api_key}"},
            timeout=15,
        )
        if resp.status_code != 200:
            print(f"  WARNING: Chai-1 poll error ({resp.status_code}): {resp.text}", file=sys.stderr)
            return None

        data = resp.json()
        status = data.get("status", "unknown")
        print(f"  Chai-1 job {job_id}: {status} ({elapsed}s)", file=sys.stderr)

        if status in ("completed", "success"):
            # The CIF URL may be nested; try common key names
            return (
                data.get("cif_url")
                or data.get("result", {}).get("cif_url")
                or data.get("download_url")
            )
        if status in ("failed", "error", "cancelled"):
            print(f"  Chai-1 job {job_id} ended with status: {status}", file=sys.stderr)
            return None

        if elapsed >= _CHAI_TIMEOUT:
            print(f"  Chai-1 job {job_id} timed out after {elapsed}s.", file=sys.stderr)
            return None

        time.sleep(_CHAI_POLL_INTERVAL)
        elapsed += _CHAI_POLL_INTERVAL


def _chai_download(cif_url: str, dest: Path, api_key: str) -> Path | None:
    """Download the CIF file from Chai-1 to dest."""
    import requests
    dest.parent.mkdir(parents=True, exist_ok=True)
    resp = requests.get(
        cif_url,
        headers={"Authorization": f"Bearer {api_key}"},
        timeout=60,
    )
    if resp.status_code != 200:
        print(f"  WARNING: CIF download failed ({resp.status_code})", file=sys.stderr)
        return None
    dest.write_bytes(resp.content)
    return dest


def _run_online(
    jobs: list[tuple[str, Path]],
    output_dir: Path,
    af3_cfg: dict[str, Any],
) -> dict[str, Path | None]:
    """
    Submit protein–DNA structure predictions to the Chai-1 API.

    Requires an API key from https://chaidiscovery.com.
    Set ECR_CHAI_API_KEY env var or af3.online.api_key in config.yaml.
    """
    online_cfg = af3_cfg.get("online", {})
    api_key = _chai_api_key(online_cfg)
    results: dict[str, Path | None] = {}

    for job_name, json_path in jobs:
        # Read protein/DNA sequences back from the already-written JSON
        job_input = json.loads(json_path.read_text(encoding="utf-8"))
        seqs = job_input.get("sequences", [])
        protein_seq = next(
            (s["protein"]["sequence"] for s in seqs if "protein" in s), None
        )
        dna_seq = next(
            (s["dna"]["sequence"] for s in seqs if "dna" in s), None
        )
        if not protein_seq or not dna_seq:
            print(f"  SKIP {job_name}: could not extract sequences from JSON.", file=sys.stderr)
            results[job_name] = None
            continue

        fasta = _chai_fasta(job_name, protein_seq, dna_seq)

        print(f"  [online] Submitting {job_name} to Chai-1...", file=sys.stderr)
        job_id = _chai_submit(fasta, api_key)
        if not job_id:
            results[job_name] = None
            continue
        print(f"  Chai-1 job ID: {job_id}", file=sys.stderr)

        cif_url = _chai_poll(job_id, api_key)
        if not cif_url:
            results[job_name] = None
            continue

        cif_path = output_dir / job_name / f"{job_name}_chai.cif"
        downloaded = _chai_download(cif_url, cif_path, api_key)
        results[job_name] = downloaded
        if downloaded:
            print(f"  Saved: {downloaded}", file=sys.stderr)

    return results


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_best_cif(job_output_dir: Path) -> Path | None:
    """
    Return the best CIF from an AF3 output directory.

    AF3 lowercases the job name, so for job "FLI1" the output is:
      fli1/fli1_model.cif
    We prefer *_model.cif; fall back to any .cif if not found.
    """
    if not job_output_dir.exists():
        return None

    top = list(job_output_dir.glob("*_model.cif"))
    if top:
        return sorted(top)[0]

    all_cifs = sorted(job_output_dir.glob("*.cif"))
    return all_cifs[0] if all_cifs else None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_af3_prediction(
    hits: pd.DataFrame,
    dna_sequence: str,
    af3_cfg: dict[str, Any] | None = None,
    output_dir: Path = AF3_OUTPUT_DIR,
    top_n: int = 5,
) -> pd.DataFrame:
    """
    Predict DBD–DNA complex structures with AlphaFold 3.

    Parameters
    ----------
    hits       : FIMO-validated predictor output; needs 'gene_name', 'sequence_aa'
    dna_sequence : original DNA query sequence
    af3_cfg    : the af3 section from config.yaml (from ecr_predictor.config)
    output_dir : local root for AF3 inputs and downloaded outputs
    top_n      : number of top hits to submit

    Returns
    -------
    hits with new column 'af3_cif_path' (path string or NA).
    """
    if af3_cfg is None:
        af3_cfg = {"backend": "hpcc"}

    backend = af3_cfg.get("backend", "hpcc")
    output_dir = Path(output_dir)
    job_dir = output_dir / "jobs"

    if "sequence_aa" not in hits.columns or hits["sequence_aa"].isna().all():
        print(
            "WARNING: 'sequence_aa' column missing or all NA. "
            "Re-run cli.py with --include-sequence.",
            file=sys.stderr,
        )
        hits = hits.copy()
        hits["af3_cif_path"] = pd.NA
        return hits

    print(f"[AF3] Backend: {backend}", file=sys.stderr)
    jobs = _write_job_jsons(hits, dna_sequence, job_dir, top_n)
    if not jobs:
        print("[AF3] No jobs to submit (all hits missing sequence_aa).", file=sys.stderr)
        hits = hits.copy()
        hits["af3_cif_path"] = pd.NA
        return hits

    if backend == "local":
        results = _run_local(jobs, output_dir, af3_cfg)
    elif backend == "hpcc":
        results = _run_hpcc(jobs, output_dir, af3_cfg)
    elif backend == "online":
        results = _run_online(jobs, output_dir, af3_cfg)
    else:
        raise ValueError(f"Unknown AF3 backend: {backend!r}. Choose local | hpcc | online.")

    hits = hits.copy()
    hits["af3_cif_path"] = hits["gene_name"].apply(
        lambda g: str(results[g.replace(" ", "_")]) if g.replace(" ", "_") in results and results[g.replace(" ", "_")] is not None else pd.NA
    )

    n_done = hits["af3_cif_path"].notna().sum()
    print(f"[AF3] {n_done}/{len(jobs)} structures retrieved.", file=sys.stderr)
    return hits
