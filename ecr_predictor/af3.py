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
        "modelSeeds": [1, 2, 3, 4, 5],
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
#SBATCH -o af3_{job_name}_%J.log
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


def _run_hpcc(
    jobs: list[tuple[str, Path]],
    output_dir: Path,
    af3_cfg: dict[str, Any],
) -> dict[str, Path | None]:
    """
    Submit all jobs to HPCC via paramiko, poll for completion, retrieve outputs.
    Returns {job_name: best_cif_path | None}
    """
    hpcc = af3_cfg.get("hpcc", {})
    if not hpcc.get("host") or not hpcc.get("user") or not hpcc.get("remote_workdir"):
        raise ValueError(
            "af3.hpcc.host, .user, and .remote_workdir must be set in config.yaml"
        )

    remote_workdir = hpcc["remote_workdir"]
    results: dict[str, Path | None] = {}

    with _ssh_client(hpcc) as client:
        # Submit all jobs
        job_ids: dict[str, str] = {}
        for job_name, json_path in jobs:
            jid = _submit_slurm_job(job_name, json_path, hpcc, client)
            if jid:
                job_ids[job_name] = jid

        # Poll all jobs to completion
        completed: set[str] = set()
        for job_name, jid in job_ids.items():
            ok = _poll_slurm_job(jid, job_name, hpcc, client)
            if ok:
                completed.add(job_name)

        # Retrieve outputs
        for job_name, _ in jobs:
            if job_name not in completed:
                results[job_name] = None
                continue

            local_out = output_dir / job_name
            remote_out = f"{remote_workdir}/{job_name}"
            try:
                _sftp_get_dir(client, remote_out, local_out)
            except Exception as e:
                print(f"  WARNING: SFTP download failed for {job_name}: {e}", file=sys.stderr)
                results[job_name] = None
                continue

            cif = _find_best_cif(local_out / job_name)
            results[job_name] = cif

    return results


# ---------------------------------------------------------------------------
# Online backend (stub)
# ---------------------------------------------------------------------------

def _run_online(
    jobs: list[tuple[str, Path]],
    output_dir: Path,
    af3_cfg: dict[str, Any],
) -> dict[str, Path | None]:
    raise NotImplementedError(
        "AlphaFold Server online backend is not yet implemented.\n"
        "Use backend: local or backend: hpcc in config.yaml."
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_best_cif(job_output_dir: Path) -> Path | None:
    """
    Return the best-ranked CIF from an AF3 output directory.

    AF3 names models like <job_name>_model.cif (single) or
    model_<seed>_<rank>.cif. We prefer a file named *_model.cif
    (the top-ranked output), falling back to any .cif file.
    """
    if not job_output_dir.exists():
        return None

    # Prefer *_model.cif (AF3's top-ranked output naming)
    top = list(job_output_dir.glob("*_model.cif"))
    if top:
        return top[0]

    # Fallback: any cif, sorted for determinism
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
        lambda g: str(results[g.replace(" ", "_")]) if g.replace(" ", "_") in results and results[g.replace(" ", "_")] else pd.NA
    )

    n_done = hits["af3_cif_path"].notna().sum()
    print(f"[AF3] {n_done}/{len(jobs)} structures retrieved.", file=sys.stderr)
    return hits
