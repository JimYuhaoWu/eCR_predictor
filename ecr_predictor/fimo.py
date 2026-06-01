"""
FIMO-based validation of motif hits.

Workflow:
  1. For each hit in the filtered table, load its JASPAR motif (from cache or API).
  2. Convert all motifs to a single MEME-format file.
  3. Write the query sequence as a FASTA file.
  4. Call `fimo --text` and parse the output TSV.
  5. Annotate the hit table with FIMO p-values.

Prerequisites:
  - MEME Suite installed and `fimo` on PATH.
  - JASPAR cache populated (run `python -m ecr_predictor.prefetch` first).

If `fimo` is not found, FIMONotAvailableError is raised with install instructions.
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from io import StringIO
from pathlib import Path

import pandas as pd

from ecr_predictor.scan import _fetch_all_parallel

# p-value threshold for a FIMO hit to be considered validated
FIMO_PVALUE_THRESHOLD = 1e-4


class FIMONotAvailableError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# MEME-format conversion
# ---------------------------------------------------------------------------

def _jaspar_to_meme_block(jaspar_id: str, motif) -> str:
    """
    Convert a BioPython JASPAR Motif object to a MEME minimal-format block.

    MEME format reference:
      https://meme-suite.org/meme/doc/meme-format.html
    """
    length = len(motif)
    nsites = sum(motif.counts["A"]) if "A" in motif.counts else 20

    lines = [
        f"MOTIF {jaspar_id} {getattr(motif, 'name', jaspar_id)}",
        f"letter-probability matrix: alength= 4 w= {length} nsites= {int(nsites)} E= 0",
    ]
    pwm = motif.counts.normalize(pseudocounts=0.5)
    for i in range(length):
        row = "\t".join(
            f"{pwm[nt][i]:.6f}" for nt in ("A", "C", "G", "T")
        )
        lines.append(row)
    lines.append("")
    return "\n".join(lines)


def _build_meme_file(motif_dict: dict[str, object]) -> str:
    """
    Build a complete MEME-format string from a {jaspar_id: motif} dict.
    Skips IDs with None motifs.
    """
    header = (
        "MEME version 4\n\n"
        "ALPHABET= ACGT\n\n"
        "strands: + -\n\n"
        "Background letter frequencies\n"
        "A 0.25 C 0.25 G 0.25 T 0.25\n\n"
    )
    blocks = []
    for jid, motif in motif_dict.items():
        if motif is None:
            continue
        blocks.append(_jaspar_to_meme_block(jid, motif))
    return header + "\n".join(blocks)


# ---------------------------------------------------------------------------
# FIMO invocation
# ---------------------------------------------------------------------------

def _check_fimo() -> None:
    """Raise FIMONotAvailableError if `fimo` is not on PATH."""
    result = subprocess.run(
        ["which", "fimo"], capture_output=True, text=True
    )
    if result.returncode != 0:
        raise FIMONotAvailableError(
            "fimo not found on PATH.\n"
            "Install MEME Suite: https://meme-suite.org/meme/doc/install.html\n"
            "  conda install -c bioconda meme   # easiest route"
        )


def _run_fimo(meme_text: str, sequence: str, pvalue_thresh: float) -> pd.DataFrame:
    """
    Write temp files, call fimo, parse results into a DataFrame.

    Returns columns: motif_id, sequence_name, start, stop, strand, score, p-value, q-value, matched_sequence
    Returns empty DataFrame if no hits.
    """
    _check_fimo()

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        meme_path = tmpdir / "motifs.meme"
        fasta_path = tmpdir / "query.fasta"
        out_dir = tmpdir / "fimo_out"

        meme_path.write_text(meme_text, encoding="utf-8")
        fasta_path.write_text(f">query\n{sequence}\n", encoding="utf-8")

        cmd = [
            "fimo",
            "--text",                        # TSV output to stdout, no HTML
            "--thresh", str(pvalue_thresh),
            "--parse-genomic-coord",         # report coords relative to FASTA header
            str(meme_path),
            str(fasta_path),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"WARNING: fimo exited with code {result.returncode}", file=sys.stderr)
            print(result.stderr, file=sys.stderr)

        stdout = result.stdout.strip()
        if not stdout:
            return pd.DataFrame()

        # fimo --text emits a header line starting with '#'
        lines = [l for l in stdout.splitlines() if not l.startswith("#")]
        if not lines:
            return pd.DataFrame()

        col_line = "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence"
        df = pd.read_csv(
            StringIO(col_line + "\n" + "\n".join(lines)),
            sep="\t",
        )
        return df


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_fimo_validation(
    hits: pd.DataFrame,
    sequence: str,
    pvalue_thresh: float = FIMO_PVALUE_THRESHOLD,
) -> pd.DataFrame:
    """
    Validate motif hits with FIMO.

    Parameters
    ----------
    hits : filtered predictor output (must have 'jaspar_id', 'gene_name' columns)
    sequence : the original DNA query sequence
    pvalue_thresh : FIMO p-value cutoff (default 1e-4)

    Returns
    -------
    hits with two new columns:
      fimo_pvalue  — best FIMO p-value for that motif (NaN if no FIMO hit)
      fimo_validated — True if fimo_pvalue <= pvalue_thresh
    """
    jaspar_ids = hits["jaspar_id"].dropna().unique().tolist()

    print(f"[FIMO] Fetching {len(jaspar_ids)} motifs for FIMO...", file=sys.stderr)
    motif_dict = _fetch_all_parallel(jaspar_ids)

    meme_text = _build_meme_file(motif_dict)
    if not any(m is not None for m in motif_dict.values()):
        print("WARNING: no motifs available for FIMO; skipping.", file=sys.stderr)
        hits = hits.copy()
        hits["fimo_pvalue"] = float("nan")
        hits["fimo_validated"] = False
        return hits

    print("[FIMO] Running fimo...", file=sys.stderr)
    fimo_results = _run_fimo(meme_text, sequence, pvalue_thresh)

    # Summarise: best (lowest) p-value per motif_id
    if fimo_results.empty:
        best_pval: dict[str, float] = {}
    else:
        best_pval = (
            fimo_results.groupby("motif_id")["p-value"].min().to_dict()
        )

    hits = hits.copy()
    hits["fimo_pvalue"] = hits["jaspar_id"].map(best_pval)
    hits["fimo_validated"] = hits["fimo_pvalue"].apply(
        lambda p: (not pd.isna(p)) and (p <= pvalue_thresh)
    )

    n_validated = hits["fimo_validated"].sum()
    print(
        f"[FIMO] {n_validated}/{len(hits)} hits FIMO-validated "
        f"(p ≤ {pvalue_thresh})",
        file=sys.stderr,
    )
    return hits
