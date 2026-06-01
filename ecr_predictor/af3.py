"""
AlphaFold 3 structure prediction stub.

TODO: implement one of the following invocation strategies:
  A) AlphaFold Server API (https://alphafoldserver.com) — POST job, poll, download CIF
  B) Local AF3 install — call `python run_alphafold.py` via subprocess
  C) HPC job submission — write JSON input, sbatch/qsub wrapper

Input contract (per hit):
  - protein_sequence : str  — DBD amino acid sequence (use --include-sequence in cli.py)
  - dna_sequence     : str  — the original query DNA sequence
  - gene_name        : str  — used to name output files

Output contract:
  - Returns a dict {gene_name: Path} mapping each hit to its predicted CIF/PDB file.
  - CIF files are written to `af3_outputs/<gene_name>/` relative to the working directory.

See: https://github.com/google-deepmind/alphafold3
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd


AF3_OUTPUT_DIR = Path("af3_outputs")


def _build_af3_json(gene_name: str, protein_sequence: str, dna_sequence: str) -> dict:
    """
    Build the AF3 JSON job input for a DBD–DNA complex.

    AF3 expects a list of sequences with type annotations:
      https://github.com/google-deepmind/alphafold3#input-format
    """
    return {
        "name": gene_name,
        "sequences": [
            {
                "proteinChain": {
                    "sequence": protein_sequence,
                    "count": 1,
                }
            },
            {
                "dnaSequence": {
                    "sequence": dna_sequence,
                    "count": 1,
                }
            },
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1,
    }


def run_af3_prediction(
    hits: pd.DataFrame,
    dna_sequence: str,
    output_dir: Path = AF3_OUTPUT_DIR,
    top_n: int = 5,
) -> pd.DataFrame:
    """
    Predict DBD–DNA complex structures with AlphaFold 3 for the top-N hits.

    Parameters
    ----------
    hits : FIMO-validated predictor output; must have 'gene_name' and 'sequence_aa'
    dna_sequence : original query DNA sequence
    output_dir : root directory for AF3 inputs/outputs
    top_n : number of top hits to submit (ranked by fimo_pvalue, then motif_score)

    Returns
    -------
    hits with new column 'af3_cif_path' (path string or NaN if not run).

    Raises
    ------
    NotImplementedError — until AF3 invocation is implemented.
    """
    if "sequence_aa" not in hits.columns or hits["sequence_aa"].isna().all():
        print(
            "WARNING: 'sequence_aa' column missing or empty. "
            "Re-run cli.py with --include-sequence to enable AF3 prediction.",
            file=sys.stderr,
        )
        hits = hits.copy()
        hits["af3_cif_path"] = pd.NA
        return hits

    # Rank: FIMO-validated first, then by fimo_pvalue, then motif_score
    ranked = hits.copy()
    if "fimo_validated" in ranked.columns:
        ranked = ranked.sort_values(
            ["fimo_validated", "fimo_pvalue", "motif_score"],
            ascending=[False, True, False],
        )
    top = ranked.head(top_n)

    output_dir.mkdir(parents=True, exist_ok=True)
    job_dir = output_dir / "jobs"
    job_dir.mkdir(exist_ok=True)

    submitted = []
    for _, row in top.iterrows():
        gene = row["gene_name"]
        protein_seq = row.get("sequence_aa", "")
        if not protein_seq or pd.isna(protein_seq):
            print(f"  SKIP {gene}: no amino acid sequence.", file=sys.stderr)
            submitted.append((gene, pd.NA))
            continue

        job_input = _build_af3_json(gene, str(protein_seq), dna_sequence)
        json_path = job_dir / f"{gene}.json"
        json_path.write_text(json.dumps(job_input, indent=2), encoding="utf-8")
        print(f"  Wrote AF3 input: {json_path}", file=sys.stderr)
        submitted.append((gene, json_path))

    # TODO: implement actual AF3 submission/polling here.
    # Options:
    #   A) Upload json_path to AlphaFold Server via REST API
    #   B) subprocess.run(["python", "/path/to/run_alphafold.py", "--json_path", ...])
    #   C) Write an HPC submission script
    raise NotImplementedError(
        "AF3 invocation not yet implemented. "
        "JSON job inputs have been written to: " + str(job_dir) + "\n"
        "See ecr_predictor/af3.py for implementation options."
    )
