"""
Structure phase — runs only on survivors of the cheap sequence gates.

Folds each surviving fusion (DBD + linker + ED) with the DNA target via the
existing AF3 stage, then estimates the DBD–DNA interaction energy with FoldX.
This is the HPCC/GPU-intensive part, deliberately deferred until the candidate
set has been pruned by the sequence-only gates (immunogenicity, stability).

Function retention = does the DBD still bind its DNA after fusion. We report the
fusion's FoldX ΔΔG and, when a per-DBD baseline ΔΔG is available (from a prior
refine.py run, carried in the DBD input TSV), the delta vs that baseline —
the comparable, interpretable signal (how much fusing perturbed binding).

Reuses ecr_predictor.af3 and ecr_predictor.foldx unchanged; nothing AF3-specific
is reimplemented here.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from ecr_predictor.fusion.assemble import FusionCandidate


def run_structure_phase(
    candidates: list[FusionCandidate],
    dna_sequence: str,
    dbd_meta: dict[str, dict],
    af3_cfg: dict,
    output_dir: Path,
) -> dict[str, dict]:
    """
    Fold survivors and compute FoldX ΔΔG.

    Parameters
    ----------
    candidates : surviving fusion candidates
    dna_sequence : DNA target the DBD must still bind (same as Steps 1–2)
    dbd_meta : {dbd_name: {"tf_family":..., "zinc_finger_count":..., "baseline_ddg":...}}
               propagates zinc-finger handling and the function-retention baseline.
    af3_cfg : the af3 config section (same one refine.py uses)
    output_dir : AF3 output root for this run

    Returns {candidate.name: {"af3_cif_path", "fusion_ddg", "function_delta_ddg"}}.
    """
    from ecr_predictor.af3 import run_af3_prediction
    from ecr_predictor.foldx import run_foldx_affinity

    # Build a hits-like frame the AF3/FoldX stages understand. The fusion is the
    # "protein"; carry the DBD's tf_family/zinc_finger_count so zinc-finger DBDs
    # still get Zn ions in the AF3 input.
    rows = []
    for c in candidates:
        meta = dbd_meta.get(c.dbd_name, {})
        rows.append({
            "gene_name": c.name,
            "sequence_aa": c.sequence,
            "tf_family": meta.get("tf_family", ""),
            "zinc_finger_count": meta.get("zinc_finger_count"),
        })
    df = pd.DataFrame(rows)

    print(f"[Structure] Folding {len(df)} survivors via AF3...", file=sys.stderr)
    df = run_af3_prediction(df, dna_sequence, af3_cfg, output_dir, top_n=len(df))

    print("[Structure] Estimating FoldX ΔΔG (function retention)...", file=sys.stderr)
    df = run_foldx_affinity(df, work_base_dir=output_dir.parent / "foldx_work")

    out: dict[str, dict] = {}
    for c, (_, row) in zip(candidates, df.iterrows()):
        cif = row.get("af3_cif_path")
        ddg = pd.to_numeric(row.get("foldx_ddg_kcal_mol"), errors="coerce")
        baseline = dbd_meta.get(c.dbd_name, {}).get("baseline_ddg")
        delta = (float(ddg) - float(baseline)) if (pd.notna(ddg) and baseline is not None) else None
        out[c.name] = {
            "af3_cif_path": None if pd.isna(cif) else str(cif),
            "fusion_ddg": None if pd.isna(ddg) else float(ddg),
            "function_delta_ddg": delta,
        }
    return out
