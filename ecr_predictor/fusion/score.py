"""
Composite risk scoring and Pareto ranking of fusion candidates.

The gates produce independent liability axes. Rather than collapse them into one
opaque number (which would hide trade-offs), we report each axis and mark the
Pareto-optimal set — candidates not strictly worse than any other on all axes.
A convenience composite `risk_score` is also provided for a quick sort, but the
Pareto flag is the decision-relevant output.

Axes (all "lower is better"):
  immuno_density      fraction of junction peptides flagged as MHC-I binders
  aggregation_risk    1.0 if a junction aggregation hotspot, else 0.0
  stability_risk      count of stability liabilities (N-end + degrons + Ub-Lys)
"""
from __future__ import annotations

import pandas as pd


def build_score_table(
    candidates,
    immuno: dict | None = None,
    aggregation: dict | None = None,
    stability: dict | None = None,
    structure: dict | None = None,
) -> pd.DataFrame:
    """Assemble a per-candidate liability table from the gate result dicts.

    `structure` maps candidate.name → {"af3_cif_path", "fusion_ddg",
    "function_delta_ddg"} from the structure phase (survivors only)."""
    immuno = immuno or {}
    aggregation = aggregation or {}
    stability = stability or {}
    structure = structure or {}

    rows = []
    for cand in candidates:
        im = immuno.get(cand.name)
        ag = aggregation.get(cand.name)
        st = stability.get(cand.name)
        sb = structure.get(cand.name, {})

        stability_risk = 0
        if st is not None:
            stability_risk = (
                int(st.nterm_destabilizing)
                + len(st.degron_hits)
                + len(st.ubiquitination_lys)
            )

        rows.append({
            "candidate": cand.name,
            "dbd": cand.dbd_name,
            "linker": cand.linker_name,
            "ed": cand.ed_name,
            "length": len(cand.sequence),
            "immuno_density": im.density if im else pd.NA,
            "immuno_flagged": im.n_flagged if im else pd.NA,
            "immuno_min_rank": (im.min_rank if im and im.min_rank != float("inf") else pd.NA),
            "aggregation_risk": (1.0 if ag and ag.junction_hotspot else 0.0) if ag else pd.NA,
            "stability_risk": stability_risk if st else pd.NA,
            # Structure phase (survivors only); function_delta_ddg is the
            # comparable function-retention signal (reported, not a Pareto axis).
            "fusion_ddg_kcal_mol": sb.get("fusion_ddg", pd.NA),
            "function_delta_ddg": sb.get("function_delta_ddg", pd.NA),
            "af3_cif_path": sb.get("af3_cif_path", pd.NA),
        })
    return pd.DataFrame(rows)


# Full axis set. The sequence-only phase ranks on a subset (no structure yet).
ALL_AXES = ["immuno_density", "aggregation_risk", "stability_risk"]
SEQUENCE_AXES = ["immuno_density", "stability_risk"]


def _dominates(a: dict, b: dict, axes: list[str]) -> bool:
    """True if a is no worse than b on all axes and strictly better on one."""
    no_worse = all(a[x] <= b[x] for x in axes)
    strictly_better = any(a[x] < b[x] for x in axes)
    return no_worse and strictly_better


def add_pareto_and_rank(df: pd.DataFrame, axes: list[str] | None = None) -> pd.DataFrame:
    """
    Add `risk_score` (sum of normalized axes) and `pareto_optimal` (bool).
    Candidates with any NA value among `axes` are excluded from the comparison.

    `axes` defaults to ALL_AXES; pass SEQUENCE_AXES for the pre-structure prune.
    """
    axes = axes or ALL_AXES
    df = df.copy()
    scored = df.dropna(subset=axes)
    if scored.empty:
        df["risk_score"] = pd.NA
        df["pareto_optimal"] = pd.NA
        return df

    # Normalize each axis to [0,1] for the convenience composite.
    norm = scored[axes].copy()
    for x in axes:
        col = pd.to_numeric(norm[x], errors="coerce")
        span = col.max() - col.min()
        norm[x] = 0.0 if span == 0 else (col - col.min()) / span
    df.loc[scored.index, "risk_score"] = norm.sum(axis=1)

    recs = {i: {x: float(scored.loc[i, x]) for x in axes} for i in scored.index}
    pareto = {}
    for i in scored.index:
        pareto[i] = not any(_dominates(recs[j], recs[i], axes) for j in scored.index if j != i)
    df["pareto_optimal"] = df.index.map(lambda i: pareto.get(i, pd.NA))

    return df.sort_values("risk_score", na_position="last").reset_index(drop=True)
