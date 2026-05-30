"""
Format and emit the results table.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

_OUTPUT_COLUMNS = [
    "gene_name",
    "species",
    "query_species_match",
    "tf_family",
    "validation_level",
    "motif_score",
    "annotation_confidence",
    "jaspar_id",
]

_COLUMN_ALIASES = {
    "gene_symbol": "gene_name",
    "organism": "species",
    "subtype": "tf_family",
}


def build_result_table(
    dbds: pd.DataFrame,
    motif_scores: pd.Series,
    annotation_confidence: pd.Series,
    include_sequence: bool = False,
) -> pd.DataFrame:
    df = dbds.rename(columns=_COLUMN_ALIASES).copy()
    df["motif_score"] = motif_scores.values
    df["annotation_confidence"] = annotation_confidence.values

    cols = _OUTPUT_COLUMNS.copy()
    if include_sequence:
        cols = cols + ["sequence_aa"]

    # Ensure all output columns exist
    for col in cols:
        if col not in df.columns:
            df[col] = pd.NA

    df = df[cols].copy()

    # Sort: exact matches first, then by motif_score descending (NA last)
    df["_match_rank"] = (df["query_species_match"] != "exact").astype(int)
    df["_score_sort"] = pd.to_numeric(df["motif_score"], errors="coerce")
    df = df.sort_values(["_match_rank", "_score_sort"], ascending=[True, False])
    df = df.drop(columns=["_match_rank", "_score_sort"]).reset_index(drop=True)

    return df


def write_output(df: pd.DataFrame, out_path: str | None) -> None:
    tsv = df.to_csv(sep="\t", index=False, na_rep="NA")
    if out_path:
        Path(out_path).write_text(tsv, encoding="utf-8")
    else:
        sys.stdout.write(tsv)
