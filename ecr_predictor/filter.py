"""
Filter the predictor output table before downstream validation.

Drops rows that are BOTH low annotation confidence AND below a motif_score
threshold. Rows with motif_score == NA are kept if annotation_confidence
is not 'low' (they may still be structurally interesting).
"""
from __future__ import annotations

import pandas as pd


def filter_hits(
    df: pd.DataFrame,
    min_motif_score: float = 0.0,
    drop_low_confidence: bool = True,
) -> pd.DataFrame:
    """
    Return a filtered copy of the predictor output DataFrame.

    Rows are dropped when BOTH conditions are true:
      - annotation_confidence == 'low'
      - motif_score < min_motif_score  (or motif_score is NA)

    Parameters
    ----------
    df : predictor output, must have 'motif_score' and 'annotation_confidence' cols
    min_motif_score : log-odds cutoff; default 0.0 (log-odds boundary)
    drop_low_confidence : set False to skip confidence filtering entirely
    """
    if not drop_low_confidence:
        return df.copy()

    score = pd.to_numeric(df["motif_score"], errors="coerce")
    is_low_conf = df["annotation_confidence"] == "low"
    below_threshold = score < min_motif_score  # NA → True after coerce
    drop_mask = is_low_conf & below_threshold

    filtered = df[~drop_mask].reset_index(drop=True)
    return filtered
