"""
Map validation_level → annotation_confidence label.
"""
from __future__ import annotations

import pandas as pd

_LEVEL_MAP = {
    "screen-validated": "high",
    "ChIP-validated": "high",
    "structurally-resolved": "high",
    "motif-only": "medium",
    "predicted": "low",
}


def assign_annotation_confidence(dbds: pd.DataFrame) -> pd.Series:
    """Return a Series of annotation_confidence labels for each DBD row."""
    return dbds["validation_level"].map(_LEVEL_MAP).fillna("low")
