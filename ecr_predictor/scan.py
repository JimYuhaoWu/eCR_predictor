"""
JASPAR PWM scanning via BioPython motifs.

Motif fetch order:
  1. BioPython JASPAR2020 local DB (fast, requires jaspar2020 package)
  2. JASPAR REST API (https://jaspar.elixir.no/api/v1/) — no install needed
"""
from __future__ import annotations

import io
import math
import sys

import numpy as np
import pandas as pd
import requests
from Bio import motifs
from Bio.Seq import Seq

_JASPAR_API = "https://jaspar.elixir.no/api/v1/matrix/{}/?format=jaspar"
_TIMEOUT = 10  # seconds per request


def _fetch_via_local_db(jaspar_id: str) -> "motifs.Motif | None":
    """Try BioPython's bundled JASPAR2020 local database."""
    try:
        with motifs.jaspar.db.JASPAR2020() as jdb:
            return jdb.fetch_motif_by_id(jaspar_id)
    except Exception:
        return None


def _fetch_via_rest(jaspar_id: str) -> "motifs.Motif | None":
    """Fetch a JASPAR motif via the public REST API."""
    try:
        resp = requests.get(_JASPAR_API.format(jaspar_id), timeout=_TIMEOUT)
        resp.raise_for_status()
        return motifs.read(io.StringIO(resp.text), "jaspar")
    except Exception:
        return None


def _fetch_jaspar_pwm(jaspar_id: str) -> "motifs.Motif | None":
    """Fetch a JASPAR motif, trying local DB then REST API."""
    motif = _fetch_via_local_db(jaspar_id)
    if motif is not None:
        return motif
    motif = _fetch_via_rest(jaspar_id)
    if motif is None:
        print(f"WARNING: could not fetch JASPAR motif {jaspar_id}", file=sys.stderr)
    return motif


def _log_odds_score(motif: "motifs.Motif", sequence: str) -> float:
    """
    Return the maximum log-odds score of a motif scanned against the sequence,
    normalized to [−1, 1] by the theoretical maximum for this motif.
    """
    seq = Seq(sequence.upper())
    pwm = motif.counts.normalize(pseudocounts=0.5)
    pssm = pwm.log_odds()

    scores = list(np.atleast_1d(pssm.calculate(seq)))
    if not scores:
        return float("nan")

    raw_max = max(scores)
    max_possible = sum(max(pssm[nt][i] for nt in "ACGT") for i in range(len(motif)))
    if max_possible <= 0:
        return float("nan")

    return raw_max / max_possible


def score_dbds(dbds: pd.DataFrame, sequence: str) -> pd.Series:
    """
    Compute motif_score for each DBD row against `sequence`.
    Returns a Series indexed like `dbds`; pd.NA where no motif is available.
    """
    scores: list[float | object] = []
    cache: dict[str, "motifs.Motif | None"] = {}

    for _, row in dbds.iterrows():
        jid = row.get("jaspar_id")
        if not jid or pd.isna(jid):
            scores.append(pd.NA)
            continue

        if jid not in cache:
            cache[jid] = _fetch_jaspar_pwm(jid)

        motif = cache[jid]
        if motif is None:
            scores.append(pd.NA)
        else:
            val = _log_odds_score(motif, sequence)
            scores.append(val if not math.isnan(val) else pd.NA)

    return pd.Series(scores, index=dbds.index, dtype=object)
