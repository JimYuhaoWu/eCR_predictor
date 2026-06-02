"""
JASPAR PWM scanning via BioPython motifs.

Motif fetch order:
  1. BioPython JASPAR2020 local DB (fast, requires jaspar2020 package)
  2. JASPAR REST API (https://jaspar.elixir.no/api/v1/) — no install needed

Unique JASPAR IDs are fetched in parallel; scoring is sequential.
"""
from __future__ import annotations

import io
import math
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from Bio import motifs
from Bio.Seq import Seq

_JASPAR_API = "https://jaspar.elixir.no/api/v1/matrix/{}/?format=jaspar"
_TIMEOUT = 10  # seconds per request
_MAX_WORKERS = 8
_CACHE_DIR = Path(__file__).parents[1] / "jaspar_cache"


def _fetch_via_cache(jaspar_id: str) -> "motifs.Motif | None":
    """Load a pre-fetched .jaspar file from the local cache directory."""
    path = _CACHE_DIR / f"{jaspar_id}.jaspar"
    if not path.exists():
        return None
    try:
        return motifs.read(path.open(), "jaspar")
    except Exception:
        return None


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
    """Fetch a JASPAR motif: local cache → BioPython DB → REST API."""
    motif = _fetch_via_cache(jaspar_id)
    if motif is not None:
        return motif
    motif = _fetch_via_local_db(jaspar_id)
    if motif is not None:
        return motif
    motif = _fetch_via_rest(jaspar_id)
    if motif is None:
        print(f"WARNING: could not fetch JASPAR motif {jaspar_id}", file=sys.stderr)
    return motif


def _fetch_all_parallel(jaspar_ids: list) -> dict:
    """Fetch all unique JASPAR motifs in parallel. Returns {jaspar_id: motif}."""
    unique = [jid for jid in set(jaspar_ids) if jid and not pd.isna(jid)]
    total = len(unique)
    cache: dict = {}

    print(f"  Fetching {total} unique JASPAR motifs ({_MAX_WORKERS} workers)...", file=sys.stderr)

    with ThreadPoolExecutor(max_workers=_MAX_WORKERS) as executor:
        futures = {executor.submit(_fetch_jaspar_pwm, jid): jid for jid in unique}
        done = 0
        failed = []
        for future in as_completed(futures):
            jid = futures[future]
            cache[jid] = future.result()
            done += 1
            if cache[jid] is None:
                failed.append(jid)
            if done % 100 == 0 or done == total:
                print(f"  {done}/{total} fetched...", file=sys.stderr)

    if failed:
        print(f"  WARNING: {len(failed)} motifs failed to fetch: {', '.join(failed)}", file=sys.stderr)
    return cache


def _log_odds_score(motif: "motifs.Motif", sequence: str) -> float:
    """
    Return the maximum log-odds score of a motif scanned against the sequence,
    normalized to [-1, 1] by the theoretical maximum.

    If the sequence is shorter than the motif, we slide the sequence against
    every sub-window of the PSSM of length len(seq) and take the best score.
    This captures cases where the input overlaps the binding site core.
    """
    seq = Seq(sequence.upper())
    n, m = len(seq), len(motif)

    pwm = motif.counts.normalize(pseudocounts=0.5)
    pssm = pwm.log_odds()

    if n >= m:
        # Standard case: slide motif over sequence
        scores = list(np.atleast_1d(pssm.calculate(seq)))
        if not scores:
            return float("nan")
        raw_max = max(scores)
        max_possible = sum(max(pssm[nt][i] for nt in "ACGT") for i in range(m))
    else:
        # Partial match: slide sequence against every n-length sub-window of the PSSM
        best_score = -math.inf
        best_max_possible = 0.0
        for start in range(m - n + 1):
            # Score seq against PSSM positions [start : start+n]
            sub_score = sum(
                pssm[str(seq[j])][start + j] for j in range(n)
                if str(seq[j]) in pssm
            )
            sub_max = sum(
                max(pssm[nt][start + j] for nt in "ACGT") for j in range(n)
            )
            if sub_score > best_score:
                best_score = sub_score
                best_max_possible = sub_max
        raw_max = best_score
        max_possible = best_max_possible

    if max_possible <= 0 or math.isinf(raw_max):
        return float("nan")

    return raw_max / max_possible


def score_dbds(dbds: pd.DataFrame, sequence: str) -> pd.Series:
    """
    Compute motif_score for each DBD row against `sequence`.
    Returns a Series indexed like `dbds`; pd.NA where no motif is available.
    """
    jaspar_ids = dbds["jaspar_id"].tolist()
    cache = _fetch_all_parallel(jaspar_ids)

    scores: list[float | object] = []
    total = len(dbds)

    for i, (_, row) in enumerate(dbds.iterrows(), 1):
        gene = row.get("gene_symbol") or row.get("name") or "?"
        jid = row.get("jaspar_id")

        if not jid or pd.isna(jid):
            scores.append(pd.NA)
            continue

        motif = cache.get(jid)
        if motif is None:
            scores.append(pd.NA)
        else:
            val = _log_odds_score(motif, sequence)
            score = val if not math.isnan(val) else pd.NA
            print(f"  [{i}/{total}] {gene} ({jid}): score = {score}", file=sys.stderr)
            scores.append(score)

    return pd.Series(scores, index=dbds.index, dtype=object)
