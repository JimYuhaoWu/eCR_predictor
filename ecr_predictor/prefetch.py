"""
Pre-fetch all JASPAR motifs referenced in the module library and store them
locally as .jaspar files under jaspar_cache/.

Run once on the server after setup:
    python -m ecr_predictor.prefetch --db path/to/module_library.db
"""
from __future__ import annotations

import argparse
import io
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests
from Bio import motifs

from scripts.schema import ModuleLibrary

_JASPAR_API = "https://jaspar.elixir.no/api/v1/matrix/{}/?format=jaspar"
_TIMEOUT = 10
_MAX_WORKERS = 8

DEFAULT_CACHE = Path(__file__).parents[1] / "jaspar_cache"
DEFAULT_DB = Path(__file__).parents[2] / "eCR_mod_lib" / "library" / "module_library.db"


def _fetch_one(jaspar_id: str) -> "tuple[str, motifs.Motif | None]":
    try:
        resp = requests.get(_JASPAR_API.format(jaspar_id), timeout=_TIMEOUT)
        resp.raise_for_status()
        return jaspar_id, motifs.read(io.StringIO(resp.text), "jaspar")
    except Exception as e:
        print(f"  WARNING: failed to fetch {jaspar_id}: {e}", file=sys.stderr)
        return jaspar_id, None


def prefetch(db_path: Path, cache_dir: Path) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)

    with ModuleLibrary(db_path) as lib:
        df = lib.to_dataframe("DBD")

    jaspar_ids = [jid for jid in df["jaspar_id"].dropna().unique() if jid]
    total = len(jaspar_ids)

    # Skip IDs already cached
    pending = [jid for jid in jaspar_ids if not (cache_dir / f"{jid}.jaspar").exists()]
    already = total - len(pending)
    if already:
        print(f"{already}/{total} motifs already cached, fetching {len(pending)} new ones.", file=sys.stderr)

    if not pending:
        print("All motifs already cached.", file=sys.stderr)
        return

    saved = 0
    with ThreadPoolExecutor(max_workers=_MAX_WORKERS) as executor:
        futures = {executor.submit(_fetch_one, jid): jid for jid in pending}
        done = 0
        for future in as_completed(futures):
            jid, motif = future.result()
            done += 1
            if motif is not None:
                out = cache_dir / f"{jid}.jaspar"
                out.write_text(format(motif, "jaspar"), encoding="utf-8")
                saved += 1
                print(f"  [{done}/{len(pending)}] {jid}: saved.", file=sys.stderr)
            else:
                print(f"  [{done}/{len(pending)}] {jid}: FAILED.", file=sys.stderr)

    print(f"\nDone. {saved}/{len(pending)} motifs saved to {cache_dir}/", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(description="Pre-fetch JASPAR motifs for the module library.")
    parser.add_argument("--db", default=str(DEFAULT_DB), help="Path to module_library.db")
    parser.add_argument("--cache-dir", default=str(DEFAULT_CACHE), help="Directory to store .jaspar files")
    args = parser.parse_args()
    prefetch(Path(args.db), Path(args.cache_dir))


if __name__ == "__main__":
    main()
