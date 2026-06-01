"""
DBD lookup and species matching against the module library.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from scripts.schema import ModuleLibrary

_DB_PATH = Path(__file__).parents[2] / "eCR_mod_lib" / "library" / "module_library.db"


def load_dbds(db_path: Path = _DB_PATH) -> pd.DataFrame:
    """Return all DBD rows from the library as a DataFrame."""
    with ModuleLibrary(db_path) as lib:
        return lib.to_dataframe("DBD")


def match_species(dbds: pd.DataFrame, species: str) -> pd.DataFrame:
    """
    Filter DBDs by species with fuzzy taxonomic fallback.

    Adds a `query_species_match` column: 'exact' for direct hits,
    'other' for genus-level fallback matches.
    """
    species = species.strip()

    exact = dbds[dbds["organism"].str.lower() == species.lower()].copy()
    if not exact.empty:
        exact["query_species_match"] = "exact"
        return exact

    # Genus fallback: first word of the query species
    genus = species.split()[0].lower()
    fallback = dbds[dbds["organism"].str.lower().str.startswith(genus)].copy()
    fallback["query_species_match"] = "other"
    return fallback
