"""Unit tests for ED loading in fuse.py (F8a — the ED name collision).

Library ED names are not unique. These tests pin the behaviour that every
distinct sequence reaches the design run, that the keys do not depend on row
order, and that --ed still accepts a bare library name.

Synthetic DataFrames only — no library, no DB, no network.
"""
import sys
import types
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import fuse  # noqa: E402


def _df(rows):
    """rows: list of (name, sequence_aa)."""
    return pd.DataFrame(rows, columns=["name", "sequence_aa"])


@pytest.fixture
def patch_library(monkeypatch):
    """Stub ecr_predictor.query, which otherwise imports the sibling
    eCR_mod_lib package to reach the library DB. `_load_eds` imports it at call
    time, so a sys.modules entry is enough and nothing touches disk."""
    def _patch(rows):
        stub = types.ModuleType("ecr_predictor.query")
        stub.load_eds = lambda *a, **k: _df(rows)
        monkeypatch.setitem(sys.modules, "ecr_predictor.query", stub)
    return _patch


def test_unique_names_keep_their_name(patch_library):
    patch_library([("VP64", "AAA"), ("KRAB", "CCC")])
    eds = fuse._load_eds(None)
    assert eds == {"VP64": "AAA", "KRAB": "CCC"}


def test_collision_keeps_every_distinct_sequence(patch_library, capsys):
    patch_library([("DUAL", "AAA"), ("DUAL", "CCC"), ("SOLO", "GGG")])
    eds = fuse._load_eds(None)

    # the colliding name is not a key on its own any more
    assert "DUAL" not in eds
    # both sequences survive; the old setdefault kept only "AAA"
    assert sorted(eds.values()) == ["AAA", "CCC", "GGG"]
    assert sum(1 for k in eds if k.startswith("DUAL~")) == 2
    assert eds["SOLO"] == "GGG"

    err = capsys.readouterr().err
    assert "1 ED name(s) carry >1 distinct sequence" in err
    assert "DUAL" in err


def test_exact_duplicate_rows_collapse(patch_library):
    patch_library([("SAME", "AAA"), ("SAME", "AAA")])
    eds = fuse._load_eds(None)
    assert eds == {"SAME": "AAA"}


def test_keys_do_not_depend_on_row_order(patch_library):
    patch_library([("DUAL", "AAA"), ("DUAL", "CCC")])
    forward = fuse._load_eds(None)
    patch_library([("DUAL", "CCC"), ("DUAL", "AAA")])
    reversed_ = fuse._load_eds(None)
    assert forward == reversed_


def test_no_row_is_lost(patch_library):
    rows = [("A", "AAA"), ("A", "CCC"), ("A", "GGG"),
            ("B", "TTT"), ("B", "TTT"), ("C", "AAA")]
    patch_library(rows)
    eds = fuse._load_eds(None)
    # 5 distinct (name, sequence) pairs; keying on name alone gave 3
    assert len(eds) == len({(n, s) for n, s in rows}) == 5


def test_subset_by_bare_name_selects_all_its_sequences(patch_library):
    patch_library([("DUAL", "AAA"), ("DUAL", "CCC"), ("OTHER", "GGG")])
    eds = fuse._load_eds("DUAL")
    assert sorted(eds.values()) == ["AAA", "CCC"]
    assert all(k.startswith("DUAL~") for k in eds)


def test_subset_by_suffixed_key_selects_one_sequence(patch_library):
    patch_library([("DUAL", "AAA"), ("DUAL", "CCC")])
    all_eds = fuse._load_eds(None)
    one_key = next(k for k, v in all_eds.items() if v == "CCC")

    patch_library([("DUAL", "AAA"), ("DUAL", "CCC")])
    eds = fuse._load_eds(one_key)
    assert eds == {one_key: "CCC"}


def test_subset_of_a_unique_name_is_unchanged(patch_library, capsys):
    patch_library([("VP64", "AAA"), ("KRAB", "CCC")])
    eds = fuse._load_eds("VP64")
    assert eds == {"VP64": "AAA"}
    # no collision among the selected keys -> no warning
    assert "carry >1 distinct sequence" not in capsys.readouterr().err
