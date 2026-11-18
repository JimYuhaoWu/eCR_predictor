"""Unit tests for DBD loading in fuse.py (F1 — refuse, don't silently fall back).

Two classes of record must be refused by default, each with its own override:
a record with no extracted `dbd_sequence_aa` (which used to be fused on the
whole TF) and a record whose span carries `dbd_review_flag = 1`.

Synthetic TSVs only — no library, no DB, no network.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import fuse  # noqa: E402

COLS = ["gene_name", "sequence_aa", "dbd_sequence_aa", "dbd_review_flag"]


def _tsv(tmp_path, rows, cols=COLS):
    """rows: list of tuples matching `cols`. None becomes an empty TSV field."""
    p = tmp_path / "dbds.tsv"
    pd.DataFrame(rows, columns=cols).to_csv(p, sep="\t", index=False)
    return p


def test_clean_span_is_used(tmp_path):
    path = _tsv(tmp_path, [("SOX2", "FULLSOX2SEQ", "SPAN", 0)])
    dbds, meta = fuse._load_dbds(path, None)
    assert dbds == {"SOX2": "SPAN"}
    # Read as strings (the TSV is loaded with dtype=str), carried unparsed.
    assert str(meta["SOX2"]["dbd_review_flag"]) == "0"


def test_flagged_span_is_skipped_by_default(tmp_path, capsys):
    path = _tsv(tmp_path, [("SOX2", "FULLSOX2SEQ", "SPAN", 0),
                           ("ZNF541", "FULLZNF541SEQ", "FLAGGEDSPAN", 1)])
    dbds, meta = fuse._load_dbds(path, None)
    assert dbds == {"SOX2": "SPAN"}
    assert "ZNF541" not in meta
    err = capsys.readouterr().err
    assert "dbd_review_flag=1" in err and "ZNF541" in err
    assert "--allow-flagged-dbd" in err


def test_flagged_span_is_used_under_the_override(tmp_path, capsys):
    path = _tsv(tmp_path, [("ZNF541", "FULLZNF541SEQ", "FLAGGEDSPAN", 1)])
    dbds, _ = fuse._load_dbds(path, None, allow_flagged=True)
    assert dbds == {"ZNF541": "FLAGGEDSPAN"}
    # Used, but not silently: the override still reports what it let through.
    assert "ZNF541" in capsys.readouterr().err


def test_missing_span_is_skipped_not_fused_full_length(tmp_path, capsys):
    """The hazard F1 exists for: flag 0, no span, previously fused whole."""
    path = _tsv(tmp_path, [("SOX2", "FULLSOX2SEQ", "SPAN", 0),
                           ("ZNF793", "FULLZNF793WITHKRAB", None, 0)])
    dbds, meta = fuse._load_dbds(path, None)
    assert dbds == {"SOX2": "SPAN"}
    assert "ZNF793" not in meta
    err = capsys.readouterr().err
    assert "no extracted span" in err and "ZNF793" in err
    assert "--allow-fulllength-dbd" in err


def test_missing_span_with_null_flag_is_also_skipped(tmp_path):
    """The 430 records whose dbd_review_flag is NULL, not 0."""
    path = _tsv(tmp_path, [("SOX2", "FULLSOX2SEQ", "SPAN", 0),
                           ("GENEX", "FULLGENEX", None, None)])
    dbds, _ = fuse._load_dbds(path, None, allow_flagged=True)
    assert dbds == {"SOX2": "SPAN"}


def test_missing_span_is_used_under_the_override(tmp_path, capsys):
    path = _tsv(tmp_path, [("ZNF793", "FULLZNF793WITHKRAB", None, 0)])
    dbds, _ = fuse._load_dbds(path, None, allow_fulllength=True)
    assert dbds == {"ZNF793": "FULLZNF793WITHKRAB"}
    assert "full-length" in capsys.readouterr().err


def test_no_span_column_at_all_is_refused(tmp_path):
    path = _tsv(tmp_path, [("SOX2", "FULLSOX2SEQ")],
                cols=["gene_name", "sequence_aa"])
    with pytest.raises(SystemExit):
        fuse._load_dbds(path, None)


def test_no_span_column_under_the_override(tmp_path):
    path = _tsv(tmp_path, [("SOX2", "FULLSOX2SEQ")],
                cols=["gene_name", "sequence_aa"])
    dbds, _ = fuse._load_dbds(path, None, allow_fulllength=True)
    assert dbds == {"SOX2": "FULLSOX2SEQ"}


def test_everything_refused_is_an_error_not_an_empty_run(tmp_path):
    path = _tsv(tmp_path, [("ZNF793", "FULL1", None, 0),
                           ("ZNF541", "FULL2", "FLAGGEDSPAN", 1)])
    with pytest.raises(SystemExit):
        fuse._load_dbds(path, None)


def test_subset_still_applies(tmp_path):
    path = _tsv(tmp_path, [("SOX2", "FULL1", "SPAN1", 0),
                           ("OCT4", "FULL2", "SPAN2", 0)])
    dbds, _ = fuse._load_dbds(path, "OCT4")
    assert dbds == {"OCT4": "SPAN2"}


def test_a_refused_gene_does_not_resurrect_from_a_later_duplicate_row(tmp_path):
    """Duplicate gene_name rows exist in the library; the first verdict holds,
    so a refusal cannot be undone by a second row for the same gene."""
    path = _tsv(tmp_path, [("SOX2", "FULLSOX2SEQ", "SPAN", 0),
                           ("ZNF541", "FULL1", "FLAGGEDSPAN", 1),
                           ("ZNF541", "FULL2", "OTHERSPAN", 0)])
    dbds, _ = fuse._load_dbds(path, None)
    assert dbds == {"SOX2": "SPAN"}
