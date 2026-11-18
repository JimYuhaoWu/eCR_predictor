"""Zn ion count in AF3 inputs (F3 — the span-local finger count).

`zinc_finger_count` counts every ZN_FING feature over the whole protein, so on
the fusion path — which folds the extracted span — it over-counts. eCR_mod_lib
writes `dbd_zinc_finger_count`, the fingers inside the span; af3.py must prefer
it and fall back loudly.

No AF3 backend: only the JSON writer runs.
"""
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ecr_predictor.af3 import _write_job_jsons  # noqa: E402


def _zinc_ids(json_path):
    doc = json.loads(Path(json_path).read_text(encoding="utf-8"))
    ligands = [s["ligand"] for s in doc["sequences"] if "ligand" in s]
    if not ligands:
        return []
    assert ligands[0]["ccdCodes"] == ["ZN"]
    return ligands[0]["id"]


def _row(**kw):
    base = {"gene_name": "ZFPM2", "sequence_aa": "ACDEFGHIK" * 5,
            "tf_family": "zf-C2H2"}
    base.update(kw)
    return pd.DataFrame([base])


def test_span_local_count_is_preferred(tmp_path):
    jobs = _write_job_jsons(_row(zinc_finger_count=8, dbd_zinc_finger_count=3),
                            "ACGTACGTACGT", tmp_path, top_n=1)
    assert len(_zinc_ids(jobs[0][1])) == 3


def test_falls_back_to_the_protein_count_with_a_warning(tmp_path, capsys):
    jobs = _write_job_jsons(_row(zinc_finger_count=8), "ACGTACGTACGT",
                            tmp_path, top_n=1)
    assert len(_zinc_ids(jobs[0][1])) == 8
    err = capsys.readouterr().err
    assert "dbd_zinc_finger_count" in err and "over-count" in err


def test_zero_span_count_is_honoured_not_treated_as_missing(tmp_path):
    """A span that kept no keep-type finger gets no ions, even though the
    protein-wide count is non-zero — 0 is an answer, not a null."""
    jobs = _write_job_jsons(_row(zinc_finger_count=8, dbd_zinc_finger_count=0),
                            "ACGTACGTACGT", tmp_path, top_n=1)
    assert _zinc_ids(jobs[0][1]) == []


def test_non_zinc_family_gets_no_ions(tmp_path):
    jobs = _write_job_jsons(_row(tf_family="HMG", zinc_finger_count=8,
                                 dbd_zinc_finger_count=3),
                            "ACGTACGTACGT", tmp_path, top_n=1)
    assert _zinc_ids(jobs[0][1]) == []
