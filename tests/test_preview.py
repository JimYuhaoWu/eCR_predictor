"""Unit tests for the accessibility-preview stage (ecr_predictor/preview.py).

Runs on a synthetic bundle + aTPM — no navigator run, no network, no external tools.
"""
import json
import tempfile
from pathlib import Path

import pytest

from ecr_predictor.preview import build_viz_data, render_html, load_atpm, main


def _bundle(tmp):
    b = Path(tmp)
    (b / "manifest.json").write_text(json.dumps(
        {"transition": {"state_a": "MEF", "state_b": "mES", "label": "MEF -> mES"}}))
    (b / "nominations.tsv").write_text(
        "chrom\tstart\tend\trank\tnomination_score\tdirection\n"
        "chr1\t100\t200\t1\t1.0\t0.8\n"       # opener
        "chr2\t300\t400\t2\t0.99\t-0.7\n")    # closer
    return str(b)


def _atpm(tmp, name="atpm.tsv"):
    p = Path(tmp) / name
    p.write_text(
        "chrom\tstart\tend\tatpm_MEF\tatpm_mES\n"
        "chr1\t100\t200\t0.10\t0.90\n"        # opener: start<end
        "chr2\t300\t400\t0.80\t0.20\n"        # closer: start>end
        "chr3\t500\t600\t0.50\t0.52\n")       # untargeted background
    return str(p)


def test_labels_from_header():
    with tempfile.TemporaryDirectory() as tmp:
        vals, s, e = load_atpm(_atpm(tmp))
        assert (s, e) == ("MEF", "mES")
        assert len(vals) == 3


def test_build_viz_data():
    with tempfile.TemporaryDirectory() as tmp:
        d = build_viz_data(_bundle(tmp), _atpm(tmp))
        assert d["start_label"] == "MEF" and d["end_label"] == "mES"
        assert "MEF→mES" in d["title"]                 # manifest state_a/state_b drive the title
        assert d["n_union"] == 3 and d["n_targets"] == 2
        assert {r["sign"] for r in d["focused"]} == {"open", "close"}
        opener = next(r for r in d["focused"] if r["sign"] == "open")
        assert opener["start"] == 0.1 and opener["end"] == 0.9   # measured endpoints preserved
        # overview: the untargeted chr3 is a background band, targeted regions are banded by sign
        bands = {r["band"] for r in d["overview"]}
        assert any("background" in b for b in bands)
        assert all(not r["targeted"] for r in d["overview"] if "background" in r["band"])


def test_render_inlines_data():
    with tempfile.TemporaryDirectory() as tmp:
        html = render_html(build_viz_data(_bundle(tmp), _atpm(tmp)))
        assert "__DATA__" not in html          # placeholder was substituted
        assert '"focused"' in html and "MEF" in html


def test_atpm_resolution_arg_overrides_bundle(tmp_path, capsys):
    # arg wins even when the bundle has accessibility.tsv
    b = Path(_bundle(tmp_path))
    (b / "accessibility.tsv").write_text(
        "chrom\tstart\tend\tatpm_A\tatpm_B\nchr1\t100\t200\t0\t1\n")
    override = _atpm(tmp_path, "override.tsv")
    out = tmp_path / "o.html"
    main(["--bundle", str(b), "--atpm", override, "--output", str(out)])
    assert "MEF" in out.read_text(encoding="utf-8")   # used override (MEF/mES), not the bundle's A/B


def test_atpm_resolution_falls_back_to_bundle(tmp_path):
    b = Path(_bundle(tmp_path))
    # no --atpm, no accessibility.tsv -> error
    with pytest.raises(SystemExit):
        main(["--bundle", str(b), "--output", str(tmp_path / "x.html")])
    # add accessibility.tsv to the bundle -> resolves without --atpm
    Path(_atpm(tmp_path)).replace(b / "accessibility.tsv")
    out = tmp_path / "o2.html"
    main(["--bundle", str(b), "--output", str(out)])
    assert out.exists() and "__DATA__" not in out.read_text(encoding="utf-8")
