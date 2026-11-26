#!/usr/bin/env python
"""
ECR_predictor accessibility-preview stage.

Predicts and visualizes what chromatin accessibility would look like AFTER the
designed eCRs are delivered, as a self-contained interactive HTML page. It is a
*preview/reporting* phase: it consumes the eCR_navigator run bundle (nominated
regions + measured `direction`) and the transition's measured per-region
accessibility (start vs end aTPM), and renders

  measured start  ->  PREDICTED post-eCR (model)  ->  measured end

The predicted column is NOT data: it is a direct-effect interpolation,
`predicted = start + efficacy * (end - start)`, applied ONLY to eCR-targeted
regions (efficacy is a user assumption exposed as a live slider). Off-target and
secondary/network effects are not modeled. Honest by construction: the two
measured endpoints flank the model so the viewer reads it as an interpolation.

Three views: highlighted delivered constructs (per-locus bars), the nominated
significant regions (heatmap), and a genome overview (banded: targeted to-open /
untargeted background / targeted to-close).

Inputs are generic — any transition. The measured-signal columns are auto-labeled
from the aTPM header (e.g. atpm_MEF / atpm_mES -> "MEF" / "mES").

Exercised by tests/test_preview.py on a synthetic bundle + aTPM.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import random
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ecr_predictor.target import read_nominations, read_manifest  # noqa: E402


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_atpm(path):
    """Read a per-region accessibility TSV: chrom,start,end,atpm_<START>,atpm_<END>.
    Returns (values, start_label, end_label) where values[(chrom,start,end)]=(s,e)."""
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        atpm_cols = [c for c in rd.fieldnames if c.lower().startswith("atpm")]
        if len(atpm_cols) < 2:
            raise SystemExit(f"{path}: need two atpm_* columns; found {atpm_cols}")
        sc, ec = atpm_cols[0], atpm_cols[1]
        s_lab = sc.split("_", 1)[1] if "_" in sc else "start"
        e_lab = ec.split("_", 1)[1] if "_" in ec else "end"
        vals = {}
        for r in rd:
            vals[(r["chrom"], int(r["start"]), int(r["end"]))] = (float(r[sc]), float(r[ec]))
    return vals, s_lab, e_lab


def load_bed_keys(path):
    """Region keys (chrom,start,end) from a BED-like file (first 3 cols)."""
    keys = set()
    with open(path) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) >= 3 and p[1].isdigit():
                keys.add((p[0], int(p[1]), int(p[2])))
    return keys


def load_gene_footprints(path, nbhd=0):
    """gene symbol -> {chrom:[min_start,max_end]} from a BED (chrom,start,end,name)."""
    fp = {}
    with open(path) as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 4 or not p[1].isdigit():
                continue
            sym = p[3]
            c, s, e = p[0], int(p[1]) - nbhd, int(p[2]) + nbhd
            d = fp.setdefault(sym, {})
            if c not in d:
                d[c] = [s, e]
            else:
                d[c][0] = min(d[c][0], s); d[c][1] = max(d[c][1], e)
    return fp


def gene_of(fp, chrom, start, end):
    for sym, chroms in fp.items():
        sp = chroms.get(chrom)
        if sp and start < sp[1] and end > sp[0]:
            return sym
    return ""


def _read_map(path):
    """Two-column TSV -> {int(col0): col1}. Tolerates a header."""
    out = {}
    if not path:
        return out
    with open(path) as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0].lstrip("-").isdigit():
                out[int(p[0])] = p[1]
    return out


def _read_ranks(path):
    if not path:
        return None
    txt = Path(path).read_text().replace(",", " ")
    return [int(x) for x in txt.split() if x.lstrip("-").isdigit()]


# ---------------------------------------------------------------------------
# viz-data assembly
# ---------------------------------------------------------------------------

def build_viz_data(bundle, atpm, *, targets=None, genes=None, highlight=None,
                   annot=None, focused_ranks=None, focused_max=30, title=None,
                   ov_open=220, ov_close=150, ov_bg=250, seed=7):
    """Assemble the JSON payload the HTML template renders."""
    rng = random.Random(seed)
    vals, s_lab, e_lab = load_atpm(atpm)
    noms = read_nominations(bundle)
    manifest = {}
    try:
        manifest = read_manifest(bundle)
    except FileNotFoundError:
        pass

    key = lambda n: (n.chrom, n.start, n.end)
    targeted = load_bed_keys(targets) if targets else {key(n) for n in noms}
    fp = load_gene_footprints(genes) if genes else {}
    annot_map = _read_map(annot)
    hl_map = _read_map(highlight)                 # rank -> label (e.g. gene) ; optional
    hl_sub = {}                                   # rank -> sublabel (col 3)
    if highlight:
        with open(highlight) as fh:
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                if len(p) >= 3 and p[0].lstrip("-").isdigit():
                    hl_sub[int(p[0])] = p[2]

    by_rank = {n.rank: n for n in noms}

    # ---- focused rows ----
    if focused_ranks:
        chosen = [by_rank[r] for r in focused_ranks if r in by_rank]
    else:
        pool = [n for n in noms if (not genes or gene_of(fp, *key(n)))]
        chosen = pool[:focused_max] if not genes else pool  # gene-mapped can exceed; capped below
    focused = []
    for n in chosen:
        k = key(n)
        if k not in vals:
            continue
        s, e = vals[k]
        d = n.direction or 0.0
        focused.append({"rank": n.rank, "coord": n.id,
                        "gene": gene_of(fp, *k) if fp else "",
                        "sign": "open" if d > 0 else "close",
                        "start": round(s, 3), "end": round(e, 3),
                        "flag": annot_map.get(n.rank, f"#{n.rank}"),
                        "hl": n.rank in hl_map})
    focused.sort(key=lambda x: (x["sign"] != "open", x["rank"]))
    op = [x for x in focused if x["sign"] == "open"][:focused_max]
    cl = [x for x in focused if x["sign"] == "close"][:focused_max]
    focused = op + cl

    # ---- highlight cards ----
    hl = []
    for rk, lab in hl_map.items():
        n = by_rank.get(rk)
        if not n or key(n) not in vals:
            continue
        s, e = vals[key(n)]
        d = n.direction or 0.0
        hl.append({"rank": rk, "label": lab, "sub": hl_sub.get(rk, ""),
                   "coord": n.id, "sign": "open" if d > 0 else "close",
                   "start": round(s, 3), "end": round(e, 3)})
    hl.sort(key=lambda x: (x["sign"] != "open", x["rank"]))

    # ---- overview bands ----
    rows = [(k in targeted, mv, ev, ev - mv) for k, (mv, ev) in vals.items()]
    open_t = sorted([r for r in rows if r[0] and r[3] > 0], key=lambda r: -r[3])[:ov_open]
    close_t = sorted([r for r in rows if r[0] and r[3] < 0], key=lambda r: r[3])[:ov_close]
    bg_pool = [r for r in rows if not r[0] and abs(r[3]) < 0.15]
    bg = rng.sample(bg_pool, min(ov_bg, len(bg_pool))) if bg_pool else []
    def pack(lst, band):
        return [{"start": round(m, 3), "end": round(e, 3), "targeted": t, "band": band}
                for (t, m, e, _d) in lst]
    overview = (pack(open_t, "to-open (targeted)") + pack(bg, "background (not targeted)")
                + pack(close_t, "to-close (targeted)"))

    tr = manifest.get("transition", {}) if isinstance(manifest, dict) else {}
    if not title:
        frm = tr.get("state_a") or s_lab
        to = tr.get("state_b") or e_lab
        title = f"Predicted chromatin accessibility after eCR delivery — {frm}→{to}"

    return {
        "title": title,
        "subtitle": ("Measured {s} and {e} accessibility flank a modeled post-delivery state. "
                     "Move the efficacy control to see how far the engineered chromatin regulators "
                     "would push their target regions from start toward target."
                     ).format(s=s_lab, e=e_lab),
        "start_label": s_lab, "end_label": e_lab,
        "focused": focused, "highlight": hl, "overview": overview,
        "n_union": len(vals), "n_targets": len(targeted),
        "provenance": [
            f"Measured accessibility: per-region aTPM ({Path(atpm).name}), {len(vals):,} regions "
            f"({s_lab} vs {e_lab}).",
            f"Targets: {len(targeted):,} eCR-delivered regions from the navigator run bundle "
            f"({Path(bundle).name}).",
            "Predicted state: direct-effect interpolation on targeted regions only; efficacy is a "
            "user assumption, not a fitted parameter.",
        ],
        "note": "Predicted = start + efficacy*(end-start) for TARGETED regions only (direct-effect). "
                "Endpoints are measured; predicted is a model, not data.",
    }


# ---------------------------------------------------------------------------
# render
# ---------------------------------------------------------------------------

def render_html(data):
    return _TEMPLATE.replace("__DATA__", json.dumps(data))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", required=True, help="eCR_navigator run bundle dir (nominations + manifest).")
    ap.add_argument("--atpm",
                    help="Per-region accessibility TSV: chrom,start,end,atpm_<START>,atpm_<END>. "
                         "Defaults to <bundle>/accessibility.tsv if present (the arg overrides the bundle).")
    ap.add_argument("--output", default="accessibility_preview.html", help="Output HTML path.")
    ap.add_argument("--targets", help="BED of eCR-delivered regions (default: all nominations).")
    ap.add_argument("--genes", help="Gene BED (chrom,start,end,name) for heatmap labels.")
    ap.add_argument("--highlight", help="TSV rank<TAB>label[<TAB>sublabel] -> per-locus cards + marks.")
    ap.add_argument("--annot", help="TSV rank<TAB>flag -> heatmap flag column (e.g. '3/3 routes').")
    ap.add_argument("--focused-ranks", help="File of ranks to show in the heatmap (default: top by rank).")
    ap.add_argument("--focused-max", type=int, default=30, help="Max regions per direction in the heatmap.")
    ap.add_argument("--title", help="Page title (default: from bundle manifest transition).")
    ap.add_argument("--data-out", help="Also write the viz JSON here.")
    args = ap.parse_args(argv)

    # accessibility source: --atpm overrides; else the bundle's accessibility.tsv; else error.
    atpm = args.atpm
    if not atpm:
        cand = Path(args.bundle) / "accessibility.tsv"
        if not cand.exists():
            ap.error(f"no --atpm given and no accessibility.tsv in the bundle ({cand}). "
                     "Pass --atpm <chrom,start,end,atpm_A,atpm_B> (e.g. the navigator's atpm_union.tsv).")
        atpm = str(cand)

    data = build_viz_data(
        args.bundle, atpm, targets=args.targets, genes=args.genes,
        highlight=args.highlight, annot=args.annot,
        focused_ranks=_read_ranks(args.focused_ranks),
        focused_max=args.focused_max, title=args.title)

    # encoding="utf-8" is required, not tidiness: write_text() otherwise uses
    # the platform default, which is GBK on a Chinese-locale Windows box, and
    # render_html() emits U+2212 (the minus sign in the "Predicted = start +
    # efficacy x (target - start)" caption). That raised UnicodeEncodeError and
    # wrote no file at all. The HTML declares utf-8 in its own <meta>, so the
    # bytes have to match.
    if args.data_out:
        Path(args.data_out).write_text(json.dumps(data), encoding="utf-8")
    Path(args.output).write_text(render_html(data), encoding="utf-8")
    print(f"[preview] {len(data['focused'])} focused rows, {len(data['overview'])} overview, "
          f"{len(data['highlight'])} highlighted -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------
# HTML template (self-contained; labels + data injected as __DATA__)
# ---------------------------------------------------------------------------
_TEMPLATE = r"""<title>eCR accessibility preview</title>
<style>
:root{
  --bg:#eef2f3; --surface:#ffffff; --surface-2:#f6f9fa; --text:#13202b; --muted:#586b78;
  --border:#d6dfe4; --accent:#0e8c8c; --accent-ink:#0a6a6a;
  --open:#2f9e6b; --close:#d1556c; --model:#8a6d1f;
  --shadow:0 1px 2px rgba(16,40,54,.06),0 8px 24px rgba(16,40,54,.06);
}
@media (prefers-color-scheme:dark){:root{
  --bg:#0b131b; --surface:#111d27; --surface-2:#0e1922; --text:#e6eef2; --muted:#8ba0ac;
  --border:#22333f; --accent:#33c4c4; --accent-ink:#7fe0e0; --open:#42c088; --close:#e77890; --model:#e2c060;
  --shadow:0 1px 2px rgba(0,0,0,.3),0 10px 30px rgba(0,0,0,.35);
}}
:root[data-theme="light"]{
  --bg:#eef2f3; --surface:#ffffff; --surface-2:#f6f9fa; --text:#13202b; --muted:#586b78;
  --border:#d6dfe4; --accent:#0e8c8c; --accent-ink:#0a6a6a; --open:#2f9e6b; --close:#d1556c; --model:#8a6d1f;
  --shadow:0 1px 2px rgba(16,40,54,.06),0 8px 24px rgba(16,40,54,.06);
}
:root[data-theme="dark"]{
  --bg:#0b131b; --surface:#111d27; --surface-2:#0e1922; --text:#e6eef2; --muted:#8ba0ac;
  --border:#22333f; --accent:#33c4c4; --accent-ink:#7fe0e0; --open:#42c088; --close:#e77890; --model:#e2c060;
  --shadow:0 1px 2px rgba(0,0,0,.3),0 10px 30px rgba(0,0,0,.35);
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--text);
  font-family:ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  line-height:1.5;-webkit-font-smoothing:antialiased;}
.mono{font-family:ui-monospace,"SF Mono","Cascadia Code",Menlo,Consolas,monospace}
.wrap{max-width:1060px;margin:0 auto;padding:32px 24px 64px}
.eyebrow{text-transform:uppercase;letter-spacing:.14em;font-size:12px;font-weight:600;color:var(--accent-ink)}
h1{font-size:30px;line-height:1.15;margin:.25em 0 .1em;font-weight:720;text-wrap:balance;letter-spacing:-.01em}
.sub{color:var(--muted);font-size:16px;max-width:64ch;margin:0}
h2{font-size:19px;margin:0 0 2px;font-weight:680;letter-spacing:-.01em}
.sec-note{color:var(--muted);font-size:13.5px;margin:0 0 16px;max-width:70ch}
section{margin-top:34px}
.card{background:var(--surface);border:1px solid var(--border);border-radius:14px;box-shadow:var(--shadow)}
.pad{padding:20px 22px}
.banner{display:flex;gap:12px;align-items:flex-start;margin-top:20px;padding:14px 16px;border-radius:12px;
  background:var(--surface-2);border:1px solid var(--border);border-left:4px solid var(--model)}
.banner .tag{flex:none;font-size:11px;font-weight:700;letter-spacing:.08em;color:var(--model);
  border:1px solid var(--model);border-radius:6px;padding:3px 7px;margin-top:1px}
.banner p{margin:0;font-size:13.5px;color:var(--muted)}
.banner b{color:var(--text)}
.controls{position:sticky;top:0;z-index:5;margin-top:22px;background:color-mix(in srgb,var(--bg) 86%,transparent);
  backdrop-filter:blur(8px);border:1px solid var(--border);border-radius:14px;padding:16px 20px;
  display:flex;flex-wrap:wrap;gap:20px 28px;align-items:center;box-shadow:var(--shadow)}
.slider-block{display:flex;flex-direction:column;gap:6px;min-width:300px;flex:1}
.slider-row{display:flex;align-items:center;gap:14px}
.slider-row label{font-size:12px;font-weight:600;text-transform:uppercase;letter-spacing:.09em;color:var(--muted);white-space:nowrap}
input[type=range]{-webkit-appearance:none;appearance:none;height:6px;border-radius:99px;flex:1;
  background:linear-gradient(90deg,var(--accent) var(--pct,70%),var(--border) var(--pct,70%));outline:none;cursor:pointer}
input[type=range]::-webkit-slider-thumb{-webkit-appearance:none;width:18px;height:18px;border-radius:50%;
  background:var(--surface);border:3px solid var(--accent);cursor:pointer;box-shadow:0 1px 4px rgba(0,0,0,.25)}
input[type=range]::-moz-range-thumb{width:16px;height:16px;border-radius:50%;background:var(--surface);border:3px solid var(--accent);cursor:pointer}
input[type=range]:focus-visible{box-shadow:0 0 0 3px color-mix(in srgb,var(--accent) 40%,transparent)}
.eff-val{font-size:26px;font-weight:720;color:var(--accent-ink);min-width:74px;text-align:right;font-variant-numeric:tabular-nums}
.presets{display:flex;gap:6px;flex-wrap:wrap}
.presets button{font:inherit;font-size:12px;font-weight:600;padding:5px 10px;border-radius:8px;cursor:pointer;
  background:var(--surface);color:var(--muted);border:1px solid var(--border)}
.presets button:hover{color:var(--text);border-color:var(--accent)}
.presets button.on{background:var(--accent);color:#fff;border-color:var(--accent)}
.legend{display:flex;flex-direction:column;gap:6px}
.legend .gl{height:12px;width:180px;border-radius:4px;border:1px solid var(--border);
  background:linear-gradient(90deg,#0d1b2a,#1b5e7e,#4fa3a1,#e8b23a,#f6e27a)}
.legend .lab{display:flex;justify-content:space-between;font-size:11px;color:var(--muted)}
.grid4{display:grid;grid-template-columns:repeat(4,1fr);gap:14px}
@media(max-width:820px){.grid4{grid-template-columns:repeat(2,1fr)}}
@media(max-width:460px){.grid4{grid-template-columns:1fr}}
.locus{padding:15px 16px}
.locus .gene{font-size:17px;font-weight:700;letter-spacing:-.01em}
.locus .chip{font-size:10.5px;font-weight:700;letter-spacing:.05em;padding:2px 7px;border-radius:999px;text-transform:uppercase}
.chip.open{color:var(--open);background:color-mix(in srgb,var(--open) 14%,transparent)}
.chip.close{color:var(--close);background:color-mix(in srgb,var(--close) 14%,transparent)}
.locus .ed{font-size:12px;color:var(--muted);margin:2px 0 12px;min-height:1em}
.bar-row{display:grid;grid-template-columns:64px 1fr 52px;align-items:center;gap:8px;margin:7px 0}
.bar-row .k{font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:.05em}
.bar-track{height:14px;border-radius:5px;background:var(--surface-2);border:1px solid var(--border);overflow:hidden;position:relative}
.bar-fill{height:100%;transition:width .25s ease}
.bar-row.pred .bar-track{border-style:dashed}
.bar-row .v{font-size:12px;text-align:right;font-variant-numeric:tabular-nums;color:var(--text)}
.coordline{font-size:11px;color:var(--muted);margin-top:10px;border-top:1px dashed var(--border);padding-top:8px}
.hm{display:grid;grid-template-columns:88px 46px repeat(3,1fr);gap:3px;align-items:center}
.hm .gene{font-size:12.5px;font-weight:600;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.hm .gene.r1::before{content:"";display:inline-block;width:7px;height:7px;border-radius:2px;background:var(--accent);margin-right:6px;vertical-align:middle}
.hm .flag{font-size:10px;color:var(--muted);text-align:center;font-variant-numeric:tabular-nums}
.cell{height:26px;border-radius:5px;position:relative;cursor:default}
.cell.model{outline:1.5px dashed color-mix(in srgb,var(--model) 70%,transparent);outline-offset:-2px}
.grouphdr{grid-column:1/-1;display:flex;align-items:center;gap:8px;margin:14px 0 3px;font-size:12px;font-weight:700;letter-spacing:.04em}
.grouphdr .dot{width:9px;height:9px;border-radius:2px}
.grouphdr .n{color:var(--muted);font-weight:500}
.colcap{display:grid;grid-template-columns:88px 46px repeat(3,1fr);gap:3px;margin-bottom:2px}
.colcap .cap{text-align:center;font-size:11px;font-weight:600;color:var(--muted)}
.colcap .cap .mm{display:block;font-size:9.5px;font-weight:600;letter-spacing:.04em;text-transform:uppercase}
.colcap .cap.meas .mm{color:var(--accent-ink)} .colcap .cap.model .mm{color:var(--model)}
.ov-cols{display:grid;grid-template-columns:150px repeat(3,1fr);gap:8px;margin-bottom:6px}
.ov-cols .cap{text-align:center;font-size:11px;font-weight:600;color:var(--muted)}
.band{display:grid;grid-template-columns:150px 1fr;gap:12px;align-items:center;margin:10px 0}
.band .blab .t{font-weight:700;font-size:12px} .band .blab .d{color:var(--muted);font-size:11px}
.band canvas{width:100%;height:auto;border-radius:6px;border:1px solid var(--border);display:block}
.tip{position:fixed;pointer-events:none;z-index:20;background:var(--text);color:var(--bg);font-size:11.5px;
  padding:7px 9px;border-radius:7px;opacity:0;transition:opacity .12s;max-width:240px;line-height:1.4;box-shadow:0 6px 18px rgba(0,0,0,.3)}
footer{margin-top:40px;padding-top:18px;border-top:1px solid var(--border);color:var(--muted);font-size:12.5px}
footer b{color:var(--text)} footer ul{margin:8px 0 0;padding-left:18px} footer li{margin:3px 0}
</style>

<div class="wrap">
  <div class="eyebrow">eCR platform · accessibility preview</div>
  <h1 id="pageTitle"></h1>
  <p class="sub" id="pageSub"></p>

  <div class="banner">
    <span class="tag">MODEL</span>
    <p><b>The "predicted" state is a hypothesis, not a measurement.</b> The
    <span class="lbl-start"></span> and <span class="lbl-end"></span> columns are measured accessibility.
    Predicted = start + efficacy × (target − start), applied <b>only to eCR-targeted regions</b>
    (direct-effect model). Off-target and secondary effects are not modeled — the wet-lab round is what
    measures the real effect.</p>
  </div>

  <div class="controls">
    <div class="slider-block">
      <div class="slider-row">
        <label for="eff">eCR efficacy</label>
        <input type="range" id="eff" min="0" max="100" value="70" step="5">
        <span class="eff-val" id="effv">70%</span>
      </div>
      <div class="presets" id="presets">
        <button data-v="0">0% (start)</button>
        <button data-v="50">50%</button>
        <button data-v="70" class="on">70%</button>
        <button data-v="100">100% (target)</button>
      </div>
    </div>
    <div class="legend">
      <div class="gl"></div>
      <div class="lab"><span>low accessibility</span><span>high</span></div>
    </div>
  </div>

  <section id="hlSection">
    <h2>Delivered constructs</h2>
    <p class="sec-note">Each highlighted eCR and its target locus. Bars are accessibility: measured start,
    the modeled post-eCR level at the current efficacy, and measured target.</p>
    <div class="grid4" id="highlight"></div>
  </section>

  <section>
    <h2>Nominated significant regions</h2>
    <p class="sec-note">Openers then closers. The dashed middle column is the modeled state; endpoints are measured.</p>
    <div class="card pad">
      <div class="colcap">
        <div></div><div></div>
        <div class="cap meas"><span class="mm">measured</span><span class="lbl-start"></span></div>
        <div class="cap model"><span class="mm">model</span>post-eCR</div>
        <div class="cap meas"><span class="mm">measured</span><span class="lbl-end"></span></div>
      </div>
      <div class="hm" id="focused"></div>
    </div>
  </section>

  <section>
    <h2>Genome overview</h2>
    <p class="sec-note">Targeted to-open and to-close regions, plus an untargeted background that stays at
    the start state — the direct-effect model acts only where an eCR is delivered.</p>
    <div class="card pad">
      <div class="ov-cols"><div></div>
        <div class="cap"><span class="lbl-start"></span> <span style="font-weight:500;color:var(--accent-ink)">measured</span></div>
        <div class="cap">post-eCR <span style="font-weight:500;color:var(--model)">model</span></div>
        <div class="cap"><span class="lbl-end"></span> <span style="font-weight:500;color:var(--accent-ink)">measured</span></div>
      </div>
      <div id="overview"></div>
    </div>
  </section>

  <footer><b>Data & method.</b><ul id="foot"></ul></footer>
</div>
<div class="tip" id="tip"></div>

<script>
const DATA = __DATA__;
const $=s=>document.querySelector(s), tip=$("#tip");
const STOPS=[[0,[13,27,42]],[0.35,[27,94,126]],[0.62,[79,163,161]],[0.83,[232,178,58]],[1,[246,226,122]]];
let CAP=1;
(function(){const v=[];DATA.overview.forEach(d=>{v.push(d.start,d.end)});DATA.focused.forEach(d=>{v.push(d.start,d.end)});
  v.sort((a,b)=>a-b);CAP=v[Math.floor(v.length*0.98)]||1;if(CAP<=0)CAP=1;})();
function col(x){let t=Math.max(0,Math.min(1,x/CAP));
  for(let i=1;i<STOPS.length;i++){if(t<=STOPS[i][0]){const[a,ca]=STOPS[i-1],[b,cb]=STOPS[i];
    const f=(t-a)/(b-a);return`rgb(${ca.map((c,j)=>Math.round(c+(cb[j]-c)*f)).join(",")})`;}}
  return`rgb(${STOPS[STOPS.length-1][1].join(",")})`;}
const pred=(s,e,eff,tg)=>tg?s+eff*(e-s):s;
const shortCoord=c=>c.replace(/(-\d+)$/,"").replace(/(\d{3})\d{3,}$/,"$1k");
let EFF=0.70;
function showTip(e,html){tip.innerHTML=html;tip.style.opacity=1;
  tip.style.left=Math.min(e.clientX+14,innerWidth-250)+"px";tip.style.top=(e.clientY+14)+"px";}
function hideTip(){tip.style.opacity=0;}

function buildHighlight(){const host=$("#highlight");if(!DATA.highlight.length){$("#hlSection").style.display="none";return;}
  host.innerHTML="";
  DATA.highlight.forEach(d=>{const el=document.createElement("div");el.className="card locus";
    el.innerHTML=`<div style="display:flex;justify-content:space-between;align-items:center;gap:8px">
      <span class="gene">${d.label}</span><span class="chip ${d.sign}">${d.sign}</span></div>
      <div class="ed">${d.sub||""}</div>
      <div class="bar-row"><span class="k lbl-start"></span><div class="bar-track"><div class="bar-fill" data-b="start"></div></div><span class="v" data-v="start"></span></div>
      <div class="bar-row pred"><span class="k">post-eCR</span><div class="bar-track"><div class="bar-fill" data-b="pred"></div></div><span class="v" data-v="pred"></span></div>
      <div class="bar-row"><span class="k lbl-end"></span><div class="bar-track"><div class="bar-fill" data-b="end"></div></div><span class="v" data-v="end"></span></div>
      <div class="coordline mono">${d.coord}</div>`;
    host.appendChild(el);d._el=el;});
  updHighlight();}
function updHighlight(){DATA.highlight.forEach(d=>{if(!d._el)return;const p=pred(d.start,d.end,EFF,true);
  const set=(b,val)=>{const f=d._el.querySelector(`[data-b=${b}]`);f.style.width=Math.min(100,val/CAP*100)+"%";f.style.background=col(val);
    d._el.querySelector(`[data-v=${b}]`).textContent=val.toFixed(2);};
  set("start",d.start);set("pred",p);set("end",d.end);});}

function buildFocused(){const host=$("#focused");host.innerHTML="";
  [["open","to-open (openers)","var(--open)"],["close","to-close (closers)","var(--close)"]].forEach(([sign,label,c])=>{
    const rows=DATA.focused.filter(d=>d.sign===sign);if(!rows.length)return;
    const gh=document.createElement("div");gh.className="grouphdr";
    gh.innerHTML=`<span class="dot" style="background:${c}"></span>${label} <span class="n">· ${rows.length} regions</span>`;
    host.appendChild(gh);
    rows.forEach(d=>{const g=document.createElement("div");g.className="gene"+(d.hl?" r1":"");
      g.textContent=d.gene||shortCoord(d.coord);g.title=d.coord;
      const fl=document.createElement("div");fl.className="flag";fl.textContent=d.flag||"";
      host.appendChild(g);host.appendChild(fl);
      const mk=(key,model)=>{const cel=document.createElement("div");cel.className="cell"+(model?" model":"");
        cel.dataset.k=key;cel._d=d;host.appendChild(cel);return cel;};
      d._cells=[mk("start",false),mk("pred",true),mk("end",false)];});
  });
  host.addEventListener("mousemove",e=>{const cel=e.target.closest(".cell");if(!cel){hideTip();return;}
    const d=cel._d,k=cel.dataset.k,val=k==="start"?d.start:k==="end"?d.end:pred(d.start,d.end,EFF,true);
    const lab=k==="pred"?`post-eCR (model, ${Math.round(EFF*100)}%)`:k==="start"?DATA.start_label+" (measured)":DATA.end_label+" (measured)";
    showTip(e,`<b>${d.gene||shortCoord(d.coord)}</b> · ${d.flag}<br><span class="mono">${d.coord}</span><br>${lab}: <b>${val.toFixed(3)}</b>`);});
  host.addEventListener("mouseleave",hideTip);updFocused();}
function updFocused(){DATA.focused.forEach(d=>{if(!d._cells)return;
  d._cells[0].style.background=col(d.start);
  d._cells[1].style.background=col(pred(d.start,d.end,EFF,true));
  d._cells[2].style.background=col(d.end);});}

const BDESC={"to-open (targeted)":"eCR opens these","background (not targeted)":"no eCR — stays start",
  "to-close (targeted)":"eCR closes these"};
let ovCanvases=[];
function buildOverview(){const host=$("#overview");host.innerHTML="";ovCanvases=[];
  ["to-open (targeted)","background (not targeted)","to-close (targeted)"].forEach(band=>{
    const rows=DATA.overview.filter(d=>d.band===band);if(!rows.length)return;
    const wrap=document.createElement("div");wrap.className="band";
    const lab=document.createElement("div");lab.className="blab";
    lab.innerHTML=`<div class="t">${band.replace(" (targeted)","").replace(" (not targeted)","")}</div><div class="d">${rows.length} regions · ${BDESC[band]}</div>`;
    const cv=document.createElement("canvas");wrap.appendChild(lab);wrap.appendChild(cv);host.appendChild(wrap);
    ovCanvases.push({cv,rows});});
  window.addEventListener("resize",drawOverview);drawOverview();}
function drawOverview(){const dpr=Math.min(devicePixelRatio||1,2);
  ovCanvases.forEach(({cv,rows})=>{const cssW=cv.clientWidth||600,rowH=Math.max(2,Math.min(6,Math.floor(360/rows.length)));
    const cssH=rows.length*rowH;cv.style.height=cssH+"px";cv.width=cssW*dpr;cv.height=cssH*dpr;
    const ctx=cv.getContext("2d");ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,cssW,cssH);
    const gap=Math.round(cssW*0.02),colW=(cssW-2*gap)/3;
    rows.forEach((d,i)=>{const y=i*rowH,p=pred(d.start,d.end,EFF,d.targeted);
      [d.start,p,d.end].forEach((v,c)=>{ctx.fillStyle=col(v);ctx.fillRect(c*(colW+gap),y,colW,Math.max(1,rowH-0.4));});});
  });}

function setEff(v){EFF=v/100;$("#eff").value=v;$("#eff").style.setProperty("--pct",v+"%");$("#effv").textContent=v+"%";
  document.querySelectorAll("#presets button").forEach(b=>b.classList.toggle("on",+b.dataset.v===v));
  updHighlight();updFocused();drawOverview();}

document.title=DATA.title;
$("#pageTitle").textContent=DATA.title;
$("#pageSub").textContent=DATA.subtitle;
document.querySelectorAll(".lbl-start").forEach(e=>e.textContent=DATA.start_label);
document.querySelectorAll(".lbl-end").forEach(e=>e.textContent=DATA.end_label);
$("#foot").innerHTML=DATA.provenance.map(x=>`<li>${x}</li>`).join("");
$("#eff").addEventListener("input",e=>setEff(+e.target.value));
$("#presets").addEventListener("click",e=>{const b=e.target.closest("button");if(b)setEff(+b.dataset.v);});
buildHighlight();buildFocused();buildOverview();setEff(70);
</script>
"""
