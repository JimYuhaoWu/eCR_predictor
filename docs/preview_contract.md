# Preview stage — contract (`preview.py`)

The **accessibility-preview** is a read-only reporting stage of `eCR_predictor`. After
the design chain (`target.py → cli.py → refine.py → fuse.py`) it renders a
**self-contained, interactive HTML** showing what chromatin accessibility would look
like *after* the designed eCRs are delivered:

```
measured start  →  PREDICTED post-eCR (model)  →  measured end
```

It runs on CPU, needs no external tools/HPCC, and consumes existing artifacts.

## The model (honest by construction)

The predicted column is **not data**. It is a direct-effect interpolation

```
predicted = start + efficacy · (end − start)     # applied ONLY to eCR-targeted regions
```

where `efficacy ∈ [0,1]` is a **user assumption** exposed as a live slider (default 0.70).
Off-target and secondary/network effects are not modeled. The two **measured** endpoints
flank the modeled column so the viewer reads it as an interpolation, and a MODEL banner
states this. Untargeted background regions stay at the start state.

## Inputs

| Flag | Required | Meaning |
|---|---|---|
| `--bundle` | yes | eCR_navigator run bundle dir (`nominations.tsv` + `manifest.json`). Supplies target regions, `direction` (open/close), and transition labels. |
| `--atpm` | see below | Per-region measured accessibility TSV (see format). |
| `--targets` | no | BED of eCR-delivered regions. **Default: all nominations** are treated as targeted. |
| `--genes` | no | Gene BED (`chrom,start,end,name`) → heatmap row labels (nearest/overlapping gene). |
| `--highlight` | no | TSV `rank⇥label[⇥sublabel]` → per-locus "delivered construct" cards + heatmap marks. |
| `--annot` | no | TSV `rank⇥flag` → a small flag column in the heatmap (e.g. `3/3 routes`). |
| `--focused-ranks` | no | File of ranks to show in the heatmap (default: top by rank). |
| `--focused-max` | no | Max regions per direction in the heatmap (default 30). |
| `--title` | no | Page title (default derived from `manifest.transition.state_a/state_b`). |
| `--output` / `--data-out` | no | HTML path (default `accessibility_preview.html`) / optional JSON dump. |

### Accessibility source resolution (`--atpm` **overrides** the bundle)

1. `--atpm PATH` if given → use it;
2. else `<bundle>/accessibility.tsv` if present → use it;
3. else **error**.

> The bundle does not carry accessibility today (it holds `driver_score` + `direction`=ΔaTPM,
> not the absolute per-state values). Passing `--atpm` (the navigator's `atpm_union.tsv`) is the
> current path. **Planned, additive:** the navigator emits an optional `accessibility.tsv` into
> the run bundle (it already computes it), after which `--atpm` becomes optional. This is an
> optional 4th bundle part — it does **not** change the frozen three-part contract.

### `--atpm` / `accessibility.tsv` format

Tab-separated, header required. Two `atpm_*` columns; the suffix auto-labels the states:

```
chrom   start   end   atpm_MEF   atpm_mES
chr1    100     200   0.011      0.838
```

`atpm_MEF` / `atpm_mES` → the page labels the two measured columns **MEF** and **mES**.
Values are per-region accessibility (any comparable 0–1-ish scale).

## Output — `viz_data.json` (the payload the template renders)

```jsonc
{
  "title": "...", "subtitle": "...", "start_label": "MEF", "end_label": "mES",
  "focused":  [ {"rank","coord","gene","sign":"open|close","start","end","flag","hl"} ],
  "highlight":[ {"rank","label","sub","coord","sign","start","end"} ],
  "overview": [ {"start","end","targeted":true,"band":"to-open (targeted)|background (not targeted)|to-close (targeted)"} ],
  "n_union": <int>, "n_targets": <int>,
  "provenance": ["...","...","..."],
  "note": "Predicted = start + efficacy*(end-start) for TARGETED regions only ..."
}
```

The HTML is a single file with the template + `__DATA__` inlined (theme-aware, self-contained;
also publishes cleanly as an Artifact). JS keys mirror this schema 1:1.

## Usage

```bash
python preview.py \
  --bundle ../eCR_navigator/examples/run_bundle/<run> \
  --atpm   atpm_union.tsv \
  --genes  markers.bed  --highlight round1.tsv  --annot routes.tsv \
  --output accessibility_preview.html --data-out viz_data.json
```

Generic across transitions — nothing is hardcoded to MEF→mES; labels come from the aTPM
header and the bundle manifest.

## Extension points (future features — keep documented)

- **Real predicted track.** Swap the in-page direct-effect interpolation for a measured or
  model-predicted post-eCR track via a future `--predicted-track` TSV; the `viz_data` schema and
  template need no change (the model column just reads supplied values instead of interpolating).
- **Enrichment annotations** (gene-set tiers, experimental ATAC-GT concordance) are produced
  **outside** this stage and passed via `--annot` / `--highlight`, keeping the core stage lean.
