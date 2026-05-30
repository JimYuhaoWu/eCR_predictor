# ECR_predictor

Given a DNA sequence and a species, predict which DNA-Binding Domains (DBDs) from the [`eCR_mod_lib`](https://github.com/JimYuhaoWu/eCR_mod_lib) library are likely to bind that sequence. Returns a ranked TSV table with two independent confidence scores.

---

## Installation

Both repos must sit in the same parent directory:

```
parent_dir/
├── eCR_mod_lib/
└── eCR_predictor/
```

```bash
pip install -e eCR_mod_lib
pip install -e eCR_predictor
```

---

## Server setup (one-time)

```bash
bash eCR_predictor/server_setup.sh
```

This will:
1. Install both packages
2. Seed `module_library.db` (skipped if the DB already exists)
3. Pre-fetch all JASPAR motifs referenced in the library into `eCR_predictor/jaspar_cache/`

After this, predictions run **fully offline** — no network calls needed.

---

## Usage

**On the server:**
```bash
bash server_run.sh <sequence> "<species>" [output.tsv]
```

**Directly:**
```bash
python cli.py --sequence <DNA> --species "<species>" [--output results.tsv] [--db path/to/module_library.db]
```

| Argument | Required | Description |
|---|---|---|
| `--sequence` | Yes | DNA sequence to scan (A/C/G/T/N only) |
| `--species` | Yes | Query species, e.g. `"Homo sapiens"` |
| `--output` | No | Output TSV path; defaults to stdout |
| `--db` | No | Path to `module_library.db`; auto-detected if omitted |

**Example:**
```bash
python cli.py \
  --sequence ACAGGAAGTGACAGGAAGTGACAGGAAGTG \
  --species "Homo sapiens" \
  --output predictions.tsv
```

---

## Output columns

| Column | Description |
|---|---|
| `gene_name` | DBD gene symbol |
| `species` | Species of the library record |
| `query_species_match` | `exact` — direct match; `other` — genus-level fallback |
| `tf_family` | Transcription factor family / subtype |
| `validation_level` | Raw validation level from the library |
| `motif_score` | Normalized PWM log-odds score (−1 to 1), or `NA` if no JASPAR motif |
| `annotation_confidence` | `high` / `medium` / `low` derived from `validation_level` |
| `jaspar_id` | JASPAR motif ID if available |

Results are sorted: exact species matches first, then by `motif_score` descending (`NA` last).

---

## Confidence scores

The two scores are **intentionally independent**:

**`motif_score`** — sequence-level evidence from JASPAR PWM scanning.
- Computed as max log-odds score normalized by theoretical maximum, range −1 to 1.
- `NA` if the DBD has no associated JASPAR motif.

**`annotation_confidence`** — curation-level evidence:
| Value | `validation_level` sources |
|---|---|
| `high` | `screen-validated`, `ChIP-validated`, `structurally-resolved` |
| `medium` | `motif-only` |
| `low` | `predicted` |

A high `motif_score` + `low` annotation confidence is a very different hit from a high `motif_score` + `high` annotation confidence.

---

## JASPAR motif resolution

Motifs are resolved in this order (fastest first):
1. **Local cache** — `jaspar_cache/<id>.jaspar` files written by `server_setup.sh`
2. **BioPython JASPAR2020 DB** — if the `jaspar2020` package is installed
3. **JASPAR REST API** — `https://jaspar.elixir.no/api/v1/` (requires internet)

Run `server_setup.sh` once on the server to populate the cache for fully offline use.

To refresh the cache manually:
```bash
python -m ecr_predictor.prefetch --db path/to/module_library.db --cache-dir jaspar_cache/
```

---

## Species matching

Exact species match is tried first. If no records match, genus-level fallback is used (first word of the species name). Non-exact hits are flagged with `query_species_match = other`.

---

## Project structure

```
ECR_predictor/
├── ecr_predictor/
│   ├── query.py      # DBD lookup + species matching
│   ├── scan.py       # JASPAR PWM scanning (parallel fetch)
│   ├── score.py      # validation_level → annotation_confidence
│   ├── output.py     # table formatting and TSV output
│   └── prefetch.py   # pre-download all JASPAR motifs to local cache
├── jaspar_cache/     # populated by server_setup.sh (gitignored)
├── cli.py            # argparse entry point
├── server_setup.sh   # one-time server setup
└── server_run.sh     # run a prediction on the server
```

---

## Dependencies

- [eCR_mod_lib](https://github.com/JimYuhaoWu/eCR_mod_lib) (sibling editable install)
- [biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
- [requests](https://requests.readthedocs.io/)
- numpy (via biopython/pandas)
