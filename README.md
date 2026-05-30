# ECR_predictor

Given a DNA sequence and a species, predict which DNA-Binding Domains (DBDs) from the [`eCR_mod_lib`](https://github.com/JimYuhaoWu/eCR_mod_lib) library are likely to bind that sequence. Returns a ranked TSV table with two independent confidence scores.

---

## Installation

**1. Install the module library as a local editable package**

```bash
pip install -e ../module_library
```

> The library lives at `D:\module_library` (or wherever you cloned it). It is imported as `ecr_mod_lib`.

**2. Install ECR_predictor**

```bash
pip install -e .
```

---

## Usage

```bash
python cli.py --sequence <DNA> --species "<species>" [--output results.tsv] [--db path/to/module_library.db]
```

| Argument | Required | Description |
|---|---|---|
| `--sequence` | Yes | DNA sequence to scan (A/C/G/T/N only) |
| `--species` | Yes | Query species, e.g. `"Homo sapiens"` |
| `--output` | No | Output TSV path; defaults to stdout |
| `--db` | No | Path to `module_library.db`; auto-detected if omitted |

**Example**

```bash
python cli.py \
  --sequence ATCGATCGATCG \
  --species "Homo sapiens" \
  --output predictions.tsv
```

---

## Output columns

| Column | Description |
|---|---|
| `gene_name` | DBD gene symbol |
| `species` | Species of the library record |
| `query_species_match` | `exact` — species matched directly; `other` — genus-level fallback |
| `tf_family` | Transcription factor family / subtype |
| `validation_level` | Raw validation level from the library |
| `motif_score` | Normalized PWM log-odds score (−1 to 1), or `NA` if no JASPAR motif |
| `annotation_confidence` | `high` / `medium` / `low` derived from `validation_level` |
| `jaspar_id` | JASPAR motif ID if available |

Results are sorted: exact species matches first, then by `motif_score` descending (`NA` last).

---

## Confidence scores

The two scores are **intentionally independent**:

- **`motif_score`** — sequence-level evidence. Computed from JASPAR PWM scanning (BioPython). `NA` if the DBD has no associated JASPAR motif. Fetched from the local BioPython JASPAR2020 DB if installed, otherwise via the [JASPAR REST API](https://jaspar.elixir.no).

- **`annotation_confidence`** — curation-level evidence mapped from `validation_level`:
  - `high` → `screen-validated`, `ChIP-validated`, `structurally-resolved`
  - `medium` → `motif-only`
  - `low` → `predicted`

A high `motif_score` + `low` annotation confidence is a very different hit from a high `motif_score` + `high` annotation confidence.

---

## Species matching

Exact species match is tried first. If no records match, genus-level fallback is used (first word of the species name). Non-exact hits are flagged with `query_species_match = other`.

---

## Project structure

```
ECR_predictor/
├── ecr_predictor/
│   ├── query.py    # DBD lookup + species matching
│   ├── scan.py     # JASPAR PWM scanning
│   ├── score.py    # validation_level → annotation_confidence
│   └── output.py   # table formatting and TSV output
├── cli.py          # argparse entry point
└── pyproject.toml
```

---

## Dependencies

- [eCR_mod_lib](https://github.com/JimYuhaoWu/eCR_mod_lib) (local editable install)
- [biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
- [requests](https://requests.readthedocs.io/)
