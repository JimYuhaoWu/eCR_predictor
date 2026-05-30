# CLAUDE.md — ECR_predictor

## What this project does

CLI tool that predicts which DNA-Binding Domains (DBDs) from the `eCR_mod_lib` library are likely to bind a given DNA sequence and species. Returns a ranked TSV with two independent confidence scores.

## Repo layout

```
ECR_predictor/
├── ecr_predictor/
│   ├── query.py      # DBD lookup + species matching (exact → genus fallback)
│   ├── scan.py       # JASPAR PWM scoring via BioPython
│   ├── score.py      # validation_level → annotation_confidence label
│   ├── output.py     # build result table, write TSV
│   └── prefetch.py   # pre-download JASPAR motifs to jaspar_cache/
├── jaspar_cache/     # .jaspar files stored here after prefetch (gitignored)
├── cli.py            # argparse entry point
├── server_setup.sh   # one-time server setup (install + seed DB + prefetch motifs)
└── server_run.sh     # run a prediction on the server
```

## Sibling repo dependency

`eCR_mod_lib` (`scripts.schema.ModuleLibrary`) must be installed as an editable package alongside this repo:

```
parent_dir/
├── eCR_mod_lib/    ← pip install -e here first
└── eCR_predictor/  ← pip install -e here second
```

## Development workflow

- **Develop and test locally on PC** — run `cli.py` directly:
  ```bash
  python cli.py --sequence <SEQ> --species "Homo sapiens" --db path/to/module_library.db
  ```
- **Deploy to server** — `git pull` in both repos, then re-run `server_setup.sh` if the DB or motif cache needs updating.

## JASPAR motif fetch order

`scan.py` resolves motifs in this priority order:
1. `jaspar_cache/<id>.jaspar` — local file written by `prefetch.py` (fastest, no network)
2. BioPython `JASPAR2020` local DB — if the `jaspar2020` package is installed
3. JASPAR REST API (`https://jaspar.elixir.no/api/v1/`) — fallback, requires internet

On the server, run `server_setup.sh` once to populate the cache. Subsequent runs are fully offline.

## Confidence scores

Two columns, intentionally independent:

| Column | Source | Meaning |
|---|---|---|
| `motif_score` | JASPAR PWM log-odds, normalized −1 to 1 | Sequence-level evidence |
| `annotation_confidence` | Mapped from `validation_level` | Curation-level evidence |

`annotation_confidence` mapping:
- `high` → `screen-validated`, `ChIP-validated`, `structurally-resolved`
- `medium` → `motif-only`
- `low` → `predicted`

## Key implementation notes

- DBDs with no `jaspar_id` are **not dropped** — they appear with `motif_score = NA`.
- Species matching: exact first, then genus-level fallback (first word of species name). Flagged in `query_species_match` column.
- `pssm.calculate()` returns a scalar when sequence length == motif length; wrapped with `np.atleast_1d`.
- Sequences shorter than a motif return `motif_score = NA` (not an error).
- Motif fetches are parallelized with `ThreadPoolExecutor` (8 workers).
