# CLAUDE.md — ECR_predictor

## Coding principles

- **No features beyond what was asked.** No speculative abstractions, configurability, or error handling for impossible scenarios.
- **Surgical changes.** Touch only what the task requires. Don't improve adjacent code, comments, or formatting. Match existing style.
- **Surface tradeoffs before coding.** If multiple interpretations exist, present them — don't pick silently. If something is unclear, ask.
- **If you notice unrelated dead code or issues, mention them — don't silently fix them.**
- **Stubs are stubs.** `af3.py` and `foldx.py` are intentional placeholders — don't flesh them out unless explicitly asked.
- **Scores are intentionally independent.** `motif_score` and `annotation_confidence` must never be merged into a single composite score without explicit instruction.

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
│   ├── prefetch.py   # pre-download JASPAR motifs to jaspar_cache/
│   ├── filter.py     # drop low-confidence hits (annotation + motif_score)
│   ├── fimo.py       # FIMO motif validation (requires MEME Suite)
│   ├── af3.py        # AlphaFold 3 structure prediction (stub — see TODOs)
│   └── foldx.py      # FoldX binding affinity estimation (stub — see TODOs)
├── jaspar_cache/     # .jaspar files stored here after prefetch (gitignored)
├── af3_outputs/      # AF3 JSON inputs + structure outputs (created by refine.py)
├── cli.py            # prediction entrypoint
├── refine.py         # refinement entrypoint (filter → FIMO → AF3 → FoldX)
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

## Refinement pipeline (refine.py)

Runs downstream validation on `cli.py` output. Use `--include-sequence` in `cli.py` to enable AF3.

```bash
python cli.py --sequence <SEQ> --species "Homo sapiens" --output results.tsv --include-sequence
python refine.py --input results.tsv --sequence <SEQ> --config config.yaml [--stop-after fimo]
```

| Stage | Flag | Status | Prerequisite |
|---|---|---|---|
| Filter | always on | done | — |
| FIMO | `--fimo-pvalue` | done | MEME Suite on PATH |
| AF3 | `--top-n-af3` | done | configure `config.yaml` |
| FoldX | automatic | stub | see `ecr_predictor/foldx.py` TODOs |

### config.yaml
`config.yaml` is gitignored (contains passwords/hosts). Copy from the template and edit:
```bash
cp config.example.yaml config.yaml
```
```yaml
af3:
  backend: hpcc   # local | hpcc | online
  hpcc:
    host: hpcc.example.edu
    user: wuyuhao
    remote_workdir: /scratch/wuyuhao/ecr_af3_jobs
    slurm_partition: a40-tmp
    slurm_qos: gpu
    slurm_module: alphafold/3_a40-tmp
    ...
```

### AF3 backends
- **local** — calls `run_alphafold3.sh input_dir output_dir json_file` via subprocess
- **hpcc** — scp JSON → SSH sbatch → poll squeue → scp CIF back. Requires key-based SSH auth (`ssh user@host echo ok`)
- **online** — Chai-1 API (https://chaidiscovery.com). Submits FASTA → polls → downloads CIF. Requires `ECR_CHAI_API_KEY` env var or `af3.online.api_key` in config.yaml.

AF3 JSON format: protein chain A + single-stranded DNA chain B. Job name = gene name. 5 model seeds.
Output CIF stored in `af3_outputs/<gene_name>/`. `af3_cif_path` column added to TSV.

### Filter logic
Drops rows where **both** are true: `annotation_confidence == 'low'` AND `motif_score < --min-motif-score` (default 0.0; NA counts as below threshold).

### FIMO
- Converts JASPAR cache → MEME format on the fly.
- Adds `fimo_pvalue` and `fimo_validated` columns. Best (lowest) p-value per motif is reported.

### FoldX
- RepairPDB → AnalyseComplex on each AF3 CIF. Output: `foldx_ddg_kcal_mol` (lower = stronger binding).
- Output parsing (`Interaction_*_AC.fxout`) marked TODO in `ecr_predictor/foldx.py`.

## Key implementation notes

- DBDs with no `jaspar_id` are **not dropped** — they appear with `motif_score = NA`.
- Species matching: exact first, then genus-level fallback (first word of species name). Flagged in `query_species_match` column.
- `pssm.calculate()` returns a scalar when sequence length == motif length; wrapped with `np.atleast_1d`.
- Sequences shorter than a motif return `motif_score = NA` (not an error).
- Motif fetches are parallelized with `ThreadPoolExecutor` (8 workers).
