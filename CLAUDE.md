# CLAUDE.md — ECR_predictor

## Coding principles

- **No features beyond what was asked.** No speculative abstractions, configurability, or error handling for impossible scenarios.
- **Surgical changes.** Touch only what the task requires. Don't improve adjacent code, comments, or formatting. Match existing style.
- **Surface tradeoffs before coding.** If multiple interpretations exist, present them — don't pick silently. If something is unclear, ask.
- **If you notice unrelated dead code or issues, mention them — don't silently fix them.**
- **Scores are intentionally independent.** `motif_score` and `annotation_confidence` must never be merged into a single composite score without explicit instruction.

## What this project does

Three-stage CLI pipeline:
1. **`cli.py`** — predicts which DBDs from `eCR_mod_lib` are likely to bind a query DNA sequence. Returns a ranked TSV with two independent confidence scores.
2. **`refine.py`** — filters, FIMO-validates, runs AF3 structure prediction, and estimates binding affinity with FoldX.
3. **`fuse.py`** — assembles DBD + linker + ED fusion candidates and screens them for developability liabilities (MHC-I junction immunogenicity, aggregation, stability) before wet-lab synthesis. Target modality: intracellular expression (viral vector / mRNA).

## Repo layout

```
ECR_predictor/
├── ecr_predictor/
│   ├── query.py          # DBD lookup + species matching (exact → genus fallback)
│   ├── scan.py           # JASPAR PWM scoring via BioPython (parallel fetch)
│   ├── score.py          # validation_level → annotation_confidence label
│   ├── output.py         # build result table, write TSV
│   ├── prefetch.py       # pre-download JASPAR motifs to jaspar_cache/
│   ├── filter.py         # drop low-confidence hits (annotation + motif_score)
│   ├── fimo.py           # FIMO motif validation (requires MEME Suite)
│   ├── af3.py            # AF3 structure prediction (local / hpcc / online backends)
│   ├── foldx.py          # FoldX binding affinity estimation
│   ├── fusion/           # Step 3 fusion-design developability gates
│   │   ├── backends.py       # shared local-CLI-or-API tool dispatch
│   │   ├── assemble.py       # build DBD+linker+ED candidates, track junctions
│   │   ├── junction.py       # junction-peptide enumeration + self-filtering
│   │   ├── immunogenicity.py # Gate 1: MHC-I neoepitope screen
│   │   ├── aggregation.py    # Gate 2: built-in AGGRESCAN scorer (a3v + freesasa) / CamSol via api
│   │   ├── stability.py      # Gate 3: N-end rule / degron / UbPred
│   │   └── score.py          # composite risk + Pareto ranking
│   └── config.py         # load/validate config.yaml
├── jaspar_cache/          # .jaspar files after prefetch (gitignored)
├── af3_outputs/           # AF3 JSON inputs, CIF outputs, run_log.json (gitignored)
│   └── jobs/
│       ├── <gene>.json       # AF3 input per job
│       ├── <gene>.slurm      # generated Slurm script (hpcc backend)
│       └── run_log.json      # job state log for resume on interruption
├── cli.py                 # Step 1 entrypoint
├── refine.py              # Step 2 entrypoint
├── fuse.py                # Step 3 entrypoint (fusion-design developability)
├── config.example.yaml    # config template — copy to config.yaml and edit
├── environment.yml        # shared conda environment (covers both repos)
├── server_setup.sh        # one-time server setup (install + seed DB + prefetch)
└── server_run.sh          # run a prediction on the server
```

## Environment & sibling repo dependency

Shared conda environment (`ecr`) covers both repos:

```bash
conda env create -f environment.yml   # once per machine
conda activate ecr
python -m pip install -e ../eCR_mod_lib   # must install mod_lib first
python -m pip install -e .
```

**Always use `python -m pip`, not bare `pip`.** On multi-Python servers the `pip`
on PATH is often a system/user pip bound to a different interpreter, so
`pip install -e ../eCR_mod_lib` silently installs into the wrong Python and the
env can't import the `scripts` package (`No module named 'scripts'`) even though
the install "succeeded". Verify with `python -m pip -V` vs `pip -V`. Because
`ecr_predictor` is picked up from the working directory when run from the repo
root, this surfaces first as a missing `scripts` (the sibling eCR_mod_lib
package), not a missing `ecr_predictor`. `check_fusion_env.sh` detects this.

Both repos must sit as siblings:
```
parent_dir/
├── eCR_mod_lib/
└── eCR_predictor/
```

Default DB path resolves to `../eCR_mod_lib/library/module_library.db`. Override with `--db`.

## Development workflow

```bash
# Step 1 — predict
python cli.py \
  --sequence ACAGGAAGTGACAGGAAGTGACAGGAAGTG \
  --species "Homo sapiens" \
  --output predictions.tsv \
  --include-sequence          # required for AF3 stage

# Step 2 — refine
python refine.py \
  --input predictions.tsv \
  --sequence ACAGGAAGTGACAGGAAGTGACAGGAAGTG \
  --config config.yaml \
  [--stop-after fimo]         # omit to run through AF3
```

Deploy to server: `git pull` in both repos, re-run `server_setup.sh` if DB or motif cache needs updating.

## JASPAR motif fetch order

`scan.py` resolves motifs in this priority order:
1. `jaspar_cache/<id>.jaspar` — local file (fastest, no network)
2. BioPython `JASPAR2020` local DB — if the `jaspar2020` package is installed
3. JASPAR REST API (`https://jaspar.elixir.no/api/v1/`) — fallback, requires internet

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

### Parameters

| Flag | Default | Description |
|---|---|---|
| `--input` | *(required)* | Predictor output TSV from `cli.py` |
| `--sequence` | *(required)* | Same DNA sequence used in `cli.py` |
| `--config` | `config.yaml` in repo root | Path to config.yaml |
| `--output` | `<input>_refined.tsv` | Output TSV path |
| `--min-motif-score` | `0.0` | Filter threshold for motif_score |
| `--fimo-pvalue` | `1e-4` | FIMO p-value cutoff for validation |
| `--top-n-af3` | `2` | Number of top hits sent to AF3 |
| `--stop-after` | *(none)* | Stop after `filter`, `fimo`, `af3`, or `foldx` |
| `--af3-output-dir` | `af3_outputs/` | Local dir for AF3 inputs/outputs |

### Stages

| Stage | Status | Prerequisite |
|---|---|---|
| Filter | done | — |
| FIMO | done | MEME Suite on PATH |
| AF3 | done | configure `config.yaml` |
| FoldX | done | `FOLDX_PATH` env var; biopython for CIF→PDB |

### Filter logic
Drops rows where **both** are true: `annotation_confidence == 'low'` AND `motif_score < --min-motif-score` (default 0.0; NA counts as below threshold).

### FIMO
- Converts JASPAR cache → MEME format on the fly (no extra files needed).
- Calls `fimo --text` on the query sequence.
- Adds `fimo_pvalue` and `fimo_validated` columns. Best (lowest) p-value per motif is reported.

### AF3 backends

Selected via `af3.backend` in `config.yaml`:

- **`local`** — calls `run_alphafold3.sh input_dir output_dir json_file` via `bash -c`. Handles `module load` via `af3.local.module_load`. Use when running `refine.py` interactively on the HPCC.
- **`hpcc`** — uploads JSON via SFTP, submits via `sbatch` over SSH, polls `squeue`, downloads CIF back. Supports two authentication methods (see below). Requires `af3.hpcc.auth_method` configured.
- **`online`** — Chai-1 API (https://chaidiscovery.com). Submits protein+DNA FASTA, polls, downloads CIF. Requires `ECR_CHAI_API_KEY` env var or `af3.online.api_key` in config.yaml.

AF3 JSON: protein chain A + double-stranded DNA — chain B = query (sense) strand, chain C = auto-generated reverse complement, folded as a base-paired duplex. For C2H2 zinc-finger DBDs (`tf_family` ∈ `zf-C2H2`, `ZBTB`), a ligand entity adds one `ZN` ion per UniProt-annotated finger (`zinc_finger_count`, chain ids `ZA`, `ZB`, …); no bond constraints. `version: 2`, `modelSeeds: [1]`. Job name = gene name (spaces → underscores). AF3 lowercases the job name for its output dir: `FLI1` → `fli1/fli1_model.cif`. Note: the Chai-1 (`online`) backend does not forward Zn ions — its FASTA path carries protein + DNA only.

### AF3 run log (resume on interruption)

State is persisted to `af3_outputs/jobs/run_log.json` after every transition:

```json
{
  "FLI1": {
    "job_id": "7888566",
    "status": "completed",
    "local_cif": "af3_outputs/fli1/fli1_model.cif",
    "submitted_at": "2026-06-01T08:00:00Z"
  }
}
```

Statuses: `submitted` → `running` → `completed` | `failed` | `cancelled` | `download_failed`

Resume behaviour on re-run:
- `completed` + CIF on disk → returned immediately, no SSH
- `submitted` / `running` + job_id → re-attaches and re-polls
- `failed` / `cancelled` / no entry → re-submits
- `download_failed` → re-polls then re-downloads

To force re-submission of a specific job, delete its entry from `run_log.json`.

### config.yaml

`config.yaml` is gitignored. Copy from the template:
```bash
cp config.example.yaml config.yaml
```

Key fields for the HPCC backend:
```yaml
af3:
  backend: hpcc
  hpcc:
    host: hpcc.example.edu
    port: 22
    user: your_username
    
    # Authentication: choose 'password' or 'totp'
    auth_method: password
    
    # For password-based auth (static password):
    ssh_password: ""          # or set ECR_HPCC_PASSWORD env var
    
    # For TOTP-based auth (time-based one-time passwords):
    totp_secret: ""           # or set ECR_HPCC_TOTP_SECRET env var
    
    remote_workdir: /scratch/your_username/ecr_af3_jobs
    slurm_partition: a40-tmp
    slurm_qos: gpu
    slurm_module: alphafold/3_a40-tmp
    poll_interval: 60
    timeout: 7200
```

**Authentication methods:**
- **`password`** — Static password from `ECR_HPCC_PASSWORD` env var or `af3.hpcc.ssh_password` in config.yaml.
- **`totp`** — Time-based one-time passwords (FreeOTP). Provide the shared secret via `ECR_HPCC_TOTP_SECRET` env var or `af3.hpcc.totp_secret` in config.yaml. The script generates codes automatically.

### FoldX
- Optionally trims low-confidence terminal loops from the protein (default threshold: pLDDT ≥ 70).
- CIF → PDB conversion (BioPython) → RepairPDB → AnalyseComplex.
- Output column: `foldx_ddg_kcal_mol` in kcal/mol (lower = stronger binding).
- Intermediate files written to `foldx_work/<gene>/` (persistent — RepairPDB is skipped on re-run if already done, as it takes ~4 min per structure).
- FoldX binary: set `FOLDX_PATH` env var (e.g. `export FOLDX_PATH=/path/to/foldx/foldx`).
- FoldX cannot resolve `~` in paths — all paths are fully resolved before use.

**Confidence-based trimming:**
- Uses sliding window (default window size = 3) over pLDDT scores (B-factors in AF3 CIF).
- Removes residues from both termini until reaching the confidence threshold.
- Set `confidence_threshold=None` in code to disable trimming.
- Configurable via `refine.py` parameters (future: add CLI flags if needed).

## Fusion-design stage (fuse.py — Step 3)

Separate entrypoint (not a `refine.py` stage) because the input model differs:
it operates on **fusion candidates** (DBD + linker + ED), not DBD–DNA hits.
DBDs come from a TSV (`gene_name` + `sequence_aa`); **EDs come from the library**
(`query.load_eds()`, `type='ED'`); linkers come from the `fusion:` config section.
Stages: `assemble → sequence → prune → structure → score` (`--stop-after <stage>`).

**Ordering is load-bearing:** the cheap SEQUENCE gates (immunogenicity, stability)
run on the whole library, then `prune` keeps the Pareto-optimal survivors
(capped by `--top-n-structure`), and only then does the STRUCTURE phase spend
HPCC/GPU (AF3 + FoldX) on survivors. Never fold the full DBD×linker×ED product.

**Target modality is intracellular expression** (viral vector / mRNA). This is
the other load-bearing assumption: the dominant immune liability is MHC-I
presentation of junction neoepitopes → CD8 killing of transduced cells.
Antibody/B-cell and serum-protease screening are deliberately omitted.

### Tool backend abstraction
Each external tool is configured under `fusion.tools.<tool>` with
`backend: local | api | disabled`, mirroring AF3. `fusion/backends.py` provides
the shared plumbing: `resolve_backend()`, `run_cli()` (local binary; `{key}`
tokens in args map to temp input files), and `submit_and_poll()` (generic
API submit/poll with caller-supplied id/status/result callbacks). Gate modules
build tool-specific commands/payloads and parse tool-specific output.

### Gates
- **Function retention** — not a separate gate; the structure phase
  (`fusion/structure.py`) reuses `af3.run_af3_prediction` + `foldx.run_foldx_affinity`
  on the *fused* construct (distinct from the bare DBD–DNA complex `refine.py`
  folds). Reports `fusion_ddg_kcal_mol` and, when the input TSV carries a per-DBD
  `foldx_ddg_kcal_mol` baseline, `function_delta_ddg`. Zinc-finger handling is
  preserved by propagating the DBD's `tf_family`/`zinc_finger_count` into the AF3
  input. Reported, not a Pareto axis (absolute ΔΔG isn't comparable across DBDs).
- **Gate 1 immunogenicity** — `junction_peptides()` enumerates k-mers (8–11) that
  straddle a domain boundary; `partition_self()` drops any matching the
  `self_proteome` FASTA (exact substring). Remaining peptides go to NetMHCpan
  (binding by %rank; NetCTLpan — proteasome+TAP+MHC — is discontinued/legacy) across an
  HLA-I panel. Flag by `%rank ≤ rank_threshold` (default 2.0); report epitope
  **density** (flagged/tested), not single hits.
- **Gate 2 aggregation** — built-in AGGRESCAN-style scorer (`local` backend, no
  external binary): published a3v scale + length-adaptive window, weighted by
  per-residue relative SASA (`freesasa`) from the AF3 survivor structure; pure
  sequence fallback if no structure/freesasa. CamSol is web-only (`api` backend
  only). Flags hotspots within `junction_window` residues of a boundary.
- **Gate 3 stability** — N-end-rule + degron regex scan run with NO external tool
  (always available); UbPred optional via backend. NOTE: degradation feeds Gate-1
  presentation — the two axes conflict; the pipeline surfaces both, doesn't resolve.
- **Score** — `score.py` reports each liability axis and marks the **Pareto-optimal**
  set; `risk_score` (normalized axis sum) is a convenience sort only.

### Parsers are version-sensitive
The NetMHCpan/NetCTLpan stdout parsers target specific tool versions (NetMHCpan
4.2, NetCTLpan 1.1) and locate columns by header name. They are the most likely
thing to break against a different install — validate output parsing before
trusting scores. (Gate 2 aggregation no longer parses external stdout — it's the
built-in a3v + freesasa scorer.) The tool-independent logic (assembly, junction
enumeration, self-filtering, N-end rule, degron scan, Pareto) is unit-testable
without any external binary.

## Key implementation notes

- DBDs with no `jaspar_id` are **not dropped** — they appear with `motif_score = NA`.
- Species matching: exact first, then genus-level fallback (first word of species name). Flagged in `query_species_match` column.
- `pssm.calculate()` returns a scalar when sequence length == motif length; wrapped with `np.atleast_1d`.
- Sequences shorter than a motif return `motif_score = NA` (not an error).
- Motif fetches are parallelized with `ThreadPoolExecutor` (8 workers).
- Python 3.8 compatibility required (server constraint) — `with_stem()` unavailable, use `with_name()`.
- `paramiko` and `pyyaml` are required for the refinement pipeline (included in `environment.yml`).
