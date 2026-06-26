# ECR_predictor

Given a DNA sequence and a species, predict which DNA-Binding Domains (DBDs) from the [`eCR_mod_lib`](https://github.com/JimYuhaoWu/eCR_mod_lib) library are likely to bind that sequence. Optionally refine top hits through FIMO motif validation, AlphaFold 3 structure prediction, and FoldX binding affinity estimation.

---

## Installation

Both repos must sit in the same parent directory:

```
parent_dir/
├── eCR_mod_lib/
└── eCR_predictor/
```

```bash
# Clone both repos
git clone https://github.com/JimYuhaoWu/eCR_mod_lib.git
git clone https://github.com/JimYuhaoWu/eCR_predictor.git

# Create and activate the shared conda environment (once per machine)
conda env create -f eCR_predictor/environment.yml
conda activate ecr

# Install both packages in editable mode
pip install -e eCR_mod_lib
pip install -e eCR_predictor
```

The `ecr` conda environment covers all dependencies for both projects, including MEME Suite (`fimo`) for the refinement pipeline. If the environment already exists (e.g. created from `eCR_mod_lib`), skip the `conda env create` step.

**One-time server setup** (seeds the DB and pre-fetches JASPAR motifs for offline use):
```bash
bash server_setup.sh
```

---

## Step 1 — Predict

```bash
python cli.py \
  --sequence ACAGGAAGTGACAGGAAGTGACAGGAAGTG \
  --species "Homo sapiens" \
  --output predictions.tsv \
  --include-sequence
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--sequence` | Yes | — | DNA sequence (A/C/G/T/N only) |
| `--species` | Yes | — | Query species, e.g. `"Homo sapiens"` |
| `--output` | No | stdout | Output TSV path |
| `--db` | No | auto-detected | Path to `module_library.db` |
| `--include-sequence` | No | off | Add `sequence_aa` column (required for AF3) |

### Output columns

| Column | Description |
|---|---|
| `gene_name` | DBD gene symbol |
| `species` | Species of the library record |
| `query_species_match` | `exact` — direct match; `other` — genus-level fallback |
| `tf_family` | Transcription factor family / subtype |
| `validation_level` | Raw curation level from the library |
| `motif_score` | Normalized PWM log-odds score (−1 to 1); `NA` if no JASPAR motif |
| `annotation_confidence` | `high` / `medium` / `low` derived from `validation_level` |
| `jaspar_id` | JASPAR motif ID if available |
| `sequence_aa` | DBD amino acid sequence (only with `--include-sequence`) |
| `zinc_finger_count` | UniProt-annotated zinc-finger count for C2H2 DBDs (only with `--include-sequence`); drives Zn ions in AF3 |

Results are sorted: exact species matches first, then by `motif_score` descending (`NA` last).

---

## Step 2 — Refine

Takes the `predictions.tsv` from Step 1 and runs four stages in sequence:

```
[1] Filter → [2] FIMO → [3] AF3 → [4] FoldX
```

```bash
python refine.py \
  --input predictions.tsv \
  --sequence ACAGGAAGTGACAGGAAGTGACAGGAAGTG \
  --config config.yaml
```

Use `--stop-after fimo` to stop before AF3 (e.g. while setting up the HPCC connection).

### Parameters

| Flag | Default | Description |
|---|---|---|
| `--input` | *(required)* | Predictor output TSV from Step 1 |
| `--sequence` | *(required)* | Same DNA sequence used in Step 1 |
| `--config` | `config.yaml` | Path to config file |
| `--output` | `<input>_refined.tsv` | Output TSV path |
| `--min-motif-score` | `0.0` | Filter: drop hits below this motif_score |
| `--fimo-pvalue` | `1e-4` | FIMO p-value cutoff for validation |
| `--top-n-af3` | `2` | Number of top hits submitted to AF3 |
| `--stop-after` | *(none)* | Stop after `filter`, `fimo`, `af3`, or `foldx` |
| `--af3-output-dir` | `af3_outputs/` | Local dir for AF3 inputs and outputs |

### Stage 1 — Filter

Drops hits where **both** conditions are true:
- `annotation_confidence == 'low'`
- `motif_score < --min-motif-score` (NA counts as below threshold)

### Stage 2 — FIMO

Validates remaining motif hits using [MEME Suite](https://meme-suite.org/meme/doc/install.html) `fimo`. Converts JASPAR motifs to MEME format on the fly — no extra files needed.

Adds columns: `fimo_pvalue`, `fimo_validated`.

### Stage 3 — AlphaFold 3 (AF3)

Predicts DBD–DNA complex structures for the top-N FIMO-validated hits. Three backends — configure in `config.yaml`:

| Backend | When to use |
|---|---|
| `local` | Running `refine.py` interactively on the HPCC (AF3 installed there) |
| `hpcc` | Submitting from a separate Linux server to the HPCC via SSH |
| `online` | [Chai-1 API](https://chaidiscovery.com) — sign up for an API key |

DNA is modelled as a double-stranded duplex (query strand + auto-generated reverse complement). For C2H2-type zinc-finger DBDs (`tf_family` ∈ `zf-C2H2`, `ZBTB`), one Zn²⁺ ion is added per UniProt-annotated finger (`zinc_finger_count`) to improve accuracy.

Adds column: `af3_cif_path`.

**Resume on interruption:** job state is persisted to `af3_outputs/jobs/run_log.json`. Re-running `refine.py` automatically resumes — re-polling running jobs, re-downloading failed downloads, skipping completed ones. To force re-submission of a specific job, delete its entry from the log file.

### Stage 4 — FoldX

Estimates binding affinity from AF3 structures. By default, trims low-confidence terminal loops before analysis to improve prediction accuracy.

**Trimming:** Removes residues from protein termini with pLDDT < 70 (configurable). Uses a sliding window strategy to identify high-confidence regions. This typically removes 5–45% of residues at the termini, depending on protein length and prediction quality.

**FoldX workflow:** Trimmed CIF → PDB conversion → RepairPDB → AnalyseComplex. Requires FoldX 5 installed on the server:

```bash
export FOLDX_PATH=/path/to/foldx/foldx
```

Adds column: `foldx_ddg_kcal_mol` in kcal/mol (lower = stronger predicted binding). Intermediate files are written to `foldx_work/<gene>/` and reused on re-runs (RepairPDB takes ~4 min per structure). If `FOLDX_PATH` is not set, this stage is skipped with a warning.

---

## Configuration (config.yaml)

`config.yaml` is gitignored (it contains credentials). Copy the template and edit:

```bash
cp config.example.yaml config.yaml
```

**HPCC backend** — submits Slurm jobs from the Linux server to the HPCC via SSH:
```yaml
af3:
  backend: hpcc
  hpcc:
    host: hpcc.example.edu
    port: 22
    user: your_username
    
    # Authentication: choose 'password' or 'totp'
    auth_method: password
    
    # For password-based auth:
    ssh_password: ""          # or set ECR_HPCC_PASSWORD env var
    
    # For TOTP-based auth (FreeOTP):
    totp_secret: ""           # or set ECR_HPCC_TOTP_SECRET env var
    
    remote_workdir: /scratch/your_username/ecr_af3_jobs
    slurm_partition: a40-tmp
    slurm_qos: gpu
    slurm_module: alphafold/3_a40-tmp
    poll_interval: 60         # seconds between squeue polls
    timeout: 7200             # max wait per job (0 = wait forever)
```

**Local backend** — run directly on the HPCC:
```yaml
af3:
  backend: local
  local:
    af3_script: run_alphafold3.sh
    module_load: alphafold/3_a40-tmp
```

**Online backend** — Chai-1 API:
```yaml
af3:
  backend: online
  online:
    api_key: ""               # or set ECR_CHAI_API_KEY env var
```

---

## Project structure

```
ECR_predictor/
├── ecr_predictor/
│   ├── query.py          # DBD lookup + species matching
│   ├── scan.py           # JASPAR PWM scanning (parallel fetch)
│   ├── score.py          # validation_level → annotation_confidence
│   ├── output.py         # table formatting and TSV output
│   ├── prefetch.py       # pre-download JASPAR motifs to local cache
│   ├── filter.py         # confidence + score filtering
│   ├── fimo.py           # FIMO validation (JASPAR → MEME format + fimo call)
│   ├── af3.py            # AF3 prediction (local / hpcc / online backends)
│   ├── foldx.py          # FoldX affinity estimation
│   └── config.py         # load/validate config.yaml
├── jaspar_cache/          # populated by server_setup.sh (gitignored)
├── af3_outputs/           # AF3 inputs, CIF outputs, run_log.json (gitignored)
├── cli.py                 # Step 1 entrypoint
├── refine.py              # Step 2 entrypoint
├── config.example.yaml    # config template
├── environment.yml        # shared conda environment
├── server_setup.sh        # one-time server setup
└── server_run.sh          # run a prediction on the server
```

---

## Dependencies

All managed via `environment.yml`:

- [eCR_mod_lib](https://github.com/JimYuhaoWu/eCR_mod_lib) — sibling editable install
- [biopython](https://biopython.org/) — JASPAR PWM scoring
- [pandas](https://pandas.pydata.org/)
- [requests](https://requests.readthedocs.io/)
- [paramiko](https://www.paramiko.org/) — SSH/SFTP for HPCC backend
- [pyyaml](https://pyyaml.org/) — config file parsing
- [MEME Suite](https://meme-suite.org/) — FIMO stage (external tool, included in conda env)
- [FoldX](https://foldxsuite.crg.eu/) — affinity stage (external tool, free academic licence)
