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

## Step 3 — Fuse (fusion-design developability screen)

Once you've selected DBDs (Steps 1–2), `fuse.py` assembles **fusion candidates** (DBD + linker + effector domain) and screens each for liabilities that would sink a wet-lab synthesis. The **target modality is intracellular expression** (viral vector / mRNA), so the dominant immune risk is presentation of *junction-spanning neoepitopes* on **MHC class I** → CD8 killing of transduced cells. Antibody / serum-protease routes are not screened (the product isn't circulating).

```bash
python fuse.py \
  --dbd-input predictions.tsv \
  --sequence ACAGGAAGTGACAGGAAGTGACAGGAAGTG \
  --config config.yaml \
  --output fusion_candidates.tsv
```

DBDs are read from the `gene_name` + `sequence_aa` columns of the Step-1/2 output (use `--dbd FLI1,ETV6` to subset). **Effector domains are read from the eCR_mod_lib library** (`type='ED'`; use `--ed VP64,KRAB` to subset). The linker library and tool backends come from the `fusion:` section of `config.yaml`. `--sequence` is the DNA target the DBD must still bind (only needed for the structure phase).

### Pipeline

The pipeline runs **cheap sequence-based gates first**, prunes to the best candidates, and only then spends HPCC/GPU on the structure phase — so AF3 folding never runs on the full combinatorial library:

```
[1] Assemble → [2] Sequence gates → [3] Prune → [4] Structure → [5] Score
                (immuno + stability)  (Pareto)   (AF3 + FoldX + aggregation)
```

| Stage | Tool(s) | What it checks | Cost |
|---|---|---|---|
| Assemble | — | DBDs × linkers × EDs; tracks each domain junction | — |
| Immunogenicity (Gate 1) | NetMHCpan (NetCTLpan deprecated) | MHC-I binders among junction-spanning peptides | sequence |
| Stability (Gate 3) | N-end rule, degron scan, *UbPred* | proteasomal degradation liabilities | sequence |
| Prune | — | Pareto over the sequence axes → survivors (`--top-n-structure`) | — |
| Function + Aggregation | AF3, FoldX, AGGRESCAN3D / CamSol | DBD–DNA binding retained; junction aggregation | **HPCC/GPU** |
| Score | — | per-axis liabilities + Pareto-optimal flag | — |

**Function retention** is not a separate gate — it reuses the AF3 + FoldX stages on the *fused* construct (a different molecule from the bare DBD–DNA complex that `refine.py` folds). The fusion's FoldX ΔΔG, compared against a per-DBD baseline ΔΔG if present in the input TSV (`function_delta_ddg`), tells you whether fusing perturbed binding. A low-liability but non-functional fusion is useless, so this runs on the survivors before you trust the ranking.

Use `--stop-after prune` to get the sequence-only ranking without any HPCC work.

### Tool backends (local CLI or API)

Every external tool is invoked through a **local-CLI-or-API backend**, selected per tool in `config.yaml` exactly like the AF3 stage:

```yaml
fusion:
  tools:
    netmhcpan:
      backend: local        # local | api | disabled
      local: { command: netMHCpan }
      api:   { url: "", api_key_env: ECR_NETMHCPAN_API_KEY }
```

A gate whose tool is `disabled` (or whose binary is missing) is skipped gracefully and its axis is reported as `NA`. The **N-end-rule and degron scans in Gate 3 need no external tool** and always run.

### Installing the fusion tools

None of these ship with the conda environment — install only the ones for the gates you intend to run. The structure phase additionally needs **FoldX** (same binary as Stage 4) and AF3 access.

| Tool | Gate | Licence | Install |
|---|---|---|---|
| **AGGRESCAN3D** | Gate 2 (aggregation) | open | `pip install aggrescan3d freesasa` |
| **NetMHCpan 4.2** | Gate 1 (MHC-I) | free academic | [DTU download](https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/) → tarball |
| **NetCTLpan 1.1** | Gate 1 (deprecated) | free academic | discontinued by DTU — superseded by NetMHCpan; legacy installs only |
| **FoldX 5** | structure phase | free academic | [foldxsuite.crg.eu](https://foldxsuite.crg.eu/) → tarball |
| **CamSol** | Gate 2 (alt) | web server | no CLI — use AGGRESCAN3D, or the `api` backend |
| **UbPred** | Gate 3 (optional) | web server | no CLI — Gate 3 runs without it |

#### Local install (helper script)

`install_fusion_tools.sh` automates the whole thing. It `pip`-installs AGGRESCAN3D, and for the licence-gated tools (**NetMHCpan 4.2** and **FoldX**) you first **register and download the tarballs** from the links above and drop them into a `vendor/` directory — the script then extracts them, patches the NetMHCpan tcsh wrapper (`NMHOME` + `TMPDIR`), fetches NetMHCpan's data files, and symlinks the binaries into one bin dir. It's idempotent and skips any tool whose tarball isn't present, so you can re-run it as you obtain each licence. (NetCTLpan is discontinued by DTU; if you still have a legacy `netCTLpan-*` tarball in `vendor/`, the script will configure it too, but new setups should use NetMHCpan.)

```bash
conda activate ecr
mkdir vendor                       # drop the downloaded *.tar.gz files here
bash install_fusion_tools.sh       # installs into ~/opt/ecr_tools by default

# override locations if you like:
VENDOR_DIR=/data/tarballs INSTALL_DIR=/opt/ecr_tools bash install_fusion_tools.sh
```

Then follow the printed summary — add the bin dir to `PATH`, export `FOLDX_PATH`, and enable the installed tools in `config.yaml`:

```bash
export PATH="$HOME/opt/ecr_tools/bin:$PATH"
export FOLDX_PATH="$HOME/opt/ecr_tools/bin/foldx"
```

```yaml
fusion:
  tools:
    netmhcpan:   { backend: local, local: { command: netMHCpan } }   # one MHC-I tool (NetMHCpan 4.2)
    aggrescan3d: { backend: local, local: { command: aggrescan3d } }  # one aggregation tool
    # leave camsol / ubpred as 'disabled'
```

#### Using API backends instead

If you'd rather not install the binaries, set a tool's `backend: api` and point it at an HTTP endpoint — the same submit → poll → fetch flow the AF3 `online` backend uses. Supply the URL in `config.yaml` and the key via the named env var (never commit a real key):

```yaml
fusion:
  tools:
    netmhcpan:
      backend: api
      api:
        url: https://your-mhc-service.example/api/v1/predict
        api_key_env: ECR_NETMHCPAN_API_KEY   # read from this env var
        poll_interval: 10                     # seconds between status polls
        timeout: 600                          # give up after this many seconds
    aggrescan3d:
      backend: api
      api:
        url: https://your-a3d-service.example/api/v1/predict
        api_key_env: ECR_A3D_API_KEY
```

```bash
export ECR_NETMHCPAN_API_KEY='your_key'
export ECR_A3D_API_KEY='your_key'
```

The client POSTs the job payload to `url`, then polls `url/<job_id>` until the status reaches a done/fail state. The expected response shapes are: for MHC-I tools, `{ "peptide": rank, ... }` or a list of `{peptide, rank}` records; for aggregation tools, a per-residue `[{index, score}, ...]` list. DTU and CamSol do not publish such a REST API themselves — the `api` backend is for a service you host or wrap. A key set in the env var always overrides one written in `config.yaml`.

### Key design points

- **Sequence-first ordering.** Cheap gates prune the library before any AF3/FoldX run, so the expensive structure phase only ever sees `--top-n-structure` survivors.
- **Self-tolerance filtering.** Because both domains are endogenous, only junction-spanning peptides are potential neoepitopes — and any that occur verbatim in the human proteome are still self. Provide `fusion.self_proteome` (a proteome FASTA) to subtract them before the MHC-I scan.
- **%rank, not raw IC₅₀.** Binders are flagged by `%rank ≤ rank_threshold` (default 2.0) across an HLA-I panel; the gate reports **epitope density** (flagged / tested), not single hits.
- **Pareto over composite.** The output marks the Pareto-optimal set across the liability axes rather than collapsing them into one opaque score (a `risk_score` is also provided for quick sorting).
- **Degradation ↔ presentation tension.** Proteasomal degradation (Gate 3) is what *generates* MHC-I peptides (Gate 1); the pipeline surfaces both and does not auto-resolve the trade-off.

### Output columns

`candidate`, `dbd`, `linker`, `ed`, `length`, `immuno_density`, `immuno_flagged`, `immuno_min_rank`, `aggregation_risk`, `stability_risk`, `fusion_ddg_kcal_mol`, `function_delta_ddg`, `af3_cif_path`, `risk_score`, `pareto_optimal`.

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
