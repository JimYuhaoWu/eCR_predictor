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
python -m pip install -e eCR_mod_lib
python -m pip install -e eCR_predictor
```

The `ecr` conda environment covers all dependencies for both projects, including MEME Suite (`fimo`) for the refinement pipeline. If the environment already exists (e.g. created from `eCR_mod_lib`), skip the `conda env create` step.

> **Use `python -m pip`, not bare `pip`.** On servers with multiple Python installs, `pip` on PATH may belong to a different interpreter than the `ecr` env, so a bare `pip install -e eCR_mod_lib` installs into the wrong Python — the install "succeeds" but the env then fails with `No module named 'scripts'`. Confirm they match with `python -m pip -V` vs `pip -V`. Run `bash eCR_predictor/check_fusion_env.sh` to verify the whole environment.

**One-time server setup** (seeds the DB and pre-fetches JASPAR motifs for offline use):
```bash
bash server_setup.sh
```

---

## Step 0 — Target (nominated region → binding site)

`cli.py` takes a **short** DNA site (~30 bp). [`eCR_navigator`](https://github.com/JimYuhaoWu/eCR_navigator) nominates **regions** (hundreds of bp). `target.py` is the bridge: it turns a navigator **run bundle** into a table of AF3-safe sites, one per library DBD that can bind inside a nominated region, ranked by on-target strength × specificity. This step is **optional** — if you already have a target site, hand it to Step 1 directly.

```bash
python target.py \
  --bundle    ../eCR_navigator/examples/run_bundle/in_gse299923 \
  --genome    /path/to/hg38.fa \
  --offtarget offtarget.tsv \
  --top-n     50 \
  --output    sites.tsv \
  --gaps-out  gaps.tsv
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--bundle` | Yes | — | eCR_navigator run-bundle dir (`manifest.json` + `nominations.tsv`) |
| `--genome` | Yes | — | Genome FASTA; assembly **must** match `manifest.transition.assembly` |
| `--offtarget` | No | — | `offtarget.py` output TSV, joined on `jaspar_id` for the specificity axis. Omit to rank on on-target strength alone |
| `--output` | No | `sites.tsv` | Ranked target-site table |
| `--gaps-out` | No | — | Library-gap report (regions no library DBD can target) |
| `--top-n` | No | `50` | Design budget: emit sites for the top-N regions **by rank** |
| `--site-length` | No | `30` | AF3 duplex length in bp |
| `--direction-min` | No | `0.05` | `\|direction\|` below this is labelled `ed_class=ambiguous` (labels only, never drops) |
| `--species` | No | from manifest | Override the library species filter |

The output `sequence` column feeds straight into Step 1:
```bash
python cli.py --sequence <sequence-from-sites.tsv> --species "Homo sapiens" --include-sequence
```

**How sites are chosen (library-first).** The library is small enough to scan exhaustively, so `target.py` takes the best PWM hit per library DBD inside each nominated region, then filters by off-target burden — not the reverse. Both strands are scanned (a target is folded as a duplex), and `strand` records which one the motif reads on.

**Four load-bearing rules:**
- **Order by `rank`, never `nomination_score`.** The navigator's score is a percentile, so a top-1% band is ~0.99–1.0 by construction and carries no within-band prioritization. `--top-n` cuts on `rank`.
- **`direction` picks the effector domain, not the DBD** (`+` → activator e.g. VP64; `−` → repressor e.g. KRAB), and flows through to `fuse.py`. Validated trust cutoff: `|direction| ≥ 0.05` (eCR_navigator Claim 2B); below it, `ed_class=ambiguous`.
- **Refusal is a valid answer.** A Gate-1-refused bundle (e.g. an anti-informative transition) yields an empty `sites.tsv` and exit 0 — not an error.
- **The gap report is the expansion signal.** Regions where *no* library DBD clears the motif threshold are undesignable today; `gaps.tsv` lists them with the best near-miss DBD — the direct input to `eCR_mod_lib` expansion. With only 9 of 38 library DBDs carrying a motif (8 human, 1 mouse), expect this to be the larger output.

### Output columns (`sites.tsv`)

| Column | Meaning |
|---|---|
| `site_rank` | Rank by `design_score` (on-target × specificity) |
| `region_id` | `chrom:start-end` of the source nomination |
| `nomination_rank` | `rank` from `nominations.tsv` (the orderable column) |
| `chrom`, `site_start`, `site_end` | Genomic coordinates of the emitted site (BED, + strand) |
| `sequence` | The site → `cli.py --sequence`. Always the + strand reference at `[site_start, site_end)` |
| `strand` | Strand the DBD's motif reads on (`+` / `−`) |
| `motif_start`, `motif_end` | The PWM hit inside the site |
| `gene_name`, `module_id`, `jaspar_id` | Library DBD that motivated the site |
| `motif_score` | On-target PWM strength (−1 to 1) |
| `n_offtarget_sites`, `offtarget_effect_score` | Per-DBD burden, joined from `--offtarget` |
| `specificity` | Library-relative rescaling of the burden (0–1) |
| `design_score` | `motif_score × specificity` — convenience sort only |
| `pareto_optimal` | On the (on-target ↑, burden ↓) Pareto front |
| `direction` | Carried verbatim from the bundle (empty = unmeasured) |
| `ed_class` | `activator` / `repressor` / `ambiguous` / `unmeasured` |

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
| `motif_score` | Normalized PWM log-odds score (−1 to 1), best of both strands; `NA` if no JASPAR motif |
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

**Expression host (`--species`).** Gate-1 self-tolerance needs the proteome of the organism the fusion will be expressed in. This follows the same `--species "Homo sapiens"` style as `cli.py`, but is **optional**: if omitted, the organism is **auto-detected** from the `species` column of `--dbd-input` (which `cli.py` populated from your Step-1 `--species`). Pass `--species` explicitly only to override — e.g. a human DBD expressed in a murine model. See [Self-tolerance filtering](#key-design-points) for where the FASTA files go.

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
| Function + Aggregation | AF3, FoldX, built-in AGGRESCAN-style scorer | DBD–DNA binding retained; junction aggregation | **HPCC/GPU** |
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
| **AGGRESCAN scorer** (built-in) | Gate 2 (aggregation) | open | none — uses `freesasa` (already in `environment.yml`) |
| **NetMHCpan 4.2** | Gate 1 (MHC-I) | free academic | [DTU download](https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/) → tarball |
| **NetCTLpan 1.1** | Gate 1 (deprecated) | free academic | discontinued by DTU — superseded by NetMHCpan; legacy installs only |
| **FoldX 5** | structure phase | free academic | [foldxsuite.crg.eu](https://foldxsuite.crg.eu/) → tarball |
| **CamSol** | Gate 2 (alt) | web server | no CLI — use the built-in scorer, or the `api` backend |
| **UbPred** | Gate 3 (optional) | web server | no CLI — Gate 3 runs without it |

#### Local install (helper script)

`install_fusion_tools.sh` automates the licence-gated tools. (Gate 2 aggregation needs no install — it's the built-in AGGRESCAN-style scorer, which only uses `freesasa` from `environment.yml`.) For the licence-gated tools (**NetMHCpan 4.2** and **FoldX**) you first **register and download the tarballs** from the links above and drop them into a `vendor/` directory — the script then extracts them, patches the NetMHCpan tcsh wrapper (`NMHOME` + `TMPDIR`), fetches NetMHCpan's data files, and symlinks the binaries into one bin dir. It's idempotent and skips any tool whose tarball isn't present, so you can re-run it as you obtain each licence. (NetCTLpan is discontinued by DTU; if you still have a legacy `netCTLpan-*` tarball in `vendor/`, the script will configure it too, but new setups should use NetMHCpan.)

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
    aggrescan3d: { backend: local }                                  # built-in scorer (freesasa)
    # leave camsol / ubpred as 'disabled'
```

> **tcsh note:** the NetMHCpan/NetCTLpan wrapper scripts start with `#! /bin/tcsh`. If your `tcsh` is from conda (e.g. `~/miniconda3/envs/ecr/bin/tcsh`) rather than `/bin/tcsh`, the wrapper fails with `bad interpreter: /bin/tcsh`. The installer repoints the shebang automatically; if you unpacked manually, either `sudo apt-get install tcsh` (puts it at `/bin/tcsh`) or edit line 1 of the `netMHCpan` script to your tcsh path. Verify the install with the bundled self-test: `cd test && ../netMHCpan -p test.pep` should reproduce `test.pep.out`.

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
- **Self-tolerance filtering.** Because both domains are endogenous, only junction-spanning peptides are potential neoepitopes — and any that occur verbatim in the host proteome are still self, so they're subtracted before the MHC-I scan. The proteome is resolved per run in this priority order: an explicit `fusion.self_proteome` path > `--species` > the organism auto-detected from `--dbd-input`. The `--species` form maps to `<proteome_dir>/<slug>.fasta` (default `data/proteomes/`, e.g. `"Homo sapiens"` → `data/proteomes/homo_sapiens.fasta`), so you pick the organism per run without editing config. **Drop the FASTA files there yourself** (the dir is gitignored — proteomes are large); download the reviewed canonical set from UniProt:
  ```bash
  mkdir -p data/proteomes
  # Human (UP000005640)
  curl -o data/proteomes/homo_sapiens.fasta.gz \
    'https://rest.uniprot.org/uniprotkb/stream?query=(proteome:UP000005640)+AND+(reviewed:true)&format=fasta&compressed=true'
  # Mouse (UP000000589)
  curl -o data/proteomes/mus_musculus.fasta.gz \
    'https://rest.uniprot.org/uniprotkb/stream?query=(proteome:UP000000589)+AND+(reviewed:true)&format=fasta&compressed=true'
  gunzip data/proteomes/*.fasta.gz
  ```
  The reviewed canonical set (~13 MB human, ~13–15 MB mouse) is the right tier — don't pull the full TrEMBL proteome (~60 MB, mostly redundant predicted entries) for exact-match filtering. `check_fusion_env.sh` reports which `<species>.fasta` files it finds.
- **%rank, not raw IC₅₀.** Binders are flagged by `%rank ≤ rank_threshold` (default 2.0) across an HLA-I panel; the gate reports **epitope density** (flagged / tested), not single hits.
- **Pareto over composite.** The output marks the Pareto-optimal set across the liability axes rather than collapsing them into one opaque score (a `risk_score` is also provided for quick sorting).
- **Degradation ↔ presentation tension.** Proteasomal degradation (Gate 3) is what *generates* MHC-I peptides (Gate 1); the pipeline surfaces both and does not auto-resolve the trade-off.

### Output columns

`candidate`, `dbd`, `linker`, `ed`, `length`, `immuno_density`, `immuno_flagged`, `immuno_min_rank`, `aggregation_risk`, `stability_risk`, `fusion_ddg_kcal_mol`, `function_delta_ddg`, `af3_cif_path`, `risk_score`, `pareto_optimal`.

---

## Step 4 — Preview (post-eCR accessibility)

`preview.py` is a read-only reporting stage: after the design chain it renders a
**self-contained, interactive HTML** of the predicted post-eCR chromatin accessibility —
`measured start -> predicted (model) -> measured end` — with a live efficacy slider. The
predicted column is a transparent direct-effect interpolation
(`start + efficacy*(end - start)`, targeted regions only), flanked by the two measured
endpoints; no HPCC/external tools.

```bash
python preview.py \
  --bundle ../eCR_navigator/examples/run_bundle/<run> \
  --atpm   atpm_union.tsv \
  --output accessibility_preview.html
```

Accessibility resolves `--atpm` > `<bundle>/accessibility.tsv` > error. Optional `--genes` /
`--highlight` / `--annot` add gene labels, per-locus construct cards, and a flag column. Full
input/output contract: [docs/preview_contract.md](docs/preview_contract.md).

## Off-target evaluation (`offtarget.py`)

A DBD-level **specificity** screen. For each DBD's JASPAR motif, scan the *accessible* genomic regions of the reprogramming trajectory for unintended strong binding sites and summarize them as an **off-target effect score**. A TF can only bind open chromatin, so the background is the **union of accessibility peaks across states** (e.g. MEF ∪ mES ATAC peaks), not the whole genome.

```bash
python offtarget.py \
  --peaks MEF=MEF.e7_peaks.bed mES=mES.e7_peaks.bed \
  --from MEF --to mES \
  --genome /path/to/mm10.fa \
  --dbd-input predictions.tsv \
  --output offtarget.tsv
```

Build/inspect the background region set alone (no genome needed):
```bash
python offtarget.py --peaks MEF=MEF.e7_peaks.bed mES=mES.e7_peaks.bed \
  --from MEF --to mES --stop-after-regions --regions-out union_regions.bed
```

**Two-tier severity weighting** — the score is Σ over off-target sites of (site strength × region weight):
- **Tier 1 (default):** region weight from accessibility *dynamics* — regions that **open** (target-state only) or **close** (start-state only) during the transition outweigh **shared/constitutive** background (default 1.0 / 1.0 / 0.3). A provisional proxy for driver importance, derived from the peak sets alone.
- **Tier 2 (`--weights weights.tsv`):** swap in per-region driver weights from an **eCR_navigator** run bundle. Pass the bundle's `weights.tsv`; each off-target region takes the max `driver_score` of the navigator rows it overlaps, falling back to its Tier-1 dynamics weight where the bundle has no coverage. The scoring code is unchanged — only the weight source:
  ```bash
  python offtarget.py \
    --peaks MEF=MEF.e7_peaks.bed mES=mES.e7_peaks.bed --from MEF --to mES \
    --genome /path/to/mm10.fa --dbd-input predictions.tsv \
    --weights ../eCR_navigator/examples/run_bundle/<run>/weights.tsv \
    --output offtarget.tsv
  ```
  The `weights.tsv` assembly must also match the peak coordinates. This is the same `offtarget.tsv` that `target.py --offtarget` consumes, so running it with `--weights` makes the specificity axis navigator-driver-weighted.

| Output column | Meaning |
|---|---|
| `n_offtarget_sites` | # accessible regions with a site above `--threshold` (default 0.8) |
| `offtarget_effect_score` | Σ (site strength × region weight) — severity-weighted |
| `offtarget_unweighted` | Σ site strength (weights ignored) — for comparison |

The genome FASTA assembly **must match the peak coordinates** (e.g. mm10). DBDs without a JASPAR motif return `NA`. Reuses the Step-1 PWM core (`scan.fetch_pwm` / `scan.best_site_score`), which scans both strands.

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
│   ├── scan.py           # JASPAR PWM scanning, both strands (parallel fetch)
│   ├── target.py         # region → binding site (consumes eCR_navigator bundle)
│   ├── offtarget.py      # DBD off-target / specificity screen (Tier-1/Tier-2)
│   ├── score.py          # validation_level → annotation_confidence
│   ├── output.py         # table formatting and TSV output
│   ├── prefetch.py       # pre-download JASPAR motifs to local cache
│   ├── filter.py         # confidence + score filtering
│   ├── fimo.py           # FIMO validation (JASPAR → MEME format + fimo call)
│   ├── af3.py            # AF3 prediction (local / hpcc / online backends)
│   ├── foldx.py          # FoldX affinity estimation
│   ├── fusion/           # Step 3 fusion-design developability gates
│   └── config.py         # load/validate config.yaml
├── jaspar_cache/          # populated by server_setup.sh (gitignored)
├── af3_outputs/           # AF3 inputs, CIF outputs, run_log.json (gitignored)
├── target.py              # Step 0 entrypoint (region → site)
├── cli.py                 # Step 1 entrypoint
├── refine.py              # Step 2 entrypoint
├── fuse.py                # Step 3 entrypoint (fusion-design developability)
├── offtarget.py           # off-target evaluation entrypoint
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
- [pyfaidx](https://github.com/mdshw5/pyfaidx) — genome FASTA extraction (`target.py`, `offtarget.py`)
- [pandas](https://pandas.pydata.org/)
- [requests](https://requests.readthedocs.io/)
- [paramiko](https://www.paramiko.org/) — SSH/SFTP for HPCC backend
- [pyyaml](https://pyyaml.org/) — config file parsing
- [MEME Suite](https://meme-suite.org/) — FIMO stage (external tool, included in conda env)
- [FoldX](https://foldxsuite.crg.eu/) — affinity stage (external tool, free academic licence)
