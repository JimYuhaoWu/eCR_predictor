#!/usr/bin/env bash
#
# check_fusion_env.sh — read-only readiness check for the fusion stage (fuse.py, Step 3).
#
# Reports what is working and what is missing WITHOUT running a full pipeline:
#   - Python env + required packages
#   - eCR_predictor / eCR_mod_lib install + library DB (ED/DBD/CR counts)
#   - config.yaml, and for each gate the backend it will ACTUALLY use
#     (local binary present? api url/key set? disabled?)
#   - self-proteome, freesasa, structure-phase prereqs (AF3 / FoldX / MEME)
#   - a tiny SYNTHETIC dry-run through assemble + sequence gates (no HPCC/GPU)
#
# Universal: no hard-coded usernames/paths. Run from anywhere; it locates the repo
# from its own location. Diagnostic only — always exits 0 so it never blocks a
# provisioning script. Meaning is in the PASS/WARN/FAIL lines and the summary.
#
# Usage:
#   bash check_fusion_env.sh [--config path/to/config.yaml]
#
set -u

CONFIG_ARG=""
while [ $# -gt 0 ]; do
  case "$1" in
    --config) CONFIG_ARG="${2:-}"; shift 2 ;;
    -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

# Repo root = directory of this script.
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
export ECR_CHECK_REPO_DIR="$REPO_DIR"
export ECR_CHECK_CONFIG="$CONFIG_ARG"

PY="$(command -v python3 || command -v python || true)"
if [ -z "$PY" ]; then
  echo "FAIL  No python interpreter found on PATH. Activate the 'ecr' conda env first." >&2
  exit 0
fi

echo "eCR_predictor — fusion-stage readiness check"
echo "repo: $REPO_DIR"
echo "python: $PY"
echo

cd "$REPO_DIR" || exit 0
"$PY" - <<'PYEOF'
import os, sys, shutil, importlib, tempfile, subprocess
from pathlib import Path

REPO = Path(os.environ["ECR_CHECK_REPO_DIR"]).resolve()
CFG_ARG = os.environ.get("ECR_CHECK_CONFIG", "").strip()

# ----- pretty printing -----------------------------------------------------
USE_COLOR = sys.stdout.isatty() and os.environ.get("NO_COLOR") is None
def _c(code, s): return f"\033[{code}m{s}\033[0m" if USE_COLOR else s
GREEN, YELLOW, RED, BOLD = "32", "33", "31", "1"

counts = {"PASS": 0, "WARN": 0, "FAIL": 0}
def line(status, msg, hint=""):
    counts[status] = counts.get(status, 0) + 1
    tag = {"PASS": _c(GREEN, "PASS"), "WARN": _c(YELLOW, "WARN"), "FAIL": _c(RED, "FAIL")}[status]
    print(f"  {tag}  {msg}")
    if hint:
        print(f"        -> {hint}")

def section(title):
    print(_c(BOLD, title))

# ----- 1. Python + packages ------------------------------------------------
section("[1] Python environment")
v = sys.version_info
if v[:2] >= (3, 8):
    line("PASS", f"python {v.major}.{v.minor}.{v.micro}")
else:
    line("FAIL", f"python {v.major}.{v.minor}", "need >= 3.8 (server constraint)")

# (module, pip/conda name, required?)  required -> FAIL if missing, else WARN
PKGS = [
    ("pandas", "pandas", True),
    ("yaml", "pyyaml", True),
    ("Bio", "biopython", True),       # junction structure parsing / FoldX CIF
    ("requests", "requests", False),  # only for api backends
    ("paramiko", "paramiko", False),  # only for hpcc structure phase
    ("pyotp", "pyotp", False),        # only for hpcc TOTP auth
    ("freesasa", "freesasa", False),  # Gate 2 structure-weighted aggregation
    ("tqdm", "tqdm", False),
]
for mod, name, required in PKGS:
    try:
        importlib.import_module(mod)
        line("PASS", f"import {mod}")
    except Exception as e:
        if required:
            line("FAIL", f"import {mod} ({name})", f"pip/conda install {name} — {e}")
        else:
            line("WARN", f"import {mod} ({name}) not available",
                 f"needed only for some backends; conda install -c conda-forge {name}")

# ----- 2. Package installs (both repos) ------------------------------------
section("[2] Package installs")
try:
    import ecr_predictor  # noqa
    line("PASS", "ecr_predictor importable")
except Exception as e:
    line("FAIL", "ecr_predictor not importable", f"pip install -e . in repo root — {e}")
try:
    from ecr_predictor.query import load_eds  # noqa
    line("PASS", "eCR_mod_lib wired (query.load_eds present)")
    _have_query = True
except Exception as e:
    # Most common cause: bare `pip` points at a different interpreter than the
    # active env, so eCR_mod_lib ('scripts' package) got installed elsewhere.
    hint = "install it INTO THIS interpreter: python -m pip install -e ../eCR_mod_lib"
    if "No module named 'scripts'" in str(e):
        hint += "  (do NOT use bare `pip` — verify with `python -m pip -V` vs `pip -V`)"
    line("FAIL", f"eCR_mod_lib / query.load_eds unavailable — {e}", hint)
    _have_query = False

# ----- 3. Library DB (EDs come from here) ----------------------------------
section("[3] Module library DB")
ed_names = []
try:
    # DB path is derivable without importing query (which needs the eCR_mod_lib
    # 'scripts' package) so a missing install doesn't hide whether the DB exists.
    db_path = REPO.parent / "eCR_mod_lib" / "library" / "module_library.db"
    if not db_path.exists():
        line("FAIL", f"DB not found: {db_path}",
             "sync module_library.db, or rebuild it in eCR_mod_lib (server_setup.sh)")
    else:
        line("PASS", f"DB present: {db_path}")
        try:
            from ecr_predictor.query import load_dbds
        except Exception:
            load_dbds = None
        # ED count (the seed the fusion stage needs)
        if _have_query:
            try:
                eds = load_eds()
                eds = eds.dropna(subset=["sequence_aa"])
                n_ed = len(eds)
                if n_ed > 0:
                    line("PASS", f"EDs with sequence: {n_ed}")
                    name_col = "name" if "name" in eds.columns else (
                        "gene_symbol" if "gene_symbol" in eds.columns else None)
                    if name_col:
                        ed_names = [str(x) for x in eds[name_col].head(2).tolist()]
                else:
                    line("FAIL", "0 EDs with sequence in DB",
                         "rebuild the ED table in eCR_mod_lib")
            except Exception as e:
                line("FAIL", "load_eds() failed", str(e))
        if load_dbds is not None:
            try:
                n_dbd = len(load_dbds().dropna(subset=["sequence_aa"]))
                line("PASS" if n_dbd else "WARN", f"DBDs with sequence: {n_dbd}",
                     "" if n_dbd else "DBDs come from your Step-1/2 TSV, not the DB — informational")
            except Exception:
                pass
except Exception as e:
    line("FAIL", "library query layer error", str(e))

# ----- 4. Config + per-gate backend resolution -----------------------------
section("[4] Config & gate backends")
cfg = None
fusion_cfg = {}
cfg_path = None
try:
    from ecr_predictor.config import load_config, get_fusion_config
    candidates = [CFG_ARG] if CFG_ARG else []
    candidates += [str(REPO / "config.yaml")]
    picked = next((c for c in candidates if c and Path(c).exists()), None)
    if picked is None:
        example = REPO / "config.example.yaml"
        if example.exists():
            line("WARN", "no config.yaml — falling back to config.example.yaml",
                 "cp config.example.yaml config.yaml, then edit for your site")
            picked = str(example)
        else:
            line("FAIL", "no config.yaml and no config.example.yaml found", "")
    else:
        line("PASS", f"config: {picked}")
    if picked:
        cfg_path = picked
        cfg = load_config(picked)
        fusion_cfg = get_fusion_config(cfg)
except Exception as e:
    line("FAIL", "config load failed", str(e))

def _resolve(tool_cfg):
    try:
        from ecr_predictor.fusion import backends
        return backends.resolve_backend(tool_cfg or {})
    except Exception:
        return str((tool_cfg or {}).get("backend", "disabled")).lower()

def check_tool(gate, keys, builtin_local=False):
    """keys: list of tool names to consider in priority order for this gate."""
    tools = (fusion_cfg or {}).get("tools", {})
    chosen = None
    for name in keys:
        tcfg = tools.get(name, {})
        if _resolve(tcfg) != "disabled":
            chosen = (name, tcfg); break
    if chosen is None:
        line("WARN", f"{gate}: no tool enabled ({'/'.join(keys)})",
             "gate will be SKIPPED — set a backend to local|api")
        return
    name, tcfg = chosen
    backend = _resolve(tcfg)
    if backend == "local":
        if builtin_local and name in ("aggrescan3d",):
            line("PASS", f"{gate}: {name} (local, built-in scorer — no binary)")
        else:
            cmd = (tcfg.get("local", {}) or {}).get("command", name)
            if shutil.which(cmd):
                line("PASS", f"{gate}: {name} (local: {cmd} on PATH)")
            else:
                line("FAIL", f"{gate}: {name} local binary '{cmd}' NOT on PATH",
                     "install it (see install_fusion_tools.sh) or set backend: api/disabled — "
                     "gate is SKIPPED as-is")
    elif backend == "api":
        api = tcfg.get("api", {}) or {}
        url = api.get("url", "")
        env_name = api.get("api_key_env", "")
        key = os.environ.get(env_name, "") if env_name else api.get("api_key", "")
        if not url:
            line("FAIL", f"{gate}: {name} (api) has empty url", "set fusion.tools." + name + ".api.url")
        elif env_name and not key:
            line("WARN", f"{gate}: {name} (api) url set but ${env_name} unset",
                 f"export {env_name}=... if the service needs auth")
        else:
            line("PASS", f"{gate}: {name} (api: {url})")

if fusion_cfg:
    check_tool("Gate 1 immunogenicity", ["netmhcpan", "netctlpan"])
    check_tool("Gate 2 aggregation", ["aggrescan3d", "camsol"], builtin_local=True)
    check_tool("Gate 3 stability (UbPred, optional)", ["ubpred"])
    # linkers
    linkers = fusion_cfg.get("linkers", [])
    line("PASS" if linkers else "FAIL", f"linkers configured: {len(linkers)}",
         "" if linkers else "add fusion.linkers in config.yaml")
    # self proteome — explicit config path wins, else species resolves to
    # <proteome_dir>/<slug>.fasta at run time; report which organism files exist.
    sp = fusion_cfg.get("self_proteome", "")
    if sp:
        if Path(sp).exists():
            line("PASS", f"self_proteome (explicit config): {sp}")
        else:
            line("FAIL", f"self_proteome path missing: {sp}", "fix the path or clear it to use --species resolution")
    else:
        pdir = Path(fusion_cfg.get("proteome_dir", "data/proteomes"))
        if not pdir.is_absolute():
            pdir = REPO / pdir
        fastas = sorted(pdir.glob("*.fasta")) if pdir.is_dir() else []
        if fastas:
            names = ", ".join(f.stem for f in fastas)
            line("PASS", f"proteome files in {pdir}: {names}",
                 "Gate 1 self-tolerance uses <slug>.fasta matching the run's species "
                 "(--species or auto-detected from --dbd-input)")
        else:
            line("WARN", f"no <species>.fasta in {pdir}",
                 "Gate 1 self-tolerance OFF unless present — put e.g. homo_sapiens.fasta / "
                 "mus_musculus.fasta there (or set fusion.self_proteome)")

# ----- 5. Structure-phase prereqs (optional; only past --stop-after prune) --
section("[5] Structure phase (optional)")
af3_backend = "?"
try:
    from ecr_predictor.config import get_af3_config
    if cfg is not None:
        af3_backend = str(get_af3_config(cfg).get("backend", "?"))
except Exception:
    pass
line("PASS", f"af3.backend = {af3_backend}",
     "hpcc needs paramiko + creds; online needs ECR_CHAI_API_KEY; local needs AF3 on this host")
if os.environ.get("FOLDX_PATH") and Path(os.environ["FOLDX_PATH"]).exists():
    line("PASS", f"FOLDX_PATH = {os.environ['FOLDX_PATH']}")
else:
    line("WARN", "FOLDX_PATH unset/invalid",
         "export FOLDX_PATH=/path/to/foldx — needed for function-retention ΔΔG")
if shutil.which("fimo"):
    line("PASS", "MEME 'fimo' on PATH (Steps 1/2)")
else:
    line("WARN", "MEME 'fimo' not on PATH", "conda install -c bioconda meme — needed for Step 2 FIMO")

# ----- 6. Synthetic dry-run (assemble + sequence gates, no HPCC/GPU) --------
section("[6] Synthetic dry-run (assemble + sequence gates)")
if cfg_path and ed_names:
    tmp = Path(tempfile.mkdtemp(prefix="ecr_fuse_check_"))
    dbd_tsv = tmp / "dbd.tsv"
    dbd_tsv.write_text(
        "gene_name\tsequence_aa\ttf_family\n"
        "CHECK_DBD\tMKQLEDKVEELLSKNYHLENEVARLKKLVGER\tbZIP\n",
        encoding="utf-8",
    )
    out_tsv = tmp / "out.tsv"
    cmd = [sys.executable, str(REPO / "fuse.py"),
           "--dbd-input", str(dbd_tsv),
           "--ed", ",".join(ed_names),
           "--config", cfg_path,
           "--stop-after", "sequence",
           "--output", str(out_tsv)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if out_tsv.exists() and proc.returncode == 0:
            nrows = max(0, len(out_tsv.read_text().splitlines()) - 1)
            line("PASS", f"dry-run OK — assembled + scored {nrows} candidates "
                         f"(1 synthetic DBD x linkers x {len(ed_names)} EDs)")
            tail = [l for l in proc.stderr.splitlines() if l.strip().startswith(("SKIP", "[Gate", "      SKIP"))]
            for l in tail:
                print(f"        {l.strip()}")
        else:
            line("FAIL", f"dry-run failed (exit {proc.returncode})",
                 (proc.stderr.strip().splitlines() or ["no stderr"])[-1])
    except subprocess.TimeoutExpired:
        line("FAIL", "dry-run timed out (>300s)",
             "unexpected for the sequence stage — check for a hung tool call")
    except Exception as e:
        line("FAIL", "dry-run could not start", str(e))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
elif not ed_names:
    line("WARN", "dry-run skipped — no ED names available (see [3])")
else:
    line("WARN", "dry-run skipped — no config available (see [4])")

# ----- summary -------------------------------------------------------------
print()
print(_c(BOLD, "Summary: ")
      + f"{_c(GREEN, str(counts['PASS']) + ' pass')}, "
      + f"{_c(YELLOW, str(counts['WARN']) + ' warn')}, "
      + f"{_c(RED, str(counts['FAIL']) + ' fail')}")
if counts["FAIL"]:
    print("  FAIL items block a run (or a specific gate). WARN items degrade quality/optional phases.")
else:
    print("  No blockers. WARN items are optional or affect only the structure phase / gate quality.")
PYEOF

exit 0
