#!/usr/bin/env bash
# server_setup.sh — one-time setup on the server
# Run from the directory that contains both eCR_mod_lib/ and eCR_predictor/.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"   # parent dir containing both repos

MOD_LIB="$ROOT/eCR_mod_lib"
PREDICTOR="$SCRIPT_DIR"
DB="$MOD_LIB/library/module_library.db"

# ── 1. Install Python deps ────────────────────────────────────────────────────
pip install -e "$MOD_LIB"
pip install -e "$PREDICTOR"

# ── 2. Populate the database (skip if already built) ─────────────────────────
if [[ -f "$DB" ]]; then
    echo "Database already exists at $DB — skipping seeding."
else
    echo "Building database..."
    cd "$MOD_LIB"
    python scripts/01_fetch_dbd.py
    python scripts/02_seed_ed.py
    python scripts/03_fetch_cr.py
    python scripts/04_build_library.py
    python scripts/05_validate.py
    cd "$SCRIPT_DIR"
fi

# ── 3. Pre-fetch JASPAR motifs into local cache ───────────────────────────────
echo "Pre-fetching JASPAR motifs..."
python -m ecr_predictor.prefetch --db "$DB" --cache-dir "$PREDICTOR/jaspar_cache"

echo ""
echo "Setup complete. Run predictions with:"
echo "  bash $PREDICTOR/server_run.sh <sequence> \"<species>\" [output.tsv]"
