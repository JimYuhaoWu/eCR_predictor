#!/usr/bin/env bash
# server_setup.sh — one-time setup on the server
# Run from the directory that contains both eCR_mod_lib/ and eCR_predictor/.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"   # parent dir containing both repos

MOD_LIB="$ROOT/eCR_mod_lib"
PREDICTOR="$SCRIPT_DIR"

# ── 1. Install Python deps ────────────────────────────────────────────────────
pip install -e "$MOD_LIB"
pip install -e "$PREDICTOR"

# ── 2. Populate the database ──────────────────────────────────────────────────
# Run the seeding scripts in order; each fetches data and writes to
# eCR_mod_lib/library/module_library.db.
cd "$MOD_LIB"
python scripts/01_fetch_dbd.py
python scripts/02_seed_ed.py
python scripts/03_fetch_cr.py
python scripts/04_build_library.py
python scripts/05_validate.py

echo ""
echo "Setup complete. Run predictions with:"
echo "  bash $PREDICTOR/server_run.sh <sequence> \"<species>\" [output.tsv]"
