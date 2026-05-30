#!/usr/bin/env bash
# server_setup.sh — one-time setup on the server
# Run once after cloning or when setting up a fresh environment.
set -euo pipefail

# ── 1. Clone repos ────────────────────────────────────────────────────────────
git clone https://github.com/JimYuhaoWu/eCR_mod_lib.git
git clone https://github.com/JimYuhaoWu/eCR_predictor.git

# ── 2. Install Python deps ────────────────────────────────────────────────────
pip install -e eCR_mod_lib
pip install -e eCR_predictor

# ── 3. Populate the database ──────────────────────────────────────────────────
# Run the seeding scripts in order from inside the library repo.
# Each script fetches data and writes to library/module_library.db.
cd eCR_mod_lib
python scripts/01_fetch_dbd.py
python scripts/02_seed_ed.py
python scripts/03_fetch_cr.py
python scripts/04_build_library.py
python scripts/05_validate.py
cd ..

echo ""
echo "Setup complete. Run predictions with:"
echo "  python eCR_predictor/cli.py --sequence <SEQ> --species \"Homo sapiens\" --db eCR_mod_lib/library/module_library.db"
