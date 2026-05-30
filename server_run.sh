#!/usr/bin/env bash
# server_run.sh — run a prediction
# Usage:
#   bash server_run.sh <sequence> <species> [output.tsv]
#
# Examples:
#   bash server_run.sh ACAGGAAGTGACAGGAAGTG "Homo sapiens"
#   bash server_run.sh ACAGGAAGTGACAGGAAGTG "Homo sapiens" results.tsv
set -euo pipefail

SEQUENCE="${1:?Usage: $0 <sequence> <species> [output.tsv]}"
SPECIES="${2:?Usage: $0 <sequence> <species> [output.tsv]}"
OUTPUT="${3:-}"

DB="eCR_mod_lib/library/module_library.db"

if [[ ! -f "$DB" ]]; then
    echo "ERROR: database not found at $DB — run server_setup.sh first." >&2
    exit 1
fi

ARGS=(--sequence "$SEQUENCE" --species "$SPECIES" --db "$DB")
if [[ -n "$OUTPUT" ]]; then
    ARGS+=(--output "$OUTPUT")
fi

python eCR_predictor/cli.py "${ARGS[@]}"
