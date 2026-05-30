"""
ECR_predictor — predict DBD binding for a DNA sequence.

Usage:
    python cli.py --sequence ATCG... --species "Homo sapiens" [--output results.tsv]
"""
from __future__ import annotations

import argparse
import sys

from ecr_predictor.output import build_result_table, write_output
from ecr_predictor.query import load_dbds, match_species
from ecr_predictor.scan import score_dbds
from ecr_predictor.score import assign_annotation_confidence


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Predict DBD binding candidates for a DNA sequence."
    )
    parser.add_argument("--sequence", required=True, help="DNA sequence to scan")
    parser.add_argument("--species", required=True, help='Query species, e.g. "Homo sapiens"')
    parser.add_argument("--output", default=None, help="Output TSV path (default: stdout)")
    parser.add_argument(
        "--db",
        default=None,
        help="Path to module_library.db (default: auto-detected)",
    )
    args = parser.parse_args(argv)

    # Validate sequence
    sequence = args.sequence.upper().strip()
    invalid = set(sequence) - set("ACGTN")
    if invalid:
        print(f"ERROR: sequence contains non-DNA characters: {invalid}", file=sys.stderr)
        sys.exit(1)

    # Load and filter DBDs
    from pathlib import Path
    db_path = Path(args.db) if args.db else None
    dbds = load_dbds(db_path) if db_path else load_dbds()
    candidates = match_species(dbds, args.species)

    if candidates.empty:
        print(
            f"WARNING: no DBDs found for species '{args.species}' (exact or genus match).",
            file=sys.stderr,
        )

    # Score
    motif_scores = score_dbds(candidates, sequence)
    confidence = assign_annotation_confidence(candidates)

    # Build and emit result table
    result = build_result_table(candidates, motif_scores, confidence)
    write_output(result, args.output)


if __name__ == "__main__":
    main()
