"""
ECR_predictor refinement pipeline.

Takes the output.tsv from cli.py and runs:
  [1] Filter  — drop low-confidence hits below the motif_score threshold
  [2] FIMO    — validate remaining hits with FIMO (requires MEME Suite)
  [3] AF3     — predict DBD–DNA complex structures for top hits
  [4] FoldX   — estimate binding affinity from AF3 structures (stub)

Usage:
    python refine.py --input results.tsv --sequence ATCG... [--config config.yaml]

AF3 backend is selected in config.yaml (local | hpcc | online).
Run with --stop-after fimo to stop before AF3.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


STAGES = ["filter", "fimo", "af3", "foldx"]


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Refine ECR_predictor hits: filter → FIMO → AF3 → FoldX."
    )
    parser.add_argument(
        "--input", required=True,
        help="Path to the predictor output TSV (from cli.py --output).",
    )
    parser.add_argument(
        "--sequence", required=True,
        help="Original DNA query sequence (same as used in cli.py).",
    )
    parser.add_argument(
        "--output", default=None,
        help="Output TSV path (default: <input stem>_refined.tsv).",
    )
    parser.add_argument(
        "--config", default=None,
        help="Path to config.yaml (default: config.yaml in repo root).",
    )
    parser.add_argument(
        "--min-motif-score", type=float, default=0.0,
        help="Drop hits with motif_score below this (default: 0.0).",
    )
    parser.add_argument(
        "--fimo-pvalue", type=float, default=1e-4,
        help="FIMO p-value threshold for a hit to be validated (default: 1e-4).",
    )
    parser.add_argument(
        "--top-n-af3", type=int, default=2,
        help="Number of top FIMO-validated hits to submit to AF3 (default: 2).",
    )
    parser.add_argument(
        "--stop-after",
        choices=STAGES,
        default=None,
        help="Stop the pipeline after this stage (inclusive).",
    )
    parser.add_argument(
        "--af3-output-dir", default="af3_outputs",
        help="Directory for AF3 JSON inputs and structure outputs (default: af3_outputs/).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"ERROR: input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output) if args.output else input_path.with_name(
        input_path.stem + "_refined" + input_path.suffix
    )

    sequence = args.sequence.upper().strip()
    invalid = set(sequence) - set("ACGTN")
    if invalid:
        print(f"ERROR: sequence contains non-DNA characters: {invalid}", file=sys.stderr)
        sys.exit(1)

    stop_after = args.stop_after

    # -------------------------------------------------------------------------
    # Load config
    # -------------------------------------------------------------------------
    from ecr_predictor.config import load_config, get_af3_config
    try:
        cfg = load_config(args.config)
    except FileNotFoundError as e:
        print(f"WARNING: {e}\nAF3 stage will use default hpcc settings.", file=sys.stderr)
        cfg = {}
    af3_cfg = get_af3_config(cfg)

    # -------------------------------------------------------------------------
    # Load
    # -------------------------------------------------------------------------
    print(f"Loading {input_path}...", file=sys.stderr)
    df = pd.read_csv(input_path, sep="\t", dtype=str)
    print(f"  {len(df)} hits loaded.", file=sys.stderr)

    # -------------------------------------------------------------------------
    # [1] Filter
    # -------------------------------------------------------------------------
    print("\n[1/4] Filtering low-confidence hits...", file=sys.stderr)
    from ecr_predictor.filter import filter_hits
    df = filter_hits(df, min_motif_score=args.min_motif_score)
    print(f"      {len(df)} hits remain after filtering.", file=sys.stderr)

    if df.empty:
        print("WARNING: all hits were filtered out. Exiting.", file=sys.stderr)
        _write(df, output_path)
        return

    if stop_after == "filter":
        _write(df, output_path)
        return

    # -------------------------------------------------------------------------
    # [2] FIMO
    # -------------------------------------------------------------------------
    print("\n[2/4] Running FIMO validation...", file=sys.stderr)
    from ecr_predictor.fimo import run_fimo_validation, FIMONotAvailableError
    try:
        df = run_fimo_validation(df, sequence, pvalue_thresh=args.fimo_pvalue)
    except FIMONotAvailableError as e:
        print(f"WARNING: FIMO not available — skipping validation.\n  {e}", file=sys.stderr)
        df["fimo_pvalue"] = float("nan")
        df["fimo_validated"] = False

    if stop_after == "fimo":
        _write(df, output_path)
        return

    # -------------------------------------------------------------------------
    # [3] AF3
    # -------------------------------------------------------------------------
    print("\n[3/4] Running AlphaFold 3 structure prediction...", file=sys.stderr)
    from ecr_predictor.af3 import run_af3_prediction
    try:
        df = run_af3_prediction(
            df,
            dna_sequence=sequence,
            af3_cfg=af3_cfg,
            output_dir=Path(args.af3_output_dir),
            top_n=args.top_n_af3,
        )
    except NotImplementedError as e:
        print(f"NOTE: AF3 backend not implemented — stopping after FIMO output.\n  {e}", file=sys.stderr)
        _write(df, output_path)
        return

    if stop_after == "af3":
        _write(df, output_path)
        return

    # -------------------------------------------------------------------------
    # [4] FoldX
    # -------------------------------------------------------------------------
    print("\n[4/4] Estimating binding affinity with FoldX...", file=sys.stderr)
    from ecr_predictor.foldx import run_foldx_affinity
    try:
        df = run_foldx_affinity(df)
    except (NotImplementedError, FileNotFoundError) as e:
        print(f"NOTE: FoldX not yet implemented — stopping after AF3 output.\n  {e}", file=sys.stderr)
        _write(df, output_path)
        return

    _write(df, output_path)


def _write(df: pd.DataFrame, path: Path) -> None:
    path.write_text(df.to_csv(sep="\t", index=False, na_rep="NA"), encoding="utf-8")
    print(f"\nOutput written to: {path}", file=sys.stderr)


if __name__ == "__main__":
    main()
