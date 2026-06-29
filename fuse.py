"""
ECR_predictor fusion-design pipeline (Step 3).

Assembles fusion candidates (DBD + linker + ED) and screens them for
developability before wet-lab synthesis. Target modality: intracellular
expression (viral vector / mRNA) → MHC-I junction neoepitopes are the primary
immune liability.

Ordering principle: run cheap SEQUENCE-based gates on the whole library first,
prune to the Pareto-optimal survivors, and only then spend HPCC/GPU on the
STRUCTURE phase (AF3 fold + FoldX) for the handful that remain.

  [1] assemble   DBDs × linkers × EDs → candidate library
  [2] sequence   Gate 1 (MHC-I immunogenicity) + Gate 3 (stability)  — all candidates
  [3] prune      Pareto over sequence axes → survivors
  [4] structure  AF3 fold survivors + FoldX ΔΔG (function) + Gate 2 aggregation
  [5] score      final composite risk + Pareto ranking

DBDs are read from a TSV with gene_name + sequence_aa (Step-1/2 output).
EDs are read from the eCR_mod_lib library (type='ED'). Linkers and tool backends
come from config.yaml (fusion section).

Usage:
    python fuse.py --dbd-input predictions.tsv --sequence ACGT... --config config.yaml
    python fuse.py --dbd-input predictions.tsv --stop-after prune
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

STAGES = ["assemble", "sequence", "prune", "structure", "score"]


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Assemble and screen DBD–ED fusion candidates."
    )
    p.add_argument("--dbd-input", required=True,
                   help="TSV with gene_name + sequence_aa columns (Step 1/2 output).")
    p.add_argument("--dbd", default=None,
                   help="Comma-separated gene_name subset to use as DBDs (default: all).")
    p.add_argument("--ed", default=None,
                   help="Comma-separated ED name subset from the library (default: all).")
    p.add_argument("--sequence", default=None,
                   help="DNA target the DBD must still bind (required for the structure phase).")
    p.add_argument("--config", default=None, help="Path to config.yaml.")
    p.add_argument("--output", default="fusion_candidates.tsv", help="Output TSV path.")
    p.add_argument("--top-n-structure", type=int, default=10,
                   help="Max survivors to fold in the structure phase (default: 10).")
    p.add_argument("--af3-output-dir", default="af3_outputs",
                   help="AF3 output root for the structure phase (default: af3_outputs/).")
    p.add_argument("--stop-after", choices=STAGES, default=None,
                   help="Stop the pipeline after this stage (inclusive).")
    return p.parse_args(argv)


def _load_dbds(input_path: Path, subset: str | None) -> tuple[dict[str, str], dict[str, dict]]:
    """Return ({gene: sequence}, {gene: meta}) where meta carries tf_family,
    zinc_finger_count, and a FoldX baseline ΔΔG if present in the TSV."""
    df = pd.read_csv(input_path, sep="\t", dtype=str)
    if "sequence_aa" not in df.columns:
        print("ERROR: --dbd-input needs a 'sequence_aa' column "
              "(re-run cli.py with --include-sequence).", file=sys.stderr)
        sys.exit(1)
    df = df.dropna(subset=["sequence_aa"])
    if subset:
        wanted = {s.strip() for s in subset.split(",")}
        df = df[df["gene_name"].isin(wanted)]

    dbds: dict[str, str] = {}
    meta: dict[str, dict] = {}
    for _, row in df.iterrows():
        gene = str(row["gene_name"])
        if gene in dbds:
            continue
        dbds[gene] = str(row["sequence_aa"])
        baseline = pd.to_numeric(row.get("foldx_ddg_kcal_mol"), errors="coerce")
        meta[gene] = {
            "tf_family": row.get("tf_family", ""),
            "zinc_finger_count": row.get("zinc_finger_count"),
            "baseline_ddg": None if pd.isna(baseline) else float(baseline),
        }
    return dbds, meta


def _load_eds(subset: str | None) -> dict[str, str]:
    """Load effector domains (type='ED') from the eCR_mod_lib library."""
    from ecr_predictor.query import load_eds
    df = load_eds()
    df = df.dropna(subset=["sequence_aa"])
    name_col = "name" if "name" in df.columns else "gene_symbol"
    if subset:
        wanted = {s.strip() for s in subset.split(",")}
        df = df[df[name_col].isin(wanted)]
    eds: dict[str, str] = {}
    for _, row in df.iterrows():
        eds.setdefault(str(row[name_col]), str(row["sequence_aa"]))
    return eds


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    input_path = Path(args.dbd_input)
    if not input_path.exists():
        print(f"ERROR: input not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    from ecr_predictor.config import load_config, get_fusion_config, get_af3_config
    try:
        cfg = load_config(args.config)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    fusion_cfg = get_fusion_config(cfg)

    # -------------------------------------------------------------------------
    # [1] Assemble
    # -------------------------------------------------------------------------
    print("[1/5] Assembling fusion candidates...", file=sys.stderr)
    from ecr_predictor.fusion.assemble import assemble_fusions
    dbds, dbd_meta = _load_dbds(input_path, args.dbd)
    eds = _load_eds(args.ed)
    linkers = fusion_cfg.get("linkers", [])
    if not dbds or not eds or not linkers:
        print(f"ERROR: need ≥1 DBD ({len(dbds)}), ED ({len(eds)}), and linker "
              f"({len(linkers)}). Check --dbd-input, the ED library, and "
              f"fusion.linkers in config.yaml.", file=sys.stderr)
        sys.exit(1)
    candidates = assemble_fusions(
        dbds, eds, linkers, orientation=fusion_cfg.get("orientation", "dbd_n")
    )
    print(f"      {len(candidates)} candidates "
          f"({len(dbds)} DBD × {len(linkers)} linker × {len(eds)} ED).",
          file=sys.stderr)

    if args.stop_after == "assemble":
        _write_assembly(candidates, Path(args.output))
        return

    # -------------------------------------------------------------------------
    # Self-proteome blob (Gate 1)
    # -------------------------------------------------------------------------
    proteome_blob = None
    self_fa = fusion_cfg.get("self_proteome", "")
    if self_fa:
        from ecr_predictor.fusion.junction import load_proteome_blob
        print(f"      Loading self proteome: {self_fa}", file=sys.stderr)
        proteome_blob = load_proteome_blob(self_fa)

    from ecr_predictor.fusion import backends
    immuno = stability = aggregation = structure = None

    # -------------------------------------------------------------------------
    # [2] Sequence gates (cheap; all candidates)
    # -------------------------------------------------------------------------
    print("\n[2/5] Sequence gates — immunogenicity + stability...", file=sys.stderr)
    from ecr_predictor.fusion.immunogenicity import screen_immunogenicity
    from ecr_predictor.fusion.stability import screen_stability
    try:
        immuno = screen_immunogenicity(candidates, fusion_cfg, proteome_blob)
    except (backends.ToolDisabled, backends.ToolNotAvailableError) as e:
        print(f"      SKIP Gate 1: {e}", file=sys.stderr)
    stability = screen_stability(candidates, fusion_cfg)

    if args.stop_after == "sequence":
        _write_scores(candidates, immuno, aggregation, stability, structure, Path(args.output))
        return

    # -------------------------------------------------------------------------
    # [3] Prune — Pareto over sequence axes → survivors
    # -------------------------------------------------------------------------
    print("\n[3/5] Pruning to Pareto-optimal survivors...", file=sys.stderr)
    from ecr_predictor.fusion.score import (
        build_score_table, add_pareto_and_rank, SEQUENCE_AXES,
    )
    seq_table = add_pareto_and_rank(
        build_score_table(candidates, immuno, None, stability), axes=SEQUENCE_AXES,
    )
    survivor_names = _select_survivors(seq_table, args.top_n_structure)
    survivors = [c for c in candidates if c.name in survivor_names]
    print(f"      {len(survivors)}/{len(candidates)} survivors advance to structure.",
          file=sys.stderr)

    if args.stop_after == "prune":
        seq_table.to_csv(Path(args.output), sep="\t", index=False, na_rep="NA")
        print(f"\nPruned table written to: {args.output}", file=sys.stderr)
        return

    # -------------------------------------------------------------------------
    # [4] Structure phase (HPCC/GPU; survivors only)
    # -------------------------------------------------------------------------
    print("\n[4/5] Structure phase — AF3 fold + FoldX + aggregation...", file=sys.stderr)
    if not args.sequence:
        print("ERROR: --sequence (DNA target) is required for the structure phase. "
              "Use --stop-after prune to skip it.", file=sys.stderr)
        sys.exit(1)
    dna = args.sequence.upper().strip()

    from ecr_predictor.fusion.structure import run_structure_phase
    from ecr_predictor.fusion.aggregation import screen_aggregation
    af3_cfg = get_af3_config(cfg)
    try:
        structure = run_structure_phase(
            survivors, dna, dbd_meta, af3_cfg, Path(args.af3_output_dir),
        )
    except (NotImplementedError, FileNotFoundError) as e:
        print(f"      NOTE: structure phase unavailable — {e}", file=sys.stderr)
        structure = None

    cif_paths = {n: d["af3_cif_path"] for n, d in (structure or {}).items()
                 if d.get("af3_cif_path")}
    try:
        aggregation = screen_aggregation(survivors, fusion_cfg, structure_paths=cif_paths)
    except (backends.ToolDisabled, backends.ToolNotAvailableError) as e:
        print(f"      SKIP Gate 2: {e}", file=sys.stderr)

    if args.stop_after == "structure":
        _write_scores(candidates, immuno, aggregation, stability, structure, Path(args.output))
        return

    # -------------------------------------------------------------------------
    # [5] Final score
    # -------------------------------------------------------------------------
    print("\n[5/5] Final composite risk + Pareto ranking...", file=sys.stderr)
    _write_scores(candidates, immuno, aggregation, stability, structure, Path(args.output))


def _select_survivors(seq_table: pd.DataFrame, top_n: int) -> set[str]:
    """Pareto-optimal candidates (sequence axes), capped at top_n by risk_score.
    Falls back to top_n by risk_score if Pareto flags are unavailable."""
    if "pareto_optimal" in seq_table.columns and seq_table["pareto_optimal"].any():
        winners = seq_table[seq_table["pareto_optimal"] == True]  # noqa: E712
    else:
        winners = seq_table
    winners = winners.sort_values("risk_score", na_position="last").head(top_n)
    return set(winners["candidate"].tolist())


def _write_assembly(candidates, path: Path) -> None:
    df = pd.DataFrame([{
        "candidate": c.name, "dbd": c.dbd_name, "linker": c.linker_name,
        "ed": c.ed_name, "length": len(c.sequence),
        "junctions": ";".join(f"{j.left}|{j.right}@{j.position}" for j in c.junctions),
        "sequence": c.sequence,
    } for c in candidates])
    path.write_text(df.to_csv(sep="\t", index=False), encoding="utf-8")
    print(f"\nAssembly written to: {path}", file=sys.stderr)


def _write_scores(candidates, immuno, aggregation, stability, structure, path: Path) -> None:
    from ecr_predictor.fusion.score import build_score_table, add_pareto_and_rank
    df = build_score_table(candidates, immuno, aggregation, stability, structure)
    df = add_pareto_and_rank(df)
    path.write_text(df.to_csv(sep="\t", index=False, na_rep="NA"), encoding="utf-8")
    print(f"\nScores written to: {path}", file=sys.stderr)


if __name__ == "__main__":
    main()
