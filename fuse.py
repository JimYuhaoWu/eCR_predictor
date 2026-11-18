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

DBDs are read from a TSV with gene_name + dbd_sequence_aa (the extracted DNA-binding
span; falls back to sequence_aa with a warning) -- Step-1/2 output.
EDs are read from the eCR_mod_lib library (type='ED'). Linkers and tool backends
come from config.yaml (fusion section).

Usage:
    python fuse.py --dbd-input predictions.tsv --sequence ACGT... --config config.yaml
    python fuse.py --dbd-input predictions.tsv --stop-after prune
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd

STAGES = ["assemble", "sequence", "prune", "structure", "score"]


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Assemble and screen DBD–ED fusion candidates."
    )
    p.add_argument("--dbd-input", required=True,
                   help="TSV with gene_name + dbd_sequence_aa (preferred) or "
                        "sequence_aa columns (Step 1/2 output).")
    p.add_argument("--dbd", default=None,
                   help="Comma-separated gene_name subset to use as DBDs (default: all).")
    p.add_argument("--ed", default=None,
                   help="Comma-separated ED name subset from the library (default: all).")
    p.add_argument("--sequence", default=None,
                   help="DNA target the DBD must still bind (required for the structure phase).")
    p.add_argument("--species", default=None,
                   help='Expression-host species for Gate-1 self-tolerance, e.g. '
                        '"Homo sapiens" (same style as cli.py). Default: auto-detected '
                        "from the --dbd-input 'species' column. Resolves to "
                        "data/proteomes/<slug>.fasta.")
    p.add_argument("--allow-fulllength-dbd", action="store_true",
                   help="Fuse on the full-length protein for records with no "
                        "extracted DBD span, instead of skipping them. Builds "
                        "constructs carrying every native domain of the TF -- "
                        "including KRAB on an opener.")
    p.add_argument("--allow-flagged-dbd", action="store_true",
                   help="Use spans with dbd_review_flag=1 (an excluded domain "
                        "lies inside the kept span) instead of skipping them. "
                        "Some flags are conservative: for REST, YY1, FOS and "
                        "HINFP the flagged region is itself required for "
                        "binding.")
    p.add_argument("--config", default=None, help="Path to config.yaml.")
    p.add_argument("--output", default="fusion_candidates.tsv", help="Output TSV path.")
    p.add_argument("--top-n-structure", type=int, default=10,
                   help="Max survivors to fold in the structure phase (default: 10).")
    p.add_argument("--af3-output-dir", default="af3_outputs",
                   help="AF3 output root for the structure phase (default: af3_outputs/).")
    p.add_argument("--stop-after", choices=STAGES, default=None,
                   help="Stop the pipeline after this stage (inclusive).")
    return p.parse_args(argv)


def _load_dbds(input_path: Path, subset: str | None,
               allow_fulllength: bool = False,
               allow_flagged: bool = False) -> tuple[dict[str, str], dict[str, dict]]:
    """Return ({gene: sequence}, {gene: meta}) where meta carries tf_family,
    zinc_finger_count, and a FoldX baseline ΔΔG if present in the TSV.

    Two classes of record are REFUSED rather than used (F1 in
    eCR_mod_lib/docs/DBD_EXTRACTION_FIXES.md), each with its own override
    because they are different decisions:

    - **No extracted span.** 692 of the library's 3,248 DBD records have no
      `dbd_sequence_aa` -- 262 where the extraction abstained, 430 never
      attempted -- and every one of them used to fall through to a full-length
      fusion here. That is the KRAB self-antagonism `06_extract_dbd.py` exists
      to prevent: the fib->beta run selected ZNF793 and ZNF850 as OPENERS while
      both carry native KRAB domains, so the construct opens and silences the
      same locus. `--allow-fulllength-dbd` restores the old behaviour.
    - **`dbd_review_flag == 1`.** 287 records where an excluded domain overlaps
      the kept span -- ZNF541 and TRERF1 carry ELM2 (NuRD recruiters, i.e.
      unintended closers). `--allow-flagged-dbd` uses them anyway.

    `dbd_dropped_review_flag` (F7, under-inclusion) is deliberately NOT a
    refusal: dropping a module is a biology question the design is meant to
    answer, not a construct that contradicts itself."""
    df = pd.read_csv(input_path, sep="\t", dtype=str)
    if "sequence_aa" not in df.columns:
        print("ERROR: --dbd-input needs a 'sequence_aa' column "
              "(re-run cli.py with --include-sequence).", file=sys.stderr)
        sys.exit(1)
    df = df.dropna(subset=["sequence_aa"])
    if subset:
        wanted = {s.strip() for s in subset.split(",")}
        df = df[df["gene_name"].isin(wanted)]

    # eCR_mod_lib stores the WHOLE UniProt protein in sequence_aa, so fusing on
    # that column builds a full-length-TF fusion carrying every native domain --
    # including KRAB on a construct meant to open chromatin. dbd_sequence_aa is
    # the extracted DNA-binding span (06_extract_dbd.py); prefer it, and say so
    # loudly when a record has none rather than silently reverting.
    has_span = "dbd_sequence_aa" in df.columns
    fell_back: list[str] = []      # no span: used full-length under the override
    skipped_nospan: list[str] = []  # no span: refused
    used_flagged: list[str] = []    # flagged span: used under the override
    skipped_flagged: list[str] = []  # flagged span: refused

    dbds: dict[str, str] = {}
    meta: dict[str, dict] = {}
    for _, row in df.iterrows():
        gene = str(row["gene_name"])
        if gene in dbds or gene in skipped_nospan or gene in skipped_flagged:
            continue
        span = str(row["dbd_sequence_aa"]) if has_span else ""
        flagged = pd.to_numeric(row.get("dbd_review_flag"), errors="coerce") == 1
        if span and span.lower() not in ("nan", "none", ""):
            if flagged and not allow_flagged:
                skipped_flagged.append(gene)
                continue
            if flagged:
                used_flagged.append(gene)
            dbds[gene] = span
        elif allow_fulllength:
            fell_back.append(gene)
            dbds[gene] = str(row["sequence_aa"])
        else:
            skipped_nospan.append(gene)
            continue
        baseline = pd.to_numeric(row.get("foldx_ddg_kcal_mol"), errors="coerce")
        meta[gene] = {
            "tf_family": row.get("tf_family", ""),
            "zinc_finger_count": row.get("zinc_finger_count"),
            "baseline_ddg": None if pd.isna(baseline) else float(baseline),
            "species": str(row.get("species", "") or ""),
            "dbd_extraction": row.get("dbd_extraction"),
            "dbd_review_flag": row.get("dbd_review_flag"),
            "dbd_zinc_finger_count": row.get("dbd_zinc_finger_count"),
        }

    def _names(genes: list[str]) -> str:
        shown = ", ".join(sorted(genes)[:10])
        return f"{shown}{' ...' if len(genes) > 10 else ''}"

    if not has_span and not allow_fulllength:
        print("ERROR: --dbd-input has no 'dbd_sequence_aa' column, so every "
              "record would be fused on its FULL-LENGTH protein -- every native "
              "domain of each TF, KRAB included. Run eCR_mod_lib "
              "scripts/06_extract_dbd.py and re-export with --include-sequence, "
              "or pass --allow-fulllength-dbd to accept that.", file=sys.stderr)
        sys.exit(1)
    if not has_span:
        print("WARNING: --dbd-input has no 'dbd_sequence_aa' column; fusing on "
              "FULL-LENGTH sequences (--allow-fulllength-dbd). These constructs "
              "carry every native domain of each TF.", file=sys.stderr)
    elif fell_back:
        print(f"WARNING: no extracted DBD span for {len(fell_back)} gene(s); "
              f"using the full-length sequence (--allow-fulllength-dbd) for: "
              f"{_names(fell_back)}", file=sys.stderr)

    # The run summary names both refusal classes, so a design run cannot
    # silently lose chassis or silently fuse a whole TF.
    if skipped_nospan:
        print(f"SKIPPED {len(skipped_nospan)} DBD(s) with no extracted span "
              f"(fusing them would build full-length-TF constructs); pass "
              f"--allow-fulllength-dbd to use them: {_names(skipped_nospan)}",
              file=sys.stderr)
    if skipped_flagged:
        print(f"SKIPPED {len(skipped_flagged)} DBD(s) with dbd_review_flag=1 "
              f"(an excluded domain lies inside the kept span); pass "
              f"--allow-flagged-dbd to use them: {_names(skipped_flagged)}",
              file=sys.stderr)
    if used_flagged:
        print(f"WARNING: using {len(used_flagged)} span(s) with "
              f"dbd_review_flag=1 (--allow-flagged-dbd): {_names(used_flagged)}",
              file=sys.stderr)
    print(f"      DBD chassis: {len(dbds)} used, "
          f"{len(skipped_nospan)} skipped (no span), "
          f"{len(skipped_flagged)} skipped (flagged).", file=sys.stderr)
    if not dbds and (skipped_nospan or skipped_flagged):
        print("ERROR: every DBD in --dbd-input was refused. Re-run the "
              "eCR_mod_lib extraction, narrow --dbd, or pass "
              "--allow-fulllength-dbd / --allow-flagged-dbd deliberately.",
              file=sys.stderr)
        sys.exit(1)

    return dbds, meta


def _detect_species(dbd_meta: dict[str, dict]) -> str:
    """Most common non-empty species across the input DBDs (auto-detect target)."""
    from collections import Counter
    counts = Counter(m.get("species", "") for m in dbd_meta.values() if m.get("species"))
    return counts.most_common(1)[0][0] if counts else ""


def _resolve_self_proteome(
    fusion_cfg: dict, species_flag: str | None, dbd_meta: dict[str, dict],
) -> tuple[str | None, str]:
    """Resolve the self-proteome FASTA. Priority: explicit config path >
    --species flag > species auto-detected from the input TSV. Returns
    (path_or_None, source_label)."""
    explicit = fusion_cfg.get("self_proteome", "")
    if explicit:
        return explicit, "config"
    species = species_flag or _detect_species(dbd_meta)
    if not species:
        return None, "none"
    from ecr_predictor.fusion.junction import resolve_proteome_path
    path = resolve_proteome_path(species, fusion_cfg.get("proteome_dir"))
    origin = "--species" if species_flag else "auto-detected"
    return str(path), f"{origin}: {species}"


def _load_eds(subset: str | None) -> dict[str, str]:
    """Load effector domains (type='ED') from the eCR_mod_lib library.

    Library ED names are not unique: a name can carry several distinct
    sequences (1,148 of 5,530 names do), so keying on name alone kept the first
    row reached and discarded 1,439 of 6,969 rows, with the survivor decided by
    row order. A name carrying one sequence keeps that name as its key; a name
    carrying several yields one key per sequence, suffixed with a short digest
    of that sequence ("NAME~a1b2c3"), so every distinct ED reaches the design
    run and the key is reproducible across runs and rebuilds. Rows that are
    exact duplicates (same name and same sequence) still collapse to one.

    `subset` (--ed) matches either a bare library name - selecting every
    sequence under it - or one suffixed key, selecting just that sequence.
    """
    from ecr_predictor.query import load_eds
    df = load_eds()
    df = df.dropna(subset=["sequence_aa"])
    name_col = "name" if "name" in df.columns else "gene_symbol"
    by_name: dict[str, list[str]] = {}
    for _, row in df.iterrows():
        seqs = by_name.setdefault(str(row[name_col]), [])
        seq = str(row["sequence_aa"])
        if seq not in seqs:
            seqs.append(seq)
    eds: dict[str, str] = {}
    for name, seqs in by_name.items():
        if len(seqs) == 1:
            eds[name] = seqs[0]
            continue
        for seq in seqs:
            digest = hashlib.sha256(seq.encode()).hexdigest()[:6]
            eds[f"{name}~{digest}"] = seq
    if subset:
        wanted = {s.strip() for s in subset.split(",")}
        eds = {k: v for k, v in eds.items()
               if k in wanted or k.split("~", 1)[0] in wanted}
    affected = sorted({k.split("~", 1)[0] for k in eds if "~" in k})
    if affected:
        shown = ", ".join(affected[:10])
        more = f" (+{len(affected) - 10} more)" if len(affected) > 10 else ""
        print(f"WARNING: {len(affected)} ED name(s) carry >1 distinct sequence "
              f"in the library; each sequence is kept under its own "
              f"'<name>~<digest>' key: {shown}{more}", file=sys.stderr)
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
    dbds, dbd_meta = _load_dbds(input_path, args.dbd,
                                allow_fulllength=args.allow_fulllength_dbd,
                                allow_flagged=args.allow_flagged_dbd)
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
    self_fa, source = _resolve_self_proteome(fusion_cfg, args.species, dbd_meta)
    if not self_fa:
        print("      Self proteome: none (no --species, no 'species' column, no "
              "config path) — Gate 1 self-tolerance filtering OFF.", file=sys.stderr)
    elif not Path(self_fa).exists():
        print(f"      WARNING: self proteome not found ({source}): {self_fa}\n"
              f"      Gate 1 self-tolerance filtering OFF. Download it there, e.g.:\n"
              f"        curl -o {self_fa}.gz 'https://rest.uniprot.org/uniprotkb/stream"
              f"?query=(organism_name:...)+AND+(reviewed:true)&format=fasta&compressed=true'"
              f" && gunzip {self_fa}.gz", file=sys.stderr)
    else:
        from ecr_predictor.fusion.junction import load_proteome_blob
        print(f"      Loading self proteome ({source}): {self_fa}", file=sys.stderr)
        proteome_blob = load_proteome_blob(self_fa)

    from ecr_predictor.fusion import backends
    immuno = stability = aggregation = structure = None

    # -------------------------------------------------------------------------
    # [2] Sequence gates (cheap; all candidates)
    # -------------------------------------------------------------------------
    print("\n[2/5] Sequence gates — immunogenicity + aggregation + stability...",
          file=sys.stderr)
    from ecr_predictor.fusion.immunogenicity import screen_immunogenicity
    from ecr_predictor.fusion.stability import screen_stability
    from ecr_predictor.fusion.aggregation import screen_aggregation
    try:
        immuno = screen_immunogenicity(candidates, fusion_cfg, proteome_blob)
    except (backends.ToolDisabled, backends.ToolNotAvailableError) as e:
        print(f"      SKIP Gate 1: {e}", file=sys.stderr)
    stability = screen_stability(candidates, fusion_cfg)
    # Gate 2 on sequence alone, over ALL candidates. It used to run only in the
    # structure phase, so the prune that decides who gets GPU time never saw the
    # aggregation axis — and with --stop-after prune it reported NA throughout.
    # The structure phase re-runs it with per-residue SASA, which is the better
    # measurement; this is the cheap version that can actually inform pruning.
    try:
        aggregation = screen_aggregation(candidates, fusion_cfg, structure_paths=None)
    except (backends.ToolDisabled, backends.ToolNotAvailableError) as e:
        print(f"      SKIP Gate 2: {e}", file=sys.stderr)

    if args.stop_after == "sequence":
        _write_scores(candidates, immuno, aggregation, stability, structure, Path(args.output))
        return

    # -------------------------------------------------------------------------
    # [3] Prune — Pareto over sequence axes → survivors
    # -------------------------------------------------------------------------
    print("\n[3/5] Pruning to Pareto-optimal survivors...", file=sys.stderr)
    from ecr_predictor.fusion.ranking import (
        build_score_table, add_pareto_and_rank, SEQUENCE_AXES,
    )
    seq_table = add_pareto_and_rank(
        build_score_table(candidates, immuno, aggregation, stability), axes=SEQUENCE_AXES,
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
    from ecr_predictor.fusion.ranking import build_score_table, add_pareto_and_rank
    df = build_score_table(candidates, immuno, aggregation, stability, structure)
    df = add_pareto_and_rank(df)
    path.write_text(df.to_csv(sep="\t", index=False, na_rep="NA"), encoding="utf-8")
    print(f"\nScores written to: {path}", file=sys.stderr)


if __name__ == "__main__":
    main()
