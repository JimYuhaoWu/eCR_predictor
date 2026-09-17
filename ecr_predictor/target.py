"""
Nominated region -> AF3-safe target site: the on-target bridge.

eCR_navigator nominates REGIONS (hundreds of bp; see docs/run_bundle_contract.md
in that repo), but refine.py folds a base-paired duplex (af3.py chain B sense /
chain C revcomp), so the design target must be a SHORT site (~30 bp). Nothing
converted one to the other; this module does.

  nominations.tsv
    | extract region sequence     (pyfaidx, manifest.transition.assembly)
    | scan library DBD motifs     [scan.fetch_pwm + best_site_hit]
    | best window per DBD, flanked to an AF3-safe length
    | rank by on-target strength x specificity
    v
  sites.tsv -> cli.py --sequence <site> -> refine.py -> fuse.py

This is offtarget.py inverted, and deliberately shares its machinery: maximize
the PWM strength of the INTENDED site while minimizing the off-target burden
offtarget.py already computes over the accessible background. That burden is a
property of the MOTIF, not of the candidate window, so it is computed once by
offtarget.py (ideally with --weights, making it Tier-2 severity-weighted by the
bundle's own driver scores) and joined here on jaspar_id. No scoring code is
duplicated, and --offtarget stays optional: the library-gap report below needs
on-target scanning only.

LIBRARY-GAP REPORT (--gaps-out). Where NO library DBD clears the motif threshold,
the region is undesignable with today's library. Those rows, with their best near
miss, are the direct eCR_mod_lib expansion signal. Expect this to be the larger
output: only 9 of 38 library DBDs carry a JASPAR motif (8 human, 1 mouse).

Sites are scanned on BOTH strands, because a duplex has two and AF3 folds both
(af3.py chain B sense / chain C revcomp); `strand` records which one the motif
reads on. `sequence` is always the + strand reference at [site_start, site_end),
so it can be checked against the genome directly.

Site selection is LIBRARY-FIRST (best PWM hit per library DBD, then filtered by
off-target burden) rather than specificity-first, which would need a genome-wide
scan per candidate window; the library is small enough to scan exhaustively.
DBDs are species-matched to the transition like cli.py, because the eCR is
expressed intracellularly -- a cross-species DBD makes the whole domain foreign,
a worse version of the junction liability fuse.py already screens for.

Two contract rules, both load-bearing:
  * ORDER BY `rank`, never threshold on `nomination_score`. The score is a
    percentile (manifest nomination.score_norm = "rank"), so a top-1% band spans
    0.99-1.0 by construction -- the iN fixture has 17 exact ties at 1.0. `rank`
    is the only column carrying prioritization signal; --top-n is the budget.
  * SMALL |direction| IS LOW CONFIDENCE. `direction` picks the ED (+ activator,
    - repressor) and flows to fuse.py, but 42% of the iN fixture falls below
    0.05. Such rows are LABELLED `ambiguous` and never dropped: the direction
    column itself is unvalidated (Claim 2B, deferred), so discarding half a
    validated band on its say-so would destroy more than it protects.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path

from ecr_predictor.offtarget import DEFAULT_THRESHOLD

# Sites are flanked to this length before AF3 folds them as a duplex. Matches the
# query length in the development workflow (CLAUDE.md).
DEFAULT_SITE_LENGTH = 30
# |direction| below this is too weak to call an ED. Provisional: it is the
# inflection reported for the iN fixture, and v1 of the contract deliberately
# sets no threshold. It LABELS rows; it never drops them.
DEFAULT_DIRECTION_MIN = 0.05
# Design budget: how many nominated regions (by rank) to emit sites for.
DEFAULT_TOP_N = 50

# manifest.transition.species -> the library's `organism` vocabulary.
_SPECIES_BY_LABEL = {"human": "Homo sapiens", "mouse": "Mus musculus"}

SITE_COLUMNS = [
    "site_rank", "region_id", "nomination_rank", "nomination_score",
    "chrom", "site_start", "site_end", "sequence", "strand",
    "motif_start", "motif_end",
    "gene_name", "module_id", "jaspar_id", "motif_score",
    "n_offtarget_sites", "offtarget_effect_score", "specificity",
    "design_score", "pareto_optimal",
    "direction", "ed_class",
]

GAP_COLUMNS = [
    "chrom", "start", "end", "nomination_rank", "nomination_score", "direction",
    "n_dbds_scanned", "best_gene_name", "best_jaspar_id", "best_motif_score",
]


@dataclass
class Nomination:
    chrom: str
    start: int
    end: int
    rank: int
    nomination_score: float
    direction: float | None   # None = unmeasured; 0.0 = measured-flat. Distinct.

    @property
    def id(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


@dataclass
class Site:
    nom: Nomination
    gene_name: str
    module_id: str
    jaspar_id: str
    motif_score: float
    strand: str
    site_start: int
    site_end: int
    sequence: str
    motif_start: int
    motif_end: int
    n_offtarget_sites: int | None = None
    burden: float | None = None
    specificity: float | None = None
    design_score: float | None = None
    pareto_optimal: bool = False


# ---------------------------------------------------------------------------
# Bundle I/O (eCR_navigator run_bundle_contract.md v1)
# ---------------------------------------------------------------------------

def read_manifest(bundle: str | Path) -> dict:
    with open(Path(bundle) / "manifest.json") as fh:
        return json.load(fh)


def read_nominations(bundle: str | Path) -> list[Nomination]:
    """Read nominations.tsv in rank order. A refusal is a header line and
    nothing else — that is a valid bundle and yields an empty list."""
    out: list[Nomination] = []
    with open(Path(bundle) / "nominations.tsv", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            raw_dir = row.get("direction")
            out.append(Nomination(
                chrom=row["chrom"],
                start=int(row["start"]),
                end=int(row["end"]),
                rank=int(row["rank"]),
                nomination_score=float(row["nomination_score"]),
                direction=float(raw_dir) if raw_dir not in (None, "") else None,
            ))
    out.sort(key=lambda n: n.rank)
    return out


def read_offtarget(path: str | Path) -> dict[str, tuple[int, float]]:
    """{jaspar_id: (n_offtarget_sites, offtarget_effect_score)} from an
    offtarget.py output TSV. Rows it could not score are written as NA and are
    skipped here, leaving that DBD's burden unknown."""
    out: dict[str, tuple[int, float]] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            jid = (row.get("jaspar_id") or "").strip()
            if not jid or row.get("offtarget_effect_score") in (None, "", "NA"):
                continue
            out[jid] = (int(row["n_offtarget_sites"]),
                        float(row["offtarget_effect_score"]))
    return out


# ---------------------------------------------------------------------------
# Site construction
# ---------------------------------------------------------------------------

def check_assembly(noms: list[Nomination], chrom_lengths: dict[str, int],
                   assembly: str) -> None:
    """Fail loudly on an assembly mismatch between the bundle and the FASTA.

    manifest.transition.assembly is load-bearing: mismatched coordinates
    silently corrupt site selection, so this is a hard error, not a warning."""
    missing = sorted({n.chrom for n in noms if n.chrom not in chrom_lengths})
    if missing:
        raise SystemExit(
            f"ERROR: {len(missing)} nominated chromosome(s) absent from the genome "
            f"FASTA: {', '.join(missing[:5])}. The bundle is {assembly} — is the "
            f"FASTA?"
        )
    over = [n for n in noms if n.end > chrom_lengths[n.chrom]]
    if over:
        n = over[0]
        raise SystemExit(
            f"ERROR: {len(over)} nominated region(s) run past the end of their "
            f"chromosome (e.g. {n.id}, {n.chrom} is {chrom_lengths[n.chrom]} bp). "
            f"The bundle is {assembly} — the FASTA is a different assembly."
        )


def site_window(region_start: int, hit_offset: int, motif_len: int,
                site_length: int, chrom_length: int) -> tuple[int, int]:
    """Centre an AF3-safe window on the motif hit, clamped to the chromosome.

    Flanks are real genomic sequence and may extend past the nominated region's
    edge — a hit near the boundary gets true context rather than a truncated
    site. Windows never shrink below the motif itself."""
    site_len = min(max(site_length, motif_len), chrom_length)
    hit_start = region_start + hit_offset
    start = hit_start - (site_len - motif_len) // 2
    start = max(0, min(start, chrom_length - site_len))
    return start, start + site_len


def ed_class(direction: float | None, direction_min: float) -> str:
    """Which effector domain the direction supports: activator (+) / repressor
    (-) / ambiguous (too weak, or measured-flat) / unmeasured (empty)."""
    if direction is None:
        return "unmeasured"
    if direction == 0.0 or abs(direction) < direction_min:
        return "ambiguous"
    return "activator" if direction > 0 else "repressor"


def compute_specificity(burdens: dict[str, float],
                        scale: float | None = None) -> tuple[dict[str, float], float]:
    """({jaspar_id: specificity in (0,1]}, scale) — 1/(1 + burden/scale).

    Specificity is library-relative, not absolute: it says which of your DBDs is
    cleaner, not that any of them is clean. Scale-free (the burden's absolute
    units cancel) and never degenerate — the worst DBD in the library still gets
    a positive weight rather than being annihilated in design_score.

    `scale` defaults to the median positive burden across the DBDs in THIS run,
    which makes design_score incomparable between runs whenever the library
    changes: adding motifs moves the median, so every existing DBD's specificity
    and design_score shift and the Pareto front reshuffles even though no
    burden moved. Pass an explicit `scale` (echoed by every run, and recorded in
    the output header) to hold rankings comparable across library versions.
    Returns the scale actually used so the caller can record it."""
    positive = [b for b in burdens.values() if b > 0]
    if scale is None:
        if not positive:
            return {jid: 1.0 for jid in burdens}, 0.0
        scale = statistics.median(positive)
    if scale <= 0:
        return {jid: 1.0 for jid in burdens}, scale
    return {jid: 1.0 / (1.0 + b / scale) for jid, b in burdens.items()}, scale


def mark_pareto(sites: list[Site]) -> None:
    """Flag the Pareto front over (motif_score up, burden down).

    design_score is a convenience sort; this is the honest ranking — it keeps the
    two axes independent instead of letting the product hide which one carried a
    row. Same pattern as fusion/ranking.py. Sites with unknown burden cannot be
    placed on the front."""
    known = [s for s in sites if s.burden is not None]
    for a in known:
        a.pareto_optimal = not any(
            b is not a
            and b.motif_score >= a.motif_score and b.burden <= a.burden
            and (b.motif_score > a.motif_score or b.burden < a.burden)
            for b in known
        )


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _fmt(value, nd: int | None = None) -> str:
    if value is None:
        return "NA"
    return str(round(value, nd)) if nd is not None else str(value)


def write_sites(sites: list[Site], path: str | Path, direction_min: float) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(SITE_COLUMNS)
        for i, s in enumerate(sites, 1):
            n = s.nom
            w.writerow([
                i, n.id, n.rank, n.nomination_score,
                n.chrom, s.site_start, s.site_end, s.sequence, s.strand,
                s.motif_start, s.motif_end,
                s.gene_name, s.module_id, s.jaspar_id, round(s.motif_score, 4),
                _fmt(s.n_offtarget_sites), _fmt(s.burden, 4), _fmt(s.specificity, 4),
                _fmt(s.design_score, 4), s.pareto_optimal,
                "" if n.direction is None else n.direction,
                ed_class(n.direction, direction_min),
            ])


def write_gaps(gaps: list[tuple], path: str | Path) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(GAP_COLUMNS)
        for nom, n_scanned, best in sorted(
            gaps, key=lambda g: (-(g[2][2] if g[2] else -2), g[0].rank)
        ):
            w.writerow([
                nom.chrom, nom.start, nom.end, nom.rank, nom.nomination_score,
                "" if nom.direction is None else nom.direction,
                n_scanned,
                best[0] if best else "NA",
                best[1] if best else "NA",
                round(best[2], 4) if best else "NA",
            ])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _load_library(species: str, db_path: Path | None):
    """Motif-bearing library DBDs for `species` (exact, then genus — as cli.py).

    A DBD with no jaspar_id cannot be placed on a locus by sequence, so unlike
    cli.py (which keeps them with motif_score = NA) target.py must drop them."""
    from ecr_predictor.query import load_dbds, match_species
    dbds = load_dbds(db_path) if db_path else load_dbds()
    candidates = match_species(dbds, species)
    if candidates.empty:
        raise SystemExit(f"ERROR: no DBDs in the library for species {species!r}.")
    scannable = candidates[
        candidates["jaspar_id"].notna() & (candidates["jaspar_id"].str.strip() != "")
    ]
    return scannable


def _parse_args(argv):
    p = argparse.ArgumentParser(
        description="Nominated regions -> AF3-safe target sites for cli.py --sequence."
    )
    p.add_argument("--bundle", required=True,
                   help="eCR_navigator run-bundle directory (manifest.json + "
                        "nominations.tsv).")
    p.add_argument("--genome", required=True,
                   help="Genome FASTA matching manifest.transition.assembly.")
    p.add_argument("--offtarget",
                   help="offtarget.py output TSV, joined on jaspar_id for the "
                        "specificity axis (run it with --weights for Tier-2). "
                        "Omit to rank on on-target strength alone — the "
                        "library-gap report does not need it.")
    p.add_argument("--output", default="sites.tsv", help="Output sites TSV.")
    p.add_argument("--gaps-out",
                   help="Library-gap report: nominated regions no library DBD can "
                        "target (the eCR_mod_lib expansion signal).")
    p.add_argument("--top-n", type=int, default=DEFAULT_TOP_N,
                   help=f"Design budget: emit sites for the top N regions by rank "
                        f"(default {DEFAULT_TOP_N}). All nominations are still "
                        f"scanned for the gap report.")
    p.add_argument("--site-length", type=int, default=DEFAULT_SITE_LENGTH,
                   help=f"AF3 duplex length in bp (default {DEFAULT_SITE_LENGTH}).")
    p.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD,
                   help=f"Normalized log-odds bound-site threshold (default "
                        f"{DEFAULT_THRESHOLD}); same scale as motif_score.")
    p.add_argument("--specificity-scale", type=float,
                   help="Burden scale for specificity = 1/(1 + burden/scale). "
                        "Default: the median positive burden in THIS run, which "
                        "shifts every design_score whenever the library changes. "
                        "Pin it to the value a previous run reported to keep "
                        "rankings comparable across library versions.")
    p.add_argument("--direction-min", type=float, default=DEFAULT_DIRECTION_MIN,
                   help=f"|direction| below this is labelled ed_class=ambiguous "
                        f"(default {DEFAULT_DIRECTION_MIN}). Labels only — no row "
                        f"is ever dropped.")
    p.add_argument("--species",
                   help="Override the library species (default: from "
                        "manifest.transition.species).")
    p.add_argument("--db", help="Path to module_library.db (default: auto-detected).")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    manifest = read_manifest(args.bundle)
    transition = manifest["transition"]
    assembly = transition["assembly"]

    print(f"[target] Bundle {manifest['run_id']}: {transition['label']} "
          f"({assembly}).", file=sys.stderr)

    # Refusal is a first-class output, not an error: a Gate-1 reject emits a valid
    # bundle with zero nominations. Emit an empty sites.tsv and exit clean.
    nomination = manifest.get("nomination", {})
    if nomination.get("refused"):
        print(f"[target] REFUSED — {nomination.get('refusal_reason')}", file=sys.stderr)
        print("[target] No trustworthy nomination exists; no sites designed.",
              file=sys.stderr)
        write_sites([], args.output, args.direction_min)
        print(f"[target] Empty site table written to {args.output}", file=sys.stderr)
        return

    noms = read_nominations(args.bundle)
    if not noms:
        raise SystemExit("ERROR: bundle is not marked refused but nominations.tsv "
                         "is empty — the bundle is inconsistent.")
    print(f"[target] {len(noms)} nominations, ordered by rank "
          f"(score_norm={nomination.get('score_norm')}; nomination_score is a "
          f"percentile — not thresholded).", file=sys.stderr)

    species = args.species or _SPECIES_BY_LABEL.get(transition.get("species"))
    if not species:
        raise SystemExit(
            f"ERROR: cannot map manifest species {transition.get('species')!r} to the "
            f"library's organism vocabulary. Pass --species explicitly."
        )
    dbds = _load_library(species, Path(args.db) if args.db else None)
    if dbds.empty:
        raise SystemExit(
            f"ERROR: no {species} DBD in the library carries a jaspar_id, so no "
            f"locus can be targeted by sequence. Every nomination is a library gap."
        )
    print(f"[target] {len(dbds)} motif-bearing {species} DBDs: "
          f"{', '.join(dbds['gene_symbol'].tolist())}", file=sys.stderr)

    from pyfaidx import Fasta
    fa = Fasta(str(args.genome), sequence_always_upper=True, as_raw=True)
    chrom_lengths = {name: len(fa[name]) for name in fa.keys()}
    check_assembly(noms, chrom_lengths, assembly)

    from ecr_predictor.scan import best_site_hit, fetch_pwm
    motifs = {}
    for jid in sorted(set(dbds["jaspar_id"].str.strip())):
        motif = fetch_pwm(jid)
        if motif is None:
            print(f"  WARNING: {jid} unavailable — its DBD is skipped.", file=sys.stderr)
            continue
        motifs[jid] = motif
    if not motifs:
        raise SystemExit("ERROR: no JASPAR motif could be fetched; nothing to scan.")

    burdens = read_offtarget(args.offtarget) if args.offtarget else {}
    if args.offtarget:
        print(f"[target] Off-target burden joined for {len(burdens)} motif(s) from "
              f"{args.offtarget}.", file=sys.stderr)
    else:
        print("[target] No --offtarget: ranking on on-target strength alone "
              "(specificity is NA).", file=sys.stderr)

    print(f"[target] Scanning {len(noms)} regions x {len(motifs)} motifs "
          f"(threshold {args.threshold})...", file=sys.stderr)

    rows = list(dbds.itertuples())
    sites: list[Site] = []
    gaps: list[tuple] = []
    n_designable = 0

    for nom in noms:
        seq = str(fa[nom.chrom][nom.start:nom.end])
        best_overall: tuple[str, str, float] | None = None
        hits = []
        for dbd in rows:
            jid = dbd.jaspar_id.strip()
            motif = motifs.get(jid)
            if motif is None:
                continue
            hit = best_site_hit(motif, seq)
            if hit is None:
                continue
            score, offset, strand = hit
            if best_overall is None or score > best_overall[2]:
                best_overall = (dbd.gene_symbol, jid, score)
            if score >= args.threshold:
                hits.append((dbd, jid, motif, score, offset, strand))

        if not hits:
            gaps.append((nom, len(motifs), best_overall))
            continue

        n_designable += 1
        if nom.rank > args.top_n:
            continue

        for dbd, jid, motif, score, offset, strand in hits:
            s_start, s_end = site_window(nom.start, offset, len(motif),
                                         args.site_length, chrom_lengths[nom.chrom])
            hit_start = nom.start + offset
            n_off, burden = burdens.get(jid, (None, None))
            sites.append(Site(
                nom=nom, gene_name=dbd.gene_symbol, module_id=dbd.module_id,
                jaspar_id=jid, motif_score=score, strand=strand,
                site_start=s_start, site_end=s_end,
                sequence=str(fa[nom.chrom][s_start:s_end]),
                motif_start=hit_start, motif_end=hit_start + len(motif),
                n_offtarget_sites=n_off, burden=burden,
            ))

    spec, spec_scale = compute_specificity(
        {j: b for j, (_, b) in burdens.items()}, args.specificity_scale)
    if burdens:
        pinned = "pinned" if args.specificity_scale is not None else "from this run"
        print(f"[target] Specificity scale ({pinned}): {spec_scale:.6g}. Pass "
              f"--specificity-scale {spec_scale:.6g} to keep later runs comparable.",
              file=sys.stderr)
    for s in sites:
        s.specificity = spec.get(s.jaspar_id)
        s.design_score = (s.motif_score * s.specificity
                          if s.specificity is not None else s.motif_score)
    mark_pareto(sites)
    sites.sort(key=lambda s: (-s.design_score, s.nom.rank))

    write_sites(sites, args.output, args.direction_min)
    pct = 100.0 * n_designable / len(noms)
    print(f"[target] {n_designable}/{len(noms)} regions designable with today's "
          f"library ({pct:.1f}%); {len(gaps)} gaps.", file=sys.stderr)
    print(f"[target] {len(sites)} sites (top-{args.top_n} by rank) -> {args.output}",
          file=sys.stderr)

    if args.gaps_out:
        write_gaps(gaps, args.gaps_out)
        print(f"[target] Library-gap report -> {args.gaps_out}", file=sys.stderr)
