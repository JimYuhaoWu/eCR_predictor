"""
DBD off-target evaluation (Extension 1).

For a DBD with a JASPAR motif, scan a background of *accessible* genomic regions
for unintended strong binding sites and summarize them as an off-target effect
score. The background is the union of chromatin-accessibility peaks across the
reprogramming trajectory (e.g. MEF ∪ mES ATAC peaks) — a TF can only bind open
chromatin, so this is a far more meaningful denominator than the whole genome.

Two-tier severity weighting (score = Σ over off-target sites of strength × weight):

  Tier 1 (now)   region weight from accessibility DYNAMICS — regions that OPEN or
                 CLOSE across the transition are weighted above constitutively
                 accessible background. A provisional stand-in for |log2FC| /
                 driver importance, derived from the peak sets alone.
  Tier 2 (later) swap in per-region driver weights from eCR_navigator through the
                 same interface (a weight function over regions); the scoring code
                 does not change — only the weight source.

Reuses the PWM core from ecr_predictor.scan (fetch_pwm + best_site_score). Genome
sequence extraction uses pyfaidx against a local FASTA whose assembly MUST match
the peak coordinates (e.g. mm10). Region construction is pure-stdlib and needs no
genome, so the background set can be built/inspected independently.
"""
from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

# Tier-1 provisional weights by accessibility dynamics (see module docstring).
DEFAULT_DYNAMICS_WEIGHTS = {"opening": 1.0, "closing": 1.0, "shared": 0.3}
# Normalized log-odds ≥ this counts as a bound site (matches Step-1 motif_score scale).
DEFAULT_THRESHOLD = 0.8


@dataclass
class Region:
    chrom: str
    start: int
    end: int
    states: frozenset          # state labels whose peaks overlap this region
    dynamics: str = "shared"   # opening | closing | shared
    weight: float = 1.0

    @property
    def id(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


# ---------------------------------------------------------------------------
# Region set (pure stdlib — no genome needed)
# ---------------------------------------------------------------------------

def read_peak_bed(path: str | Path):
    """Yield (chrom, start, end) from a BED / MACS peak file (≥3 columns)."""
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                yield f[0], int(f[1]), int(f[2])
            except ValueError:
                continue


def _dynamics(states: frozenset, from_state: str, to_state: str) -> str:
    """opening = to-state only, closing = from-state only, shared = both/other."""
    has_from, has_to = from_state in states, to_state in states
    if has_from and has_to:
        return "shared"
    if has_to:
        return "opening"
    if has_from:
        return "closing"
    return "shared"


def build_region_union(
    state_beds: dict[str, str],
    from_state: str,
    to_state: str,
    dynamics_weights: dict[str, float] | None = None,
) -> list[Region]:
    """Merge peaks from labeled state BED files into a union region set.

    Each merged region is tagged with the states that contributed to it and a
    dynamics label (opening / closing / shared), and given the corresponding
    Tier-1 weight. Overlapping/bookended intervals are merged."""
    weights = dynamics_weights or DEFAULT_DYNAMICS_WEIGHTS

    tagged: list[tuple[str, int, int, str]] = []
    for state, path in state_beds.items():
        for chrom, start, end in read_peak_bed(path):
            tagged.append((chrom, start, end, state))
    tagged.sort(key=lambda t: (t[0], t[1], t[2]))

    regions: list[Region] = []
    cur_chrom: str | None = None
    cur_start = cur_end = 0
    cur_states: set[str] = set()

    def _flush() -> None:
        if cur_chrom is None:
            return
        st = frozenset(cur_states)
        dyn = _dynamics(st, from_state, to_state)
        regions.append(Region(cur_chrom, cur_start, cur_end, st, dyn, weights.get(dyn, 1.0)))

    for chrom, start, end, state in tagged:
        if chrom != cur_chrom or start > cur_end:
            _flush()
            cur_chrom, cur_start, cur_end, cur_states = chrom, start, end, {state}
        else:
            cur_end = max(cur_end, end)
            cur_states.add(state)
    _flush()
    return regions


def region_stats(regions: list[Region]) -> dict[str, int]:
    """Count regions by dynamics label (+ total)."""
    stats = {"total": len(regions), "opening": 0, "closing": 0, "shared": 0}
    for r in regions:
        stats[r.dynamics] = stats.get(r.dynamics, 0) + 1
    return stats


def write_region_bed(regions: list[Region], path: str | Path) -> None:
    """Write the union region set as BED6-ish: chrom start end dynamics weight."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        for r in regions:
            w.writerow([r.chrom, r.start, r.end, r.dynamics, r.weight])


# ---------------------------------------------------------------------------
# Genome sequence extraction (server — needs pyfaidx + a matching FASTA)
# ---------------------------------------------------------------------------

def extract_sequences(regions: list[Region], genome_fasta: str | Path) -> dict[str, str]:
    """{region.id: uppercase DNA} from a genome FASTA. The FASTA assembly MUST
    match the peak coordinates (e.g. mm10)."""
    from pyfaidx import Fasta
    fa = Fasta(str(genome_fasta), sequence_always_upper=True, as_raw=True)
    seqs: dict[str, str] = {}
    missing_chrom = 0
    for r in regions:
        if r.chrom not in fa:
            missing_chrom += 1
            continue
        try:
            seqs[r.id] = str(fa[r.chrom][r.start:r.end])
        except Exception:
            missing_chrom += 1
    if missing_chrom:
        print(f"  WARNING: {missing_chrom} regions had no sequence "
              f"(chromosome name / assembly mismatch?).", file=sys.stderr)
    return seqs


# ---------------------------------------------------------------------------
# Off-target scoring (reuses ecr_predictor.scan PWM core)
# ---------------------------------------------------------------------------

def score_offtarget(
    jaspar_id: str,
    regions: list[Region],
    sequences: dict[str, str],
    threshold: float = DEFAULT_THRESHOLD,
    weight_fn=None,
) -> dict | None:
    """Off-target effect for one motif over the region set.

    Returns {n_offtarget_sites, offtarget_effect_score, offtarget_unweighted} or
    None if the motif can't be fetched. `weight_fn(region) -> float` selects the
    severity weight (Tier 1 default: the region's dynamics weight; Tier 2: pass an
    eCR_navigator-derived weight function — no other change needed)."""
    from ecr_predictor.scan import fetch_pwm, best_site_score
    motif = fetch_pwm(jaspar_id)
    if motif is None:
        return None
    if weight_fn is None:
        weight_fn = lambda r: r.weight

    m = len(motif)
    n_sites = 0
    weighted = 0.0
    unweighted = 0.0
    for r in regions:
        seq = sequences.get(r.id)
        if not seq or len(seq) < m:
            continue
        s = best_site_score(motif, seq)
        if s != s:  # NaN
            continue
        if s >= threshold:
            n_sites += 1
            unweighted += s
            weighted += s * weight_fn(r)
    return {
        "n_offtarget_sites": n_sites,
        "offtarget_effect_score": round(weighted, 4),
        "offtarget_unweighted": round(unweighted, 4),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_peaks(items: list[str]) -> dict[str, str]:
    """['MEF=/path/a.bed', 'mES=/path/b.bed'] -> {'MEF': ..., 'mES': ...}."""
    peaks: dict[str, str] = {}
    for it in items:
        if "=" not in it:
            print(f"ERROR: --peaks entry must be STATE=path, got {it!r}", file=sys.stderr)
            sys.exit(1)
        k, v = it.split("=", 1)
        peaks[k.strip()] = v.strip()
    return peaks


def _read_dbd_input(path: str) -> list[tuple[str, str]]:
    """Return [(gene_name, jaspar_id)] from a Step-1/2 TSV."""
    rows: list[tuple[str, str]] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = row.get("gene_name") or row.get("gene_symbol") or "?"
            jid = (row.get("jaspar_id") or "").strip()
            rows.append((gene, jid))
    return rows


def _parse_args(argv):
    p = argparse.ArgumentParser(
        description="DBD off-target evaluation over accessible genomic regions."
    )
    p.add_argument("--dbd-input",
                   help="Step-1/2 TSV with gene_name + jaspar_id (required unless "
                        "--stop-after-regions).")
    p.add_argument("--peaks", nargs="+", required=True, metavar="STATE=BED",
                   help="Labeled peak BEDs, e.g. --peaks MEF=mef.bed mES=mes.bed")
    p.add_argument("--from", dest="from_state", required=True,
                   help="State label of the starting cell type (e.g. MEF).")
    p.add_argument("--to", dest="to_state", required=True,
                   help="State label of the target cell type (e.g. mES).")
    p.add_argument("--genome", help="Genome FASTA matching the peak assembly (e.g. mm10.fa).")
    p.add_argument("--output", default="offtarget.tsv", help="Output TSV path.")
    p.add_argument("--regions-out", help="Optional path to dump the union region BED.")
    p.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD,
                   help=f"Normalized log-odds bound-site threshold (default {DEFAULT_THRESHOLD}).")
    p.add_argument("--stop-after-regions", action="store_true",
                   help="Build/inspect the region set and exit (no genome/scan).")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    peaks = _parse_peaks(args.peaks)
    for label in (args.from_state, args.to_state):
        if label not in peaks:
            print(f"ERROR: --from/--to label {label!r} not among --peaks states "
                  f"{list(peaks)}.", file=sys.stderr)
            sys.exit(1)

    print(f"[offtarget] Building union region set from {list(peaks)} "
          f"(from={args.from_state} → to={args.to_state})...", file=sys.stderr)
    regions = build_region_union(peaks, args.from_state, args.to_state)
    st = region_stats(regions)
    print(f"  {st['total']} regions — opening {st['opening']}, closing "
          f"{st['closing']}, shared {st['shared']}.", file=sys.stderr)

    if args.regions_out:
        write_region_bed(regions, args.regions_out)
        print(f"  Region set written to {args.regions_out}", file=sys.stderr)

    if args.stop_after_regions:
        return
    if not args.genome or not args.dbd_input:
        print("ERROR: --genome and --dbd-input are required to score off-target "
              "(or use --stop-after-regions).", file=sys.stderr)
        sys.exit(1)

    print(f"[offtarget] Extracting region sequences from {args.genome}...", file=sys.stderr)
    sequences = extract_sequences(regions, args.genome)

    dbds = _read_dbd_input(args.dbd_input)
    print(f"[offtarget] Scoring {len(dbds)} DBDs (threshold {args.threshold})...",
          file=sys.stderr)

    cache: dict[str, dict | None] = {}
    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene_name", "jaspar_id", "n_offtarget_sites",
                         "offtarget_effect_score", "offtarget_unweighted"])
        for gene, jid in dbds:
            if not jid:
                writer.writerow([gene, "", "NA", "NA", "NA"])
                continue
            if jid not in cache:
                cache[jid] = score_offtarget(jid, regions, sequences, args.threshold)
            res = cache[jid]
            if res is None:
                writer.writerow([gene, jid, "NA", "NA", "NA"])
            else:
                writer.writerow([gene, jid, res["n_offtarget_sites"],
                                 res["offtarget_effect_score"], res["offtarget_unweighted"]])
            print(f"  {gene} ({jid}): {res}", file=sys.stderr)

    print(f"[offtarget] Written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
