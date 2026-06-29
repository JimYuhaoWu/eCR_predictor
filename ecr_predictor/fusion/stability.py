"""
Gate 3 — intracellular stability liabilities.

For intracellularly expressed fusions, half-life is governed by the
ubiquitin-proteasome system. Two sequence-level checks run without any external
tool (fully local logic), plus an optional ubiquitination predictor via backend:

  - N-end rule    destabilizing N-terminal residue after Met excision
  - degron motifs short linear motifs that recruit E3 ligases (e.g. phospho-
                  degrons, KEN/D-box) appearing in/around the linker
  - ubpred        (optional) per-Lys ubiquitination probability via backend

CAVEAT: proteasomal degradation also FEEDS MHC-I presentation (Gate 1), so
"stabilize" and "reduce immunogenicity" can conflict. This gate reports
liabilities; it does not auto-resolve that tension.
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field

from ecr_predictor.fusion import backends
from ecr_predictor.fusion.assemble import FusionCandidate

# Type-2 (basic) and type-1 destabilizing N-terminal residues (mammalian N-end rule).
DESTABILIZING_NTERM = set("RKHFWYLIDE")

# Minimal degron motif set (illustrative; extend per target E3 repertoire).
DEGRON_MOTIFS = {
    "D-box": r"R..L..[LIVM].",
    "KEN-box": r"KEN",
    "phospho-degron(DSG)": r"DSG..S",
}


@dataclass
class StabilityResult:
    nterm_destabilizing: bool
    nterm_residue: str
    degron_hits: list[str] = field(default_factory=list)
    ubiquitination_lys: list[int] = field(default_factory=list)  # flagged Lys positions


def _nterm_residue(seq: str) -> str:
    """Effective N-terminal residue after initiator-Met excision when applicable.
    Met is excised when the second residue is small (A,C,G,P,S,T,V)."""
    if not seq:
        return ""
    if seq[0] == "M" and len(seq) > 1 and seq[1] in set("ACGPSTV"):
        return seq[1]
    return seq[0]


def _degron_hits(seq: str) -> list[str]:
    return [name for name, pat in DEGRON_MOTIFS.items() if re.search(pat, seq)]


def _ubpred_local(seq: str, cfg: dict, rank_threshold: float) -> list[int]:
    """Optional: run a local ubiquitination predictor, return flagged Lys positions."""
    local = cfg.get("local", {})
    command = local.get("command", "ubpred")
    extra = local.get("extra_args", "").split()
    fasta = f">q\n{seq}\n"
    stdout = backends.run_cli(
        command, ["-f", "{fasta}", *extra], input_files={"fasta": fasta},
        timeout=int(local.get("timeout", 300)),
    )
    flagged: list[int] = []
    for line in stdout.splitlines():
        cols = line.split()
        if len(cols) < 2:
            continue
        try:
            pos, score = int(cols[0]), float(cols[-1])
        except ValueError:
            continue
        if score >= rank_threshold:
            flagged.append(pos)
    return flagged


def screen_candidate(candidate: FusionCandidate, fusion_cfg: dict) -> StabilityResult:
    seq = candidate.sequence
    nterm = _nterm_residue(seq)

    ub_lys: list[int] = []
    ub_cfg = fusion_cfg.get("tools", {}).get("ubpred", {})
    if backends.resolve_backend(ub_cfg) != "disabled":
        thr = float(ub_cfg.get("score_threshold", 0.6))
        try:
            if backends.resolve_backend(ub_cfg) == "local":
                ub_lys = _ubpred_local(seq, ub_cfg, thr)
            # (API path intentionally omitted until a concrete UbPred service exists)
        except backends.ToolNotAvailableError as e:
            print(f"  WARNING: ubpred unavailable, skipping: {e}", file=sys.stderr)

    return StabilityResult(
        nterm_destabilizing=nterm in DESTABILIZING_NTERM,
        nterm_residue=nterm,
        degron_hits=_degron_hits(seq),
        ubiquitination_lys=ub_lys,
    )


def screen_stability(
    candidates: list[FusionCandidate],
    fusion_cfg: dict,
) -> dict[str, StabilityResult]:
    """Run Gate 3 across all candidates. {name: result}. Always available
    (N-end rule + degron scan need no external tool)."""
    print("[Gate3] Stability (N-end rule + degron scan)", file=sys.stderr)
    results: dict[str, StabilityResult] = {}
    for cand in candidates:
        res = screen_candidate(cand, fusion_cfg)
        results[cand.name] = res
        print(f"  {cand.name}: Nterm={res.nterm_residue} "
              f"destab={res.nterm_destabilizing} degrons={res.degron_hits}",
              file=sys.stderr)
    return results
