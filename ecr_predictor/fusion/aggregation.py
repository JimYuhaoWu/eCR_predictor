"""
Gate 2 — aggregation / solubility liability at the junction.

Intracellular aggregation drives inclusion-body formation and reduced expression.
Two complementary predictors, each via the local-CLI-or-API backend:

  - aggrescan3d   structure-based aggregation propensity (needs a 3D model;
                  feed the AF3 fusion structure from Gate 0)
  - camsol        sequence-based intrinsic solubility

A candidate is flagged when an aggregation hotspot (residue score beyond the
tool's cutoff) falls within `junction_window` residues of any domain boundary —
the junction is where novel, potentially aggregation-prone sequence appears.

NOTE: per-residue output formats differ by tool/version. `_parse_*` below define
the expected (residue_index, score) contract; confirm against your install.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass

from ecr_predictor.fusion import backends
from ecr_predictor.fusion.assemble import FusionCandidate

DEFAULT_JUNCTION_WINDOW = 5


@dataclass
class AggregationResult:
    tool: str
    junction_hotspot: bool          # any hotspot within the junction window
    max_junction_score: float       # worst per-residue score near a junction
    n_hotspots: int                 # total hotspots in the candidate


def _residue_scores_local(
    seq: str, name: str, cfg: dict, structure_path: str | None = None,
) -> list[tuple[int, float]]:
    """Run a local aggregation tool, return [(residue_1based, score), ...].

    Structure-based tools (AGGRESCAN3D) take a PDB/CIF via --structure; sequence
    tools (CamSol) take a FASTA via -f."""
    local = cfg.get("local", {})
    command = local.get("command", name)
    extra = local.get("extra_args", "").split()
    timeout = int(local.get("timeout", 600))

    if structure_path:
        stdout = backends.run_cli(
            command, ["--structure", structure_path, *extra], timeout=timeout,
        )
    else:
        fasta = f">{name}\n{seq}\n"
        stdout = backends.run_cli(
            command, ["-f", "{fasta}", *extra], input_files={"fasta": fasta},
            timeout=timeout,
        )
    return _parse_residue_scores(stdout)


def _residue_scores_api(seq: str, cfg: dict) -> list[tuple[int, float]]:
    api = cfg.get("api", {})
    result = backends.submit_and_poll(
        api,
        {"sequence": seq},
        parse_job_id=lambda d: d.get("id") or d.get("job_id"),
        parse_status=lambda d: d.get("status", "unknown"),
        parse_result=lambda d: d.get("result", d),
    )
    if not result:
        return []
    rows = result.get("residues", result) if isinstance(result, dict) else result
    return [(int(r["index"]), float(r["score"])) for r in rows]


def _parse_residue_scores(stdout: str) -> list[tuple[int, float]]:
    """
    Parse a two-column-ish per-residue table into [(residue_1based, score)].
    Accepts lines beginning with an integer residue index followed by a float
    score somewhere on the line (covers AGGRESCAN3D / CamSol tabular output).
    """
    out: list[tuple[int, float]] = []
    for line in stdout.splitlines():
        cols = line.split()
        if len(cols) < 2:
            continue
        try:
            idx = int(cols[0])
            score = float(cols[-1])
        except ValueError:
            continue
        out.append((idx, score))
    return out


def _run_tool(
    seq: str, name: str, cfg: dict, structure_path: str | None = None,
) -> list[tuple[int, float]]:
    backend = backends.resolve_backend(cfg)
    if backend == "api":
        return _residue_scores_api(seq, cfg)
    return _residue_scores_local(seq, name, cfg, structure_path=structure_path)


def screen_candidate(
    candidate: FusionCandidate,
    name: str,
    cfg: dict,
    window: int = DEFAULT_JUNCTION_WINDOW,
    structure_path: str | None = None,
) -> AggregationResult:
    scores = _run_tool(candidate.sequence, name, cfg, structure_path=structure_path)
    cutoff = float(cfg.get("hotspot_cutoff", 0.0))
    hotspots = [i for i, s in scores if s > cutoff]

    junction_positions = [j.position for j in candidate.junctions]
    near = [
        s for i, s in scores
        if s > cutoff and any(abs(i - jp) <= window for jp in junction_positions)
    ]
    return AggregationResult(
        tool=name,
        junction_hotspot=bool(near),
        max_junction_score=max(near) if near else 0.0,
        n_hotspots=len(hotspots),
    )


def screen_aggregation(
    candidates: list[FusionCandidate],
    fusion_cfg: dict,
    structure_paths: dict[str, str] | None = None,
) -> dict[str, AggregationResult]:
    """Run Gate 2 with the first enabled aggregation tool. {name: result}.

    `structure_paths` maps candidate.name → PDB/CIF path (from the AF3 fold of
    survivors); AGGRESCAN3D uses it, CamSol ignores it."""
    tools = fusion_cfg.get("tools", {})
    structure_paths = structure_paths or {}
    picked = None
    for name in ("aggrescan3d", "camsol"):
        cfg = tools.get(name, {})
        if backends.resolve_backend(cfg) != "disabled":
            picked = (name, cfg)
            break
    if picked is None:
        raise backends.ToolDisabled(
            "No aggregation tool enabled (fusion.tools.aggrescan3d / camsol)."
        )
    name, cfg = picked
    window = int(fusion_cfg.get("junction_window", DEFAULT_JUNCTION_WINDOW))
    print(f"[Gate2] Aggregation tool: {name} ({backends.resolve_backend(cfg)})",
          file=sys.stderr)

    results: dict[str, AggregationResult] = {}
    for cand in candidates:
        res = screen_candidate(
            cand, name, cfg, window=window,
            structure_path=structure_paths.get(cand.name),
        )
        results[cand.name] = res
        print(f"  {cand.name}: junction_hotspot={res.junction_hotspot} "
              f"(max {res.max_junction_score:.2f}, {res.n_hotspots} total)",
              file=sys.stderr)
    return results
