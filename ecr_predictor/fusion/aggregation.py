"""
Gate 2 — aggregation / solubility liability at the junction.

Intracellular aggregation drives inclusion-body formation and reduced expression.
A candidate is flagged when an aggregation hotspot (per-residue score above the
cutoff) falls within `junction_window` residues of a domain boundary — the
junction is where novel, potentially aggregation-prone sequence appears.

Backends (config: fusion.tools.aggrescan3d.backend = local | api | disabled):

  - local   built-in AGGRESCAN-style scorer (no external binary). Projects the
            published AGGRESCAN a3v aggregation-propensity scale onto the chain,
            smooths it with the length-adaptive sliding window, and — when the
            AF3 structure of a survivor is available — weights each residue by its
            relative solvent accessibility (freesasa), reproducing AGGRESCAN3D's
            core idea (aggregation propensity x exposure). Falls back to the
            sequence-only profile when no structure or freesasa is available.
  - api     POST sequence to a remote service (submit/poll); expects a per-residue
            [{index, score}, ...] result.

CamSol is web-only (no CLI) — usable here only via the `api` backend.

a3v scale: Conchillo-Solé et al., BMC Bioinformatics 2007 (AGGRESCAN). The window
is length-adaptive per that paper (5 / 7 / 9 / 11 residues); HST (the frequency-
weighted mean of a3v, ~ -0.02) is the reference hot-spot threshold for the
sequence profile — tune `hotspot_cutoff` in config to match your scale.
"""
from __future__ import annotations

import math
import os
import sys
import tempfile
from dataclasses import dataclass

from ecr_predictor.fusion import backends
from ecr_predictor.fusion.assemble import FusionCandidate

DEFAULT_JUNCTION_WINDOW = 5

# AGGRESCAN a3v aggregation-propensity scale (per residue), 3-letter keyed.
# Conchillo-Solé et al., BMC Bioinformatics 2007.
A3V = {
    "ALA": -0.036, "ARG": -1.240, "ASN": -1.302, "ASP": -1.836, "CYS": 0.604,
    "GLN": -1.231, "GLU": -1.412, "GLY": -0.535, "HIS": -0.733, "ILE": 1.822,
    "LEU": 1.380, "LYS": -0.931, "MET": 0.910, "PHE": 1.754, "PRO": -0.334,
    "SER": -0.294, "THR": -0.159, "TRP": 1.037, "TYR": 1.159, "VAL": 1.594,
}
_ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}
# Frequency-weighted mean of a3v (AGGRESCAN HST) — reference hot-spot threshold.
A3V_HST = -0.02


@dataclass
class AggregationResult:
    tool: str
    junction_hotspot: bool          # any hotspot within the junction window
    max_junction_score: float       # worst per-residue score near a junction
    n_hotspots: int                 # total hotspots in the candidate


# ---------------------------------------------------------------------------
# AGGRESCAN a3v profile helpers
# ---------------------------------------------------------------------------

def _adaptive_window(n: int) -> int:
    """AGGRESCAN length-adaptive sliding-window size (Conchillo-Solé 2007)."""
    if n <= 75:
        return 5
    if n <= 175:
        return 7
    if n <= 300:
        return 9
    return 11


def _smooth(values: list[float], window: int) -> list[float]:
    """Centered moving average; the window shrinks symmetrically at the termini."""
    n = len(values)
    half = window // 2
    out: list[float] = []
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        chunk = values[lo:hi]
        out.append(sum(chunk) / len(chunk))
    return out


def _residue_scores_sequence(seq: str) -> list[tuple[int, float]]:
    """Sequence-only AGGRESCAN a3v profile, [(residue_1based, score), ...]."""
    raw = [A3V.get(_ONE_TO_THREE.get(aa.upper(), ""), 0.0) for aa in seq]
    smoothed = _smooth(raw, _adaptive_window(len(raw)))
    return [(i + 1, s) for i, s in enumerate(smoothed)]


def _residue_scores_structure(
    cif_path: str, protein_chain: str = "A",
) -> list[tuple[int, float]]:
    """Structure-weighted AGGRESCAN profile from the AF3 model.

    Smoothed a3v per residue x relative solvent accessibility (freesasa) for the
    protein chain — high score = aggregation-prone AND solvent-exposed. Raises
    ImportError if freesasa is unavailable; the caller falls back to sequence mode.
    """
    import freesasa  # raises ImportError -> sequence fallback
    from Bio.PDB import MMCIFParser, PDBIO
    from Bio.PDB.Polypeptide import is_aa

    structure = MMCIFParser(QUIET=True).get_structure("m", cif_path)
    model = structure[0]
    if protein_chain not in model:
        raise KeyError(f"chain {protein_chain} not in {cif_path}")
    chain = model[protein_chain]
    residues = [r for r in chain.get_residues() if is_aa(r, standard=True)]
    if not residues:
        raise ValueError(f"no standard residues in chain {protein_chain}")

    # Isolate the protein chain to a temp PDB so freesasa scores the monomer
    # surface (not buried by DNA/ligand chains), as AGGRESCAN3D does.
    fd, pdb_path = tempfile.mkstemp(suffix=".pdb")
    os.close(fd)
    try:
        io = PDBIO()
        io.set_structure(chain)
        io.save(pdb_path)
        result = freesasa.calc(freesasa.Structure(pdb_path))
    finally:
        try:
            os.remove(pdb_path)
        except OSError:
            pass

    chain_areas = {
        str(k).strip(): v
        for k, v in result.residueAreas().get(protein_chain, {}).items()
    }

    raw = [A3V.get(r.resname.upper(), 0.0) for r in residues]
    smoothed = _smooth(raw, _adaptive_window(len(raw)))

    out: list[tuple[int, float]] = []
    for idx, (res, prop) in enumerate(zip(residues, smoothed), start=1):
        ra = chain_areas.get(str(res.id[1]).strip())
        rel = getattr(ra, "relativeTotal", None) if ra is not None else None
        if rel is None or not math.isfinite(rel):
            rel = 1.0  # missing reference area -> assume exposed (conservative)
        out.append((idx, prop * rel))
    return out


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


def _run_tool(
    seq: str, name: str, cfg: dict, structure_path: str | None = None,
) -> tuple[list[tuple[int, float]], str]:
    """Return (per-residue scores, mode-label) for the selected backend."""
    backend = backends.resolve_backend(cfg)
    if backend == "api":
        return _residue_scores_api(seq, cfg), f"{name}:api"
    if name == "camsol":
        raise backends.ToolNotAvailableError(
            "camsol has no local CLI — set fusion.tools.camsol.backend: api, "
            "or use aggrescan3d (built-in)."
        )
    if structure_path:
        try:
            return _residue_scores_structure(structure_path), "aggrescan-a3v+freesasa"
        except ImportError:
            print("  WARNING: freesasa not installed — sequence-only AGGRESCAN "
                  "scoring (conda install -c conda-forge freesasa).", file=sys.stderr)
        except (FileNotFoundError, KeyError, ValueError) as e:
            print(f"  WARNING: structure scoring failed ({e}) — sequence-only "
                  "AGGRESCAN scoring.", file=sys.stderr)
    return _residue_scores_sequence(seq), "aggrescan-a3v(seq)"


def screen_candidate(
    candidate: FusionCandidate,
    name: str,
    cfg: dict,
    window: int = DEFAULT_JUNCTION_WINDOW,
    structure_path: str | None = None,
) -> AggregationResult:
    scores, mode = _run_tool(
        candidate.sequence, name, cfg, structure_path=structure_path,
    )
    cutoff = float(cfg.get("hotspot_cutoff", 0.0))
    hotspots = [i for i, s in scores if s > cutoff]

    junction_positions = [j.position for j in candidate.junctions]
    near = [
        s for i, s in scores
        if s > cutoff and any(abs(i - jp) <= window for jp in junction_positions)
    ]
    return AggregationResult(
        tool=mode,
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

    `structure_paths` maps candidate.name → AF3 CIF path (from the fold of
    survivors); the built-in scorer SASA-weights with it when present."""
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
    print(f"[Gate2] Aggregation: {name} ({backends.resolve_backend(cfg)})",
          file=sys.stderr)

    results: dict[str, AggregationResult] = {}
    for cand in candidates:
        res = screen_candidate(
            cand, name, cfg, window=window,
            structure_path=structure_paths.get(cand.name),
        )
        results[cand.name] = res
        print(f"  {cand.name}: junction_hotspot={res.junction_hotspot} "
              f"(max {res.max_junction_score:.2f}, {res.n_hotspots} total, {res.tool})",
              file=sys.stderr)
    return results
