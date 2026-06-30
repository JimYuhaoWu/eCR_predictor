"""
Gate 1 — junction immunogenicity (MHC-I neoepitope screen).

For intracellular expression, the dominant immune liability is presentation of
junction-spanning neoepitopes on MHC class I → CD8 killing of transduced cells.
This gate, per candidate:

  1. enumerates junction-spanning peptides (8–11mers)
  2. subtracts self peptides (exact human-proteome matches)
  3. predicts MHC-I presentation across an HLA-I panel
  4. flags binders by %rank and reports epitope DENSITY (flagged / tested)

MHC-I tool backends (config: fusion.tools.<tool>.backend = local|api|disabled):
  - netmhcpan   NetMHCpan -p peptide mode (binding affinity / %rank) — default
  - netctlpan   NetCTLpan (proteasomal cleavage + TAP + MHC-I) — DEPRECATED

NetMHCpan is the default: NetCTLpan 1.1 was discontinued by DTU (superseded by
NetMHCpan), so netmhcpan wins if both are enabled. NetCTLpan additionally models
whether the peptide is actually generated and transported (not just whether it
binds), so the legacy path is kept for pipelines that still have it installed.
Only one MHC-I tool needs to be enabled.

NOTE: the stdout parsers below target NetMHCpan 4.2 / NetCTLpan 1.1 column
layouts. Validate against your installed version — DTU tools occasionally shift
columns between releases.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass

from ecr_predictor.fusion import backends
from ecr_predictor.fusion.assemble import FusionCandidate
from ecr_predictor.fusion.junction import junction_peptides, partition_self

DEFAULT_LENGTHS = (8, 9, 10, 11)
DEFAULT_RANK_THRESHOLD = 2.0   # %rank; <= this is a flagged binder


@dataclass
class ImmunoResult:
    n_tested: int          # unique non-self junction peptides tested
    n_flagged: int         # peptides with %rank <= threshold on any allele
    min_rank: float        # best (lowest) %rank seen (inf if none)
    density: float         # n_flagged / n_tested (0 if none tested)
    flagged_peptides: list[str]


def _pick_mhc1_tool(tools: dict) -> tuple[str, dict] | None:
    """Return (tool_name, tool_cfg) for the enabled MHC-I tool, netmhcpan first
    (NetCTLpan is discontinued; used only when explicitly enabled as a legacy override)."""
    for name in ("netmhcpan", "netctlpan"):
        cfg = tools.get(name, {})
        if backends.resolve_backend(cfg) != "disabled":
            return name, cfg
    return None


# ---------------------------------------------------------------------------
# Tool-specific runners — return {peptide: best_rank_over_alleles}
# ---------------------------------------------------------------------------

def _run_netmhcpan_local(peptides: list[str], cfg: dict) -> dict[str, float]:
    alleles = ",".join(cfg.get("alleles", []))
    lengths = ",".join(str(x) for x in cfg.get("peptide_lengths", DEFAULT_LENGTHS))
    local = cfg.get("local", {})
    command = local.get("command", "netMHCpan")
    extra = local.get("extra_args", "").split()

    pep_text = "\n".join(peptides) + "\n"
    args = ["-p", "{peptides}", "-a", alleles, "-l", lengths, "-BA", *extra]
    stdout = backends.run_cli(
        command, args, input_files={"peptides": pep_text},
        timeout=int(local.get("timeout", 600)),
    )
    return _parse_netmhc_stdout(stdout, rank_col_hint="Rank")


def _run_netctlpan_local(peptides: list[str], cfg: dict) -> dict[str, float]:
    alleles = ",".join(cfg.get("alleles", []))
    local = cfg.get("local", {})
    command = local.get("command", "netCTLpan")
    extra = local.get("extra_args", "").split()

    # NetCTLpan reads a FASTA; we wrap each peptide as its own record.
    fasta = "".join(f">{i}\n{p}\n" for i, p in enumerate(peptides))
    args = ["-a", alleles, "-f", "{fasta}", *extra]
    stdout = backends.run_cli(
        command, args, input_files={"fasta": fasta},
        timeout=int(local.get("timeout", 600)),
    )
    return _parse_netmhc_stdout(stdout, rank_col_hint="Rank")


def _run_mhc1_api(peptides: list[str], cfg: dict) -> dict[str, float]:
    api = cfg.get("api", {})
    payload = {
        "peptides": peptides,
        "alleles": cfg.get("alleles", []),
        "lengths": cfg.get("peptide_lengths", list(DEFAULT_LENGTHS)),
    }
    result = backends.submit_and_poll(
        api,
        payload,
        parse_job_id=lambda d: d.get("id") or d.get("job_id"),
        parse_status=lambda d: d.get("status", "unknown"),
        parse_result=lambda d: d.get("result", d),
    )
    if not result:
        return {}
    # Expected result shape: {"peptide": rank, ...} or list of {peptide, rank}.
    if isinstance(result, dict) and all(isinstance(v, (int, float)) for v in result.values()):
        return {k: float(v) for k, v in result.items()}
    ranks: dict[str, float] = {}
    for row in result if isinstance(result, list) else result.get("predictions", []):
        pep, rank = row.get("peptide"), row.get("rank")
        if pep is not None and rank is not None:
            ranks[pep] = min(rank, ranks.get(pep, float("inf")))
    return ranks


def _parse_netmhc_stdout(stdout: str, rank_col_hint: str = "Rank") -> dict[str, float]:
    """
    Parse DTU NetMHCpan/NetCTLpan stdout into {peptide: best_%rank}.

    The tools print fixed-width tables with a header row containing the peptide
    column and a rank column. We locate columns by header name, then read the
    data rows (lines between the '---' rules). Best (min) rank per peptide is
    kept across alleles.
    """
    ranks: dict[str, float] = {}
    pep_idx = rank_idx = None
    for line in stdout.splitlines():
        cols = line.split()
        if not cols:
            continue
        if pep_idx is None and "Peptide" in cols:
            pep_idx = cols.index("Peptide")
            # Prefer an explicit '%Rank'/'Rank' column; fall back to last col.
            rank_idx = next(
                (i for i, c in enumerate(cols) if rank_col_hint in c or c.startswith("%Rank")),
                len(cols) - 1,
            )
            continue
        if pep_idx is None or line.startswith("#") or line.startswith("-"):
            continue
        if len(cols) <= max(pep_idx, rank_idx):
            continue
        pep = cols[pep_idx]
        try:
            rank = float(cols[rank_idx])
        except ValueError:
            continue
        ranks[pep] = min(rank, ranks.get(pep, float("inf")))
    return ranks


def _run_mhc1(peptides: list[str], name: str, cfg: dict) -> dict[str, float]:
    backend = backends.resolve_backend(cfg)
    if backend == "api":
        return _run_mhc1_api(peptides, cfg)
    if name == "netctlpan":
        return _run_netctlpan_local(peptides, cfg)
    return _run_netmhcpan_local(peptides, cfg)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def screen_candidate(
    candidate: FusionCandidate,
    tool_name: str,
    tool_cfg: dict,
    proteome_blob: str | None,
    lengths: tuple[int, ...] = DEFAULT_LENGTHS,
    rank_threshold: float = DEFAULT_RANK_THRESHOLD,
) -> ImmunoResult:
    """Screen one fusion candidate's junction for MHC-I neoepitopes."""
    peps = junction_peptides(candidate, lengths=lengths)
    if proteome_blob:
        peps, _self = partition_self(peps, proteome_blob)

    unique = sorted({p[0] for p in peps})
    if not unique:
        return ImmunoResult(0, 0, float("inf"), 0.0, [])

    ranks = _run_mhc1(unique, tool_name, tool_cfg)
    flagged = [p for p in unique if ranks.get(p, float("inf")) <= rank_threshold]
    min_rank = min(ranks.values()) if ranks else float("inf")
    density = len(flagged) / len(unique)
    return ImmunoResult(
        n_tested=len(unique),
        n_flagged=len(flagged),
        min_rank=min_rank,
        density=density,
        flagged_peptides=flagged,
    )


def screen_immunogenicity(
    candidates: list[FusionCandidate],
    fusion_cfg: dict,
    proteome_blob: str | None = None,
) -> dict[str, ImmunoResult]:
    """
    Run Gate 1 across all candidates. Returns {candidate.name: ImmunoResult}.
    Raises ToolDisabled if no MHC-I tool is enabled.
    """
    tools = fusion_cfg.get("tools", {})
    picked = _pick_mhc1_tool(tools)
    if picked is None:
        raise backends.ToolDisabled(
            "No MHC-I tool enabled. Set fusion.tools.netctlpan.backend or "
            "fusion.tools.netmhcpan.backend to 'local' or 'api'."
        )
    tool_name, tool_cfg = picked
    rank_threshold = float(tool_cfg.get("rank_threshold", DEFAULT_RANK_THRESHOLD))
    lengths = tuple(tool_cfg.get("peptide_lengths", DEFAULT_LENGTHS))

    print(f"[Gate1] MHC-I tool: {tool_name} ({backends.resolve_backend(tool_cfg)})",
          file=sys.stderr)

    results: dict[str, ImmunoResult] = {}
    for cand in candidates:
        res = screen_candidate(
            cand, tool_name, tool_cfg, proteome_blob,
            lengths=lengths, rank_threshold=rank_threshold,
        )
        results[cand.name] = res
        print(f"  {cand.name}: {res.n_flagged}/{res.n_tested} flagged "
              f"(density {res.density:.2f}, min %rank "
              f"{res.min_rank if res.min_rank != float('inf') else 'NA'})",
              file=sys.stderr)
    return results
