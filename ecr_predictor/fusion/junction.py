"""
Junction-spanning peptide enumeration and self-proteome tolerance filtering.

Because both fused domains are endogenous, peptides lying entirely within either
domain are self and are immunologically tolerated. Only peptides that STRADDLE a
domain boundary are candidate neoepitopes. This module enumerates exactly those
peptides and (optionally) removes any that also occur verbatim in the human
proteome — those are self too, even if they happen to cross the boundary.
"""
from __future__ import annotations

from pathlib import Path

from ecr_predictor.fusion.assemble import FusionCandidate


def junction_peptides(
    candidate: FusionCandidate,
    lengths: tuple[int, ...] = (8, 9, 10, 11),
) -> list[tuple[str, int, int]]:
    """
    Enumerate all k-mers (k in `lengths`) that overlap at least one residue on
    each side of a domain boundary.

    Returns a list of (peptide, start_1based, junction_index).
    """
    seq = candidate.sequence
    out: list[tuple[str, int, int]] = []
    for k in lengths:
        for start in range(0, len(seq) - k + 1):
            first = start + 1            # 1-based first residue
            last = start + k             # 1-based last residue (inclusive)
            for ji, j in enumerate(candidate.junctions):
                # window straddles the boundary between j.position and j.position+1
                if first <= j.position and last >= j.position + 1:
                    out.append((seq[start:start + k], first, ji))
                    break
    return out


def load_proteome_blob(fasta_path: str | Path) -> str:
    """
    Load a proteome FASTA into a single newline-delimited blob for fast
    substring membership tests. Newlines act as record separators so peptides
    cannot span two proteins.
    """
    seqs: list[str] = []
    cur: list[str] = []
    with Path(fasta_path).open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return "\n".join(seqs)


def partition_self(
    peptides: list[tuple[str, int, int]],
    proteome_blob: str,
) -> tuple[list[tuple[str, int, int]], list[tuple[str, int, int]]]:
    """
    Split peptides into (non_self, self) by exact occurrence in the proteome blob.
    A peptide present anywhere in the proteome is treated as self/tolerated.
    """
    non_self = [p for p in peptides if p[0] not in proteome_blob]
    self_p = [p for p in peptides if p[0] in proteome_blob]
    return non_self, self_p
