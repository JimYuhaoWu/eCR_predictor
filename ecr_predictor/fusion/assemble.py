"""
Assemble fusion-protein candidates from DBDs, linkers, and effector domains.

A fusion is DBD + linker + ED (or ED + linker + DBD, set by `orientation`).
The cartesian product over {DBDs} × {linkers} × {EDs} forms the candidate
library that the developability gates score. Each candidate tracks the residue
position of every domain boundary ("junction") so downstream gates can restrict
analysis to the novel, junction-spanning sequence.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from typing import Iterable


@dataclass
class Junction:
    """A domain boundary. `position` is the 1-based index of the last residue
    on the N-terminal side; the boundary sits between `position` and `position+1`."""
    position: int
    left: str
    right: str


@dataclass
class FusionCandidate:
    name: str
    dbd_name: str
    linker_name: str
    ed_name: str
    dbd_seq: str
    linker_seq: str
    ed_seq: str
    orientation: str = "dbd_n"          # "dbd_n" (DBD at N-terminus) | "ed_n"
    sequence: str = ""
    junctions: list[Junction] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.orientation == "dbd_n":
            parts = [("DBD", self.dbd_seq), ("LINK", self.linker_seq), ("ED", self.ed_seq)]
        elif self.orientation == "ed_n":
            parts = [("ED", self.ed_seq), ("LINK", self.linker_seq), ("DBD", self.dbd_seq)]
        else:
            raise ValueError(f"Unknown orientation {self.orientation!r}. Use dbd_n | ed_n.")

        present = [(label, s) for label, s in parts if s]
        seq = ""
        junctions: list[Junction] = []
        for i, (label, s) in enumerate(present):
            if i > 0:
                junctions.append(Junction(position=len(seq), left=present[i - 1][0], right=label))
            seq += s
        self.sequence = seq
        self.junctions = junctions


def _normalize_linkers(linkers: Iterable) -> list[tuple[str, str]]:
    """Accept config linkers as a list of {name, sequence} dicts (or (name, seq)
    tuples) and return a list of (name, sequence)."""
    out: list[tuple[str, str]] = []
    for item in linkers:
        if isinstance(item, dict):
            out.append((str(item["name"]), str(item["sequence"])))
        else:
            name, seq = item
            out.append((str(name), str(seq)))
    return out


def assemble_fusions(
    dbds: dict[str, str],
    eds: dict[str, str],
    linkers: Iterable,
    orientation: str = "dbd_n",
) -> list[FusionCandidate]:
    """
    Build the full candidate library (DBDs × linkers × EDs).

    Parameters
    ----------
    dbds, eds : {name: amino-acid sequence}
    linkers   : list of {name, sequence} dicts or (name, sequence) tuples
    orientation : "dbd_n" or "ed_n"
    """
    linker_pairs = _normalize_linkers(linkers)
    candidates: list[FusionCandidate] = []
    for (dn, ds), (ln, ls), (en, es) in product(dbds.items(), linker_pairs, eds.items()):
        name = f"{dn}__{ln}__{en}"
        candidates.append(
            FusionCandidate(
                name=name,
                dbd_name=dn, linker_name=ln, ed_name=en,
                dbd_seq=ds, linker_seq=ls, ed_seq=es,
                orientation=orientation,
            )
        )
    return candidates
