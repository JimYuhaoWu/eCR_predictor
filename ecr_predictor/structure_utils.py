"""
Shared Bio.PDB helpers for AF3 CIF structures.

Used by both foldx.py (function/affinity) and fusion/aggregation.py (Gate 2), so
the CIF-parse / chain-residue / PDB-write boilerplate lives in one place. Bio.PDB
is imported lazily inside each function so importing this module (or its callers)
stays cheap and doesn't require biopython for the sequence-only paths.
"""
from __future__ import annotations

from pathlib import Path


def load_structure(cif_path: str | Path):
    """Parse an mmCIF file into a Bio.PDB Structure."""
    from Bio.PDB import MMCIFParser
    return MMCIFParser(QUIET=True).get_structure(Path(cif_path).stem, str(cif_path))


def first_model(cif_path: str | Path):
    """Return model 0 of an mmCIF structure."""
    return load_structure(cif_path)[0]


def ordered_residues(model, chain_id: str = "A", standard_only: bool = True) -> list:
    """Residues of a chain in sequence order.

    Raises KeyError if the chain is absent. With standard_only, keeps only the
    20 standard amino acids (drops HETATM/ligand/DNA residues)."""
    from Bio.PDB.Polypeptide import is_aa
    if chain_id not in model:
        raise KeyError(f"chain {chain_id} not in structure")
    chain = model[chain_id]
    if standard_only:
        return [r for r in chain.get_residues() if is_aa(r, standard=True)]
    return list(chain.get_residues())


def save_pdb(entity, pdb_path: str | Path) -> None:
    """Write any Bio.PDB entity (Structure / Model / Chain) to a PDB file."""
    from Bio.PDB import PDBIO
    io = PDBIO()
    io.set_structure(entity)
    io.save(str(pdb_path))


def cif_to_pdb(cif_path: str | Path, pdb_path: str | Path) -> None:
    """Convert a whole mmCIF structure to PDB format."""
    save_pdb(load_structure(cif_path), pdb_path)
