"""
FoldX binding affinity estimation stub.

FoldX operates on a repaired PDB/CIF structure (output from AF3 or experimental).
The standard workflow is:
  1. RepairPDB — fix missing atoms/steric clashes in the input structure
  2. AnalyseComplex — compute interaction energy (ΔΔG) between chains

Prerequisites (when implementing):
  - FoldX binary installed: https://foldxsuite.crg.eu/
  - A valid FoldX licence (free for academic use)
  - Input: PDB/CIF file of the DBD–DNA complex (from af3.py output)

Output contract:
  - Returns a dict {gene_name: ddg_kcal_mol} — interaction ΔΔG in kcal/mol.
  - Lower (more negative) values indicate stronger predicted binding.

See: https://foldxsuite.crg.eu/command/AnalyseComplex
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


def _check_foldx() -> Path:
    """
    Return the path to the FoldX binary, or raise if not found.

    FoldX is not on PATH by default — check common install locations
    or let the user set FOLDX_PATH in their environment.
    """
    import os
    foldx_env = os.environ.get("FOLDX_PATH")
    if foldx_env and Path(foldx_env).exists():
        return Path(foldx_env)

    # Common install conventions
    for candidate in [
        Path.home() / "foldx" / "foldx",
        Path("/usr/local/bin/foldx"),
        Path("/opt/foldx/foldx"),
    ]:
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        "FoldX binary not found. Set the FOLDX_PATH environment variable to its path.\n"
        "Download: https://foldxsuite.crg.eu/"
    )


def _repair_and_analyse(
    foldx_bin: Path,
    structure_path: Path,
    protein_chain: str = "A",
    dna_chain: str = "B",
) -> float | None:
    """
    Run RepairPDB then AnalyseComplex; return the interaction ΔΔG in kcal/mol.

    Returns None if FoldX fails to produce a result.

    TODO: parse the correct output file. FoldX writes:
      - Indiv_energies_<name>_Repair_AC.fxout (per-residue energies)
      - Interaction_<name>_Repair_AC.fxout   (complex interaction energy ← this one)
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        # Copy structure into working dir (FoldX expects input in its run directory)
        import shutil
        local_pdb = tmpdir / structure_path.name
        shutil.copy(structure_path, local_pdb)

        # Step 1: RepairPDB
        result = subprocess.run(
            [
                str(foldx_bin),
                "--command=RepairPDB",
                f"--pdb={local_pdb.name}",
                f"--output-dir={tmpdir}",
            ],
            cwd=tmpdir,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print(f"WARNING: FoldX RepairPDB failed:\n{result.stderr}", file=sys.stderr)
            return None

        repaired_name = local_pdb.stem + "_Repair.pdb"
        repaired_path = tmpdir / repaired_name
        if not repaired_path.exists():
            print(f"WARNING: repaired PDB not found at {repaired_path}", file=sys.stderr)
            return None

        # Step 2: AnalyseComplex
        result = subprocess.run(
            [
                str(foldx_bin),
                "--command=AnalyseComplex",
                f"--pdb={repaired_name}",
                f"--analyseComplexChains={protein_chain},{dna_chain}",
                f"--output-dir={tmpdir}",
            ],
            cwd=tmpdir,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print(f"WARNING: FoldX AnalyseComplex failed:\n{result.stderr}", file=sys.stderr)
            return None

        # TODO: parse Interaction_*_AC.fxout to extract ΔΔG
        # The interaction energy is in column "Interaction Energy" of that file.
        raise NotImplementedError(
            "FoldX output parsing not yet implemented. "
            "See ecr_predictor/foldx.py::_repair_and_analyse for the TODO."
        )


def run_foldx_affinity(
    hits: pd.DataFrame,
    protein_chain: str = "A",
    dna_chain: str = "B",
) -> pd.DataFrame:
    """
    Estimate binding affinity with FoldX for hits that have an AF3 structure.

    Parameters
    ----------
    hits : must have 'gene_name' and 'af3_cif_path' columns (from af3.py)
    protein_chain : chain ID for the DBD in the predicted structure
    dna_chain : chain ID for the DNA in the predicted structure

    Returns
    -------
    hits with new column 'foldx_ddg_kcal_mol' (float or NaN).
    Lower values = stronger predicted binding.

    Raises
    ------
    NotImplementedError — until FoldX parsing is implemented.
    FileNotFoundError — if FoldX binary is not found.
    """
    foldx_bin = _check_foldx()

    hits = hits.copy()
    hits["foldx_ddg_kcal_mol"] = pd.NA

    if "af3_cif_path" not in hits.columns:
        print(
            "WARNING: 'af3_cif_path' column not found. Run AF3 prediction first.",
            file=sys.stderr,
        )
        return hits

    for idx, row in hits.iterrows():
        gene = row["gene_name"]
        cif_path = row.get("af3_cif_path")
        if pd.isna(cif_path):
            continue
        cif_path = Path(cif_path)
        if not cif_path.exists():
            print(f"  SKIP {gene}: structure file not found at {cif_path}", file=sys.stderr)
            continue

        print(f"  Running FoldX on {gene}...", file=sys.stderr)
        ddg = _repair_and_analyse(foldx_bin, cif_path, protein_chain, dna_chain)
        hits.at[idx, "foldx_ddg_kcal_mol"] = ddg if ddg is not None else pd.NA

    return hits
