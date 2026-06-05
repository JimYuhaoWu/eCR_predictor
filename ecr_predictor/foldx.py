"""
FoldX binding affinity estimation.

Workflow per hit:
  1. Convert AF3 CIF → PDB (FoldX 5.1 does not read CIF)
  2. RepairPDB  — fix missing atoms / steric clashes
  3. AnalyseComplex — compute interaction ΔΔG between protein and DNA chains

Prerequisites:
  - FoldX binary: https://foldxsuite.crg.eu/ (free academic licence)
  - Set FOLDX_PATH env var to the binary, e.g.:
      export FOLDX_PATH=/mnt3/wuyuhao/foldx/foldx
  - biopython (already a project dependency) — used for CIF → PDB conversion

Note: FoldX cannot resolve '~' in paths. All paths passed to it are
fully resolved with Path.resolve() before use.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Binary discovery
# ---------------------------------------------------------------------------

def _check_foldx() -> Path:
    """Return the FoldX binary path, or raise FileNotFoundError."""
    foldx_env = os.environ.get("FOLDX_PATH")
    if foldx_env and Path(foldx_env).exists():
        return Path(foldx_env).resolve()

    for candidate in [
        Path.home() / "foldx" / "foldx",
        Path("/usr/local/bin/foldx"),
        Path("/opt/foldx/foldx"),
    ]:
        if candidate.exists():
            return candidate.resolve()

    raise FileNotFoundError(
        "FoldX binary not found. Set the FOLDX_PATH environment variable:\n"
        "  export FOLDX_PATH=/path/to/foldx\n"
        "Download: https://foldxsuite.crg.eu/"
    )


# ---------------------------------------------------------------------------
# CIF → PDB conversion
# ---------------------------------------------------------------------------

def _cif_to_pdb(cif_path: Path, pdb_path: Path) -> None:
    """Convert an mmCIF file to PDB format using BioPython."""
    from Bio.PDB import MMCIFParser, PDBIO
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure(cif_path.stem, str(cif_path))
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(pdb_path))


# ---------------------------------------------------------------------------
# FoldX steps
# ---------------------------------------------------------------------------

def _run_repair(foldx_bin: Path, pdb_name: str, pdb_dir: Path, out_dir: Path) -> bool:
    """Run RepairPDB. Returns True on success."""
    result = subprocess.run(
        [
            str(foldx_bin),
            "--command=RepairPDB",
            f"--pdb={pdb_name}",
            f"--pdb-dir={pdb_dir}",
            f"--output-dir={out_dir}",
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or "Your file run OK" not in result.stdout:
        print(f"  WARNING: RepairPDB failed:\n{result.stdout[-500:]}", file=sys.stderr)
        return False
    return True


def _run_analyse(
    foldx_bin: Path,
    repaired_name: str,
    pdb_dir: Path,
    out_dir: Path,
    protein_chain: str,
    dna_chain: str,
) -> bool:
    """Run AnalyseComplex. Returns True on success."""
    result = subprocess.run(
        [
            str(foldx_bin),
            "--command=AnalyseComplex",
            f"--pdb={repaired_name}",
            f"--pdb-dir={pdb_dir}",
            f"--analyseComplexChains={protein_chain},{dna_chain}",
            f"--output-dir={out_dir}",
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or "Your file run OK" not in result.stdout:
        print(f"  WARNING: AnalyseComplex failed:\n{result.stdout[-500:]}", file=sys.stderr)
        return False
    return True


def _parse_interaction_energy(fxout_path: Path) -> float | None:
    """
    Parse the Interaction Energy from an Interaction_*_AC.fxout file.

    Format (tab-separated, one header row, one data row):
      Pdb  Group1  Group2  IntraclashesGroup1  IntraclashesGroup2
      Interaction Energy  Backbone Hbond  ...

    Example data row:
      /path/to/fli1_model_Repair.pdb  A  B  24.89  9.02  -5.07  ...
    """
    if not fxout_path.exists():
        print(f"  WARNING: FoldX output not found: {fxout_path}", file=sys.stderr)
        return None

    with fxout_path.open(encoding="utf-8") as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("FoldX") and
                 not l.startswith("by ") and not l.startswith("Jesper") and
                 not l.startswith("Luis") and not l.startswith("---")]

    # Find header and data lines
    header_line = next((l for l in lines if "Interaction Energy" in l), None)
    if header_line is None:
        print(f"  WARNING: could not find header in {fxout_path}", file=sys.stderr)
        return None

    headers = [h.strip() for h in header_line.split("\t")]
    data_line = lines[lines.index(header_line) + 1].strip()
    values = data_line.split("\t")

    if len(values) != len(headers):
        print(f"  WARNING: column count mismatch in {fxout_path}", file=sys.stderr)
        return None

    row = dict(zip(headers, values))
    try:
        return float(row["Interaction Energy"])
    except (KeyError, ValueError) as e:
        print(f"  WARNING: could not parse Interaction Energy: {e}", file=sys.stderr)
        return None


# ---------------------------------------------------------------------------
# Main per-hit function
# ---------------------------------------------------------------------------

def _repair_and_analyse(
    foldx_bin: Path,
    cif_path: Path,
    work_dir: Path,
    protein_chain: str,
    dna_chain: str,
) -> float | None:
    """
    Convert CIF → PDB, run RepairPDB, run AnalyseComplex, return ΔΔG.

    work_dir is a persistent per-gene directory (not a tempdir) so outputs
    are inspectable and RepairPDB is not re-run if already done.
    """
    work_dir = work_dir.resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    # CIF → PDB
    pdb_name = cif_path.stem + ".pdb"
    pdb_path = work_dir / pdb_name
    if not pdb_path.exists():
        print(f"    Converting CIF → PDB...", file=sys.stderr)
        _cif_to_pdb(cif_path.resolve(), pdb_path)

    # RepairPDB (skip if already done)
    repaired_name = cif_path.stem + "_Repair.pdb"
    repaired_path = work_dir / repaired_name
    if not repaired_path.exists():
        print(f"    Running RepairPDB (this takes ~4 min)...", file=sys.stderr)
        ok = _run_repair(foldx_bin, pdb_name, work_dir, work_dir)
        if not ok or not repaired_path.exists():
            return None
    else:
        print(f"    RepairPDB already done, skipping.", file=sys.stderr)

    # AnalyseComplex
    fxout_name = f"Interaction_{cif_path.stem}_Repair_AC.fxout"
    fxout_path = work_dir / fxout_name
    if not fxout_path.exists():
        print(f"    Running AnalyseComplex...", file=sys.stderr)
        ok = _run_analyse(foldx_bin, repaired_name, work_dir, work_dir,
                          protein_chain, dna_chain)
        if not ok:
            return None
    else:
        print(f"    AnalyseComplex already done, skipping.", file=sys.stderr)

    return _parse_interaction_energy(fxout_path)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_foldx_affinity(
    hits: pd.DataFrame,
    protein_chain: str = "A",
    dna_chain: str = "B",
    work_base_dir: str | Path = "foldx_work",
) -> pd.DataFrame:
    """
    Estimate binding affinity with FoldX for hits that have an AF3 structure.

    Parameters
    ----------
    hits : must have 'gene_name' and 'af3_cif_path' columns (from af3.py)
    protein_chain : chain ID for the DBD in the predicted structure (default A)
    dna_chain     : chain ID for the DNA in the predicted structure (default B)
    work_base_dir : directory for FoldX intermediate files (default foldx_work/)

    Returns
    -------
    hits with new column 'foldx_ddg_kcal_mol' (float or NA).
    Lower (more negative) = stronger predicted binding.
    """
    foldx_bin = _check_foldx()
    work_base_dir = Path(work_base_dir).resolve()

    hits = hits.copy()
    hits["foldx_ddg_kcal_mol"] = pd.NA

    if "af3_cif_path" not in hits.columns:
        print("WARNING: 'af3_cif_path' column not found. Run AF3 prediction first.",
              file=sys.stderr)
        return hits

    for idx, row in hits.iterrows():
        gene = str(row["gene_name"])
        cif_path = row.get("af3_cif_path")
        if pd.isna(cif_path):
            continue
        cif_path = Path(cif_path)
        if not cif_path.exists():
            print(f"  SKIP {gene}: CIF not found at {cif_path}", file=sys.stderr)
            continue

        print(f"  [FoldX] {gene}...", file=sys.stderr)
        work_dir = work_base_dir / gene.lower()
        ddg = _repair_and_analyse(foldx_bin, cif_path, work_dir,
                                  protein_chain, dna_chain)
        if ddg is not None:
            hits.at[idx, "foldx_ddg_kcal_mol"] = ddg
            print(f"    ΔΔG = {ddg:.3f} kcal/mol", file=sys.stderr)
        else:
            print(f"    WARNING: FoldX failed for {gene}", file=sys.stderr)

    return hits
