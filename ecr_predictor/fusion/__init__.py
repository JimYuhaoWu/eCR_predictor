"""
Fusion-protein developability screening (Step 3).

Given selected DBDs (Step 1/2) plus user-supplied effector domains (EDs) and a
linker library, assemble fusion candidates and screen them through a series of
in-silico developability gates BEFORE wet-lab synthesis:

  Gate 0  Function retention   reuse af3.py + foldx.py (optional, structural)
  Gate 1  Junction immunogenicity   MHC-I junction neoepitope screen
  Gate 2  Aggregation/solubility    AGGRESCAN3D + CamSol
  Gate 3  Intracellular stability   N-end rule / degron / ubiquitination
  Score   Composite risk + Pareto ranking of candidates

Target modality is INTRACELLULAR expression (viral vector / mRNA), so the
dominant immune axis is MHC-I presentation of junction-spanning neoepitopes to
CD8 T cells (loss of transduced cells), not antibody/serum-protease routes.

Each external tool is invoked through a local-CLI-or-API backend selected in
config.yaml (fusion.tools.<tool>.backend), mirroring the AF3 stage.
"""
