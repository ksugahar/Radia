"""radia_mcp.optuna — MCP server for Optuna black-box optimization
in the context of Radia / NGSolve EM analysis.

Knowledge distilled from:
  - Sano, Akiba, Imamura, Ota, Mizuno, Yanase,
    "Optuna によるブラックボックス最適化" (Ohmsha 2023.2,
    ISBN 9784274230103) -- the official Optuna team textbook.
    Repo of book code: chapter2 / chapter3 / chapter4 (Bayesian
    cookie) / chapter5 (internals).
  - Optuna paper: Akiba, Sano, Yanase, Ohta, Koyama,
    "Optuna: A Next-generation Hyperparameter Optimization
    Framework", KDD 2019.
  - https://optuna.org (project home), https://optuna.readthedocs.io

Tied to lab practice via cross-links to:
  - `radia_mcp.topology_optimization` (SIMP/MMA -- gradient-based)
  - `radia_mcp.motor` (Wakao autoencoder + LS SynRM topology opt)
  - `radia_mcp.accelerator` (coil end-region shape optimization)
  - `radia_mcp.pcb` (compensation topology selection)
  - `radia_mcp.pinn` (BBO of PINN architecture / weights)
"""
