"""radia_mcp.data_assimilation — Kalman / Ensemble / 4D-Var for EM.

Distilled from public-safe curated corpus
  - Recent Progress of Data Assimilation Methods in Meteorology
  - トロッコ問題で理解するカルマンフィルタ

Sensor fusion + state estimation for EM applications: measure B at
a few points, infer field everywhere; track moving conductor under
noise; combine FEM forward + measurement update.

Cross-links:
  - `radia_mcp.pinn` — physics-aware filter via neural surrogate
  - `radia_mcp.mor` — reduced-order model in the filter's state space
"""
