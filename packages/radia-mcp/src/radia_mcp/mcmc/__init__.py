"""radia_mcp.mcmc — MCMC + MCTS + Sampled Pattern Matching for EM
engineering optimization.

Knowledge distilled from:
  - Wada-Sugino-Mokishi, "ゼロからできるMCMC"
    (Asakura Shoten, 2019, ISBN 9784254127292) — statistical
    physics textbook on Markov Chain Monte Carlo. Lab copy:
    W:/.../00_教科書/10_統計_ベイズ_MCMC/ゼロからできるMCMC_*.pdf
  - 須山敦志 (Suyama), "ベイズ推論による機械学習入門" (MLP series)
    + "ベイズ深層学習" — Bayesian inference / NN textbooks
  - Shadare-Sadiku-Musa, "Markov Chain Monte Carlo Solution of
    Laplace's Equation in Axisymmetric Homogeneous Domain", Open
    J Modelling and Simulation 7(4):203-216, 2019 — direct MCMC
    field solver.
  - **Sato-Igarashi 2023** (Hokkaido, IEEE T-MAG 59(5)):
    "Multi-Objective Automatic Design of Permanent Magnet Motor
    Using Monte Carlo Tree Search" — MCTS + multi-objective TO
    (NGnet + NSGA-II) for PM motor design. ★ lab focus.
  - **Yin-Sato-Igarashi 2024** (Hokkaido, IEEE T-MAG 60(3)):
    "A Comprehensive Optimal Design of Inductors Using MCTS" —
    MCTS + CMA-ES for inductor design with non-linear materials.
  - Saotome-Coulomb-Saito-Sabonnadiere 1995 (IEEE T-MAG 31(3)):
    "Magnetic Core Shape Design by the Sampled Pattern Matching
    Method" — SPM cosine-similarity inverse problem method.

Cross-links to lab MCP layers:
  - `radia_mcp.optuna` — BBO alternative (TPE / NSGA-II / CMA-ES)
  - `radia_mcp.topology_optimization` — gradient-based (SIMP, MMA)
  - `radia_mcp.motor` — PM motor analysis backend for Hayaho MCTS
  - `radia_mcp.magnetic_materials` — Play/Energy hysteresis ID via MCMC
  - `radia_mcp.pinn` — surrogate model in MCMC posterior evaluation
  - `radia_mcp.mathematica` — symbolic derivation of detailed balance
"""
