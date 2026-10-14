"""radia_mcp.bayesian_opt — Bayesian optimization, Gaussian Process
regression, FMQA, and surrogate models for EM engineering.

Distilled from public-safe curated corpus (57 lab files), including:
  - GP regression theory: ガウス過程回帰の基礎.pdf, Mochihashi GP tutorial
  - ARD kernel for nonlinear EM/seismic response: 2 papers
  - Physics-informed GP: Raissi 2017, Pförtner 2023, Henderson 2023
    (Gaussian process regression generalizes linear PDE solvers)
  - FMQA: Factorization Machines + Quantum Annealing applied to
    materials informatics + metamaterial design
  - CAE懇話会 2024 surrogate model workshop (5 sessions)
  - GP for multi-objective + robotic control
  - LLM-enhanced BO (Large Language Models to Enhance Bayesian Opt)

Cross-links:
  - official `optuna/optuna-mcp` — practical Optuna study/trial
    orchestration when the workflow needs TPE, pruning, or dashboard
    analysis outside radia-mcp
  - `radia_mcp.pinn` — Physics-Informed Neural Networks (sister to
    Physics-Informed Gaussian Processes covered here)
  - `radia_mcp.topology_optimization` — gradient-based alternative
  - `radia_mcp.metamaterial` — application target for FMQA
"""
