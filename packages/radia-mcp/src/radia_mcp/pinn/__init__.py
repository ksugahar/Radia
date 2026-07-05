"""radia_mcp.pinn — Physics-Informed Neural Networks and Gaussian Processes
for Maxwell's equations and computational EM.

Two parallel families of physics-informed AI methods:

1. PINN (Physics-Informed Neural Network) — Raissi-Karniadakis 2019 et seq.
   Neural network trained to satisfy a PDE residual + data loss.
   Mesh-free, scales to high dimensions, slow to train but fast to query.

2. PI-GP (Physics-Informed Gaussian Process) — Raissi-Karniadakis 2017
   et seq.  Gaussian process prior with kernel transformed by the
   PDE operator: k_ff = L_x L_x' k_uu.  Bayesian, provides uncertainty
   quantification, scales to ~10⁴ training points.

Both are alternatives to classical FEM for specific use cases:
- Inverse problems (parameter identification from observations)
- Multi-fidelity fusion (low + high fidelity data)
- Real-time prediction after training

Distilled from:
  - Raissi 2017 (PI-GP foundation), 2019 (PINN foundation)
  - Pförtner et al 2023 (PI-GP generalizes linear PDE solvers)
  - Henderson 2023 (wave equation PI-GP)
  - public-safe curated corpus Neural Network\\
  - public-safe curated corpus
"""
