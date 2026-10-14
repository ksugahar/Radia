# Kelvin Adaptive and TEAM7 Retirement Note

The former Kelvin adaptive-mesh and TEAM7 scripts were development runners, not
maintained examples. They mixed solver experiments, local path setup, adaptive
mesh sweeps, and unfinished TEAM7 comparisons. The durable outcomes have been
split into maintained locations:

- Kelvin material and pullback math: `src/radia/kelvin_source.py` plus
  `validation_test/kelvin_source/`.
- Cubit/Omega-Reduced Omega high-order Kelvin regression:
  `validation_test/cubit/kelvin_1_4_p_convergence/`.
- Adaptive-method notes and the retained CG smoother demonstration:
  `docs/kelvin/Supplement/CG-smoother.md`,
  `docs/kelvin/Supplement/ErrorEstimator.md`, and
  `docs/kelvin/Supplement/cg_smoother_demo.ipynb`.
- General TEAM7 knowledge remains in MCP via published TEAM benchmark notes,
  not through a Kelvin examples archive.

If TEAM7 Kelvin coupling is needed again, rebuild it as a fresh
`validation_test/` target with public input data, explicit pass criteria, and a
result JSON. Do not resurrect the old archive as source.
