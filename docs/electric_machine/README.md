# Electric-Machine Docs

This docs topic now owns the human-facing electric-machine demo material that
used to live under `examples/electric_machine`.

| Artifact | Purpose |
|---|---|
| `cogging_skew_demo.ipynb` | Result-saved cogging / skew notebook with synchronized JSON. |
| `cogging_skew_demo.py` | Notebook-coupled helper for the finite-element torque sweep. |
| `cogging_skew_demo_results.json` | Computed torque, skew, checks, and runtime versions. |
| `planar_vim_motor.ipynb` | Executed planar HDiv-VIM machine showcase (2026-07-03): nonlinear iron deep-saturation vs the analytic fixed point; ellipse reluctance torque three ways with the demag operator built once per rotation sweep; salient-bar motor torque-angle vs an exact-Newton FEM (mean 0.58%); rotating-field conducting cylinder vs the Bessel closed form (0.19%); mini cage induction machine (VIM core + reduced-FEM bars) torque-slip vs an all-in-one FEM (0.57%). |
| `planar_vim_motor_helpers.py` | Notebook-coupled reference/coupling layer (exact-Newton nonlinear FEM reference, Bessel closed form, single-valued polar conjugate potential, cage stagger, frozen-secant all-in-one). |
| `planar_vim_motor_result.json` | Sidecar with the executed outputs, versions, and notebook sha. |
| `em_reference_audit.ipynb` | Executed reference-audit methodology (2026-07-04): the diagnostics that exonerate or convict a FEM cross-validation reference — coil-disk polygon current deficit (−5.4%) + the drive-equivalence probe (uniform 4.9%); the conjugate-potential gradient gate (1e-10) and atan2 branch cuts (72% jumps) vs the single-valued polar construction (closure 4e-15); the finite-Dirichlet dipole image matched by its closed form (0.00%) and the exact n=1 open Robin; the secant-Picard plateau vs the exact Newton (9 iters). Cross-linked to `bug_patterns_lookup(topic="validation")` and MCP `hdiv_vim(topic="reference_audit")`. |
| `em_reference_audit_result.json` | Sidecar for the audit notebook. |

The solver behind `planar_vim_motor.ipynb` is the promoted 2D layer in `radia.vim`
(`PlanarDemagBody` / `Solve` on a 2D mesh / `maxwell_torque_circle`), golden-locked in
`validation_test/feec/test_hdiv_vim_2d_solve.py`.

The executable validation corpus is
`validation_test/electric_machine/`; this directory is the rendered docs layer.
