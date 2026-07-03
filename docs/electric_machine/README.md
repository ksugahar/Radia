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

The solver behind `planar_vim_motor.ipynb` is the promoted 2D layer in `radia.vim`
(`PlanarDemagBody` / `hdiv_demag_solve` on a 2D mesh / `maxwell_torque_circle`), golden-locked in
`validation_test/feec/test_hdiv_vim_2d_solve.py`.

The executable validation corpus is
`validation_test/electric_machine/`; this directory is the rendered docs layer.
