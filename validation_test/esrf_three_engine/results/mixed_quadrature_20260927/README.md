# ESRF 6 mixed Omega material quadrature

The accepted 5.0.1 wheel was used for two fresh mixed Omega solves on the same
Kelvin mesh. Both used Newton, P2 response, assembly bonus 12, source projection
order 2, surface-flux loading, and four NGSolve threads. Only the material
quadrature bonus changed, from 4 to 8. The recorded solver statistics confirm
the actual bonus; filenames alone are not used as evidence.

Both solves converged in six Newton iterations. The relative vector difference
over all 45 observation volumes is `1.6319481276745185e-7`, below the declared
`1e-4` quadrature-convergence gate. Each observation uses the same eight-point
volume average with a 20 micrometre half-width. `result.json` retains both field
arrays, coordinates, true nonlinear residuals, and input/implementation hashes.

The calculations called `solve_omega` in
`validation_test/c_type_three_engine/run_three_engine.py`, using the case-6
source, B-H table, mesh, and observation helpers from the coil-yoke runner.
The public three-engine runner exposes the same settings as
`--mixed-method newton --fem-order 2 --mixed-bonus 12 --mixed-source-order 2
--mixed-source-load surface_flux --mixed-material-bonus 4` (or `8`).
Use separate outputs for the two material rules; their checkpoints are not
interchangeable.

This establishes material-quadrature stability on this mesh. It does not
establish mesh convergence, absolute field accuracy, or full three-engine
acceptance. Read-only profiler sampling occurred during this campaign, so its
elapsed times are not published as performance evidence.
