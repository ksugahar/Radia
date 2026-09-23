# BDM2 mesh sensitivity (three levels)

HDiv-only diagnostic, not a three-engine pass, continuum accuracy certificate,
or matched-error performance comparison. Checked medium/fine meshes, BDM2,
linear mu_r1000, 8 threads, Gram eps1e-14, solve tol1e-10 and direct field.
Coarse is reused from ../hibino_20260924_order/order_comparison.json.
Legacy difference fields in each new JSON are self differences and have no
comparative meaning; the JSON note explicitly states this.

Medium: 10860 DOF, 17.2653 s, 48 iterations, residual7.512853472e-11.
Fine: 16608 DOF, 27.2219 s, 48 iterations, residual9.063174227e-11.
Both converged. Core increments using the same fine denominator and campaign
axial symmetry projection are 0.0026692% then 0.2422117% (ratio90.7432).
The second increment is smaller than BDM1's roughly1%, but increments are
not contracting. A tiny first increment on independent meshes must not be
mistaken for convergence. The fourth mesh level is still needed.

Raw archive S:/Radia/validation_artifacts/hdiv_bdm2levels_20260924/recovery.zip
SHA256 f6fd6d0261a2e617e90111ee955a5cb06196cbfcb352b9da9df265ad27b06883
verified on LAB and hibino; contains wheel, driver, checked input, logs, command.
Dedicated remote directory/archive cleanup follows this commit. Shared mesh
family and system installation are preserved. No manuscript claim changed.
