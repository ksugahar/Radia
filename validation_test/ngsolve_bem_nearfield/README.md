# BEM near-field validation

Run `python validation_test/ngsolve_bem_nearfield/run_square_sheet.py --output result.json`
in an isolated NGSolve environment. No Radia native extension is needed.

A unit-density square sheet of side 2 is evaluated over a parallel target patch
at distances 0.1, 0.001 and 1e-6. The reference is the analytic rectangle integral
of the Laplace single-layer kernel 1/(4*pi*r). Both direct evaluation and the
local-expansion path must satisfy the recorded relative L2 thresholds. The script
owns its TaskManager and fixes one thread; this is an accuracy gate, not timing.

NGSolve 6.2.2606 and 6.2.2607 both passed this flat-panel case during migration
review. It does not establish an accuracy gain from upgrading, curved-panel
accuracy, derivative-kernel accuracy, or Radia-IH application acceptance. Those
require separate representative validation. The 2607 implementation includes
near-element quadrature and tangent corrections in `bem/potentialcf.cpp`.