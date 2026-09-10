# Plain reduced-Omega Kelvin retirement

The public `solve_magnetostatic_reduced_omega_kelvin` entry point is retired.
It raises before assembly rather than silently selecting another formulation.
Use `solve_magnetostatic_mixed_total_reduced_omega_kelvin`, supplying its explicit
material partitions, source-potential traces and Kelvin interface data. This
decision does not retire reduced scalar potentials as a mathematical method.

`reduced_vs_mixed_coarse.json` is the unmodified historical research output
(LAB, Radia 4.95.81, 24,134 elements, linear mu_r=1000). Its `claim` field records
the hypothesis under test, NOT an established result. At p=2, q=22 the reported
energies disagree, with a large relative discrepancy in the Kelvin contribution.
The comparison does not isolate the iron-region variable transformation because
the exterior conventions differ. The cause has not been proved to be a sign bug.

The generating research script is preserved as `reproduce_historical.py`.
It requires the original mesh/source assets and the pre-retirement implementation;
it is not runnable against the retired API and is not a production solver.
Its printed claims are historical hypotheses, not acceptance gates.

The separate ESRF #6 p=3 upper estimate is quadrature-unconverged (approximately
7.946 at q=10 and 8.188 at q=16). Neither that upper estimate nor the derived
interval ratio is certified. Do not use these observations to certify BDM2
containment in all 17 intervals or a general asymptotic 1:2 energy-error ratio.
Re-establish common boundary conventions and quadrature convergence before
testing variable-transformation equivalence or complementary bounds.
