# BDM2 fourth-level diagnostic

Checked finer mesh, 1688 iron tetrahedra, BDM2, linear mu_r1000,
8 threads, Gram eps1e-14, solve tol1e-10, direct field. No production patches.
32580 DOF, 57.2791479 seconds, 47 iterations, native residual8.157396610e-11,
converged. Legacy self-difference key is not an inter-mesh comparison.

With the final level's core norm as common denominator, medium/fine increment
is 0.2422064% and fine/finer is 0.00429434%; ratio0.0177301.
Thus the last increment contracts strongly, unlike the unusually small
coarse/medium increment. This is numerical mesh sensitivity, not an analytic
error bound. No three-engine pass or matched-error speed claim: independent
host replay and complete comparator evidence remain required.

Archive S:/Radia/validation_artifacts/hdiv_bdm2finer_20260924/recovery.zip
SHA256 80dd394849a843411b7a410b205164d22b1454edd174afe278707f97791cc1aa
verified equal on LAB and hibino. Includes checked input, wheel, driver,
command and logs. Remote dedicated directory/archive cleanup follows commit;
shared mesh family and system installation stay untouched.
