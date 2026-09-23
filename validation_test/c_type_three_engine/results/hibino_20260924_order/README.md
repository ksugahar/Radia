# Coarse HDiv order sensitivity

HDiv-only diagnostic, not a three-engine pass or performance certificate.
Same checked coarse mesh and source, linear mu_r1000, 8 threads,
Gram eps1e-14, solve tol1e-10, direct field. Orders 1 and 2 use their
respective default quadrature/preconditioning policies, not a pure basis-only ablation.
`rule` denotes HDiv order in the JSON; the historical difference key means
same-mesh BDM1 baseline, as the JSON note records.

BDM1: 3090 DOF, 3.1011682 seconds, 46 iterations, residual9.429446368e-11.
BDM2: 8784 DOF, 13.8746463 seconds, 49 iterations, residual5.558126381e-11.
Both converged. Full observation vector difference 0.007263346024180067.
This is order sensitivity, not proof that BDM2 is accurate. Next: measure
BDM2 mesh increments. No mixed-Omega speed ratio is claimed.

Archive S:/Radia/validation_artifacts/hdiv_order_20260924/recovery.zip
SHA256 acf492ede3c295aa3e64896f01edb55b3d3946a518fda5a78a05efc52545cee0
verified equal on LAB and hibino. Includes input, wheel, driver, log, command.
Dedicated remote directory/archive cleanup follows this commit; shared family
and system installation remain unchanged.
