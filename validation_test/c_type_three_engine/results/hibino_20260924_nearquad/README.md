# HDiv near/self Gram quadrature diagnostic

Experimental HDiv-only run, NOT a three-engine pass or a matched-accuracy
performance claim. Same checked fine iron mesh, mu_r=1000, BDM1, 8 threads,
Gram eps=1e-14, CG tolerance=1e-10 and direct observation field as the campaign.
The wrapper overrides _solve.build_charge_gram intorder with None,9,13.
This changes outer quadrature from the symmetric degree-5 rule to product
rules with quad=5,7. Charge-basis assembly also receives quad, and default
inner quadrature tracks it. Far-rule settings and source RHS remain fixed.
Therefore this is not an outer-quadrature-only ablation.

Relative full-vector changes from the previous fine result:
- default: 2.2859223789941344e-10 (tighter solve tolerance control);
- intorder 9: 0.0009490163388221704;
- intorder 13: 0.0009566546961459431.

All three native solves converged; residuals are retained in JSON.
Approximately 0.096% sensitivity is observable, but this single-mesh test
does not establish the origin of the roughly 1% medium/fine increment.
Other mesh levels must be checked before interpreting convergence changes.

Raw recovery (wheel, driver, inputs, log, command):
S:/Radia/validation_artifacts/hdiv_nearquad_20260924/recovery.zip
SHA256 7878a33feeb08de93ea23e92fab84ab16e8cae37bdb18cd9ff0f44c50cfd1b14
verified equal on LAB and hibino. Remote dedicated directory/archive will be
removed after this evidence commit. Shared mesh family and system runtime
are preserved. No manuscript or production source changed.
