# Far-rule sensitivity diagnostic

HDiv-only diagnostic, not a three-engine validation pass. The checked fine mesh, BDM1, material, coil and direct observation evaluation are unchanged. Both solves request tolerance 1e-10; second requests ho_far_factor=infinity, which bypasses the lower-order far-pair rule. Other integration rules are unchanged, so this is not a general integration convergence certificate.

Relative vector field difference between the two runs: 1.3082353336935784e-8. Reported residuals: 9.601969136684058e-11 (production), 9.601997411887672e-11 (no low-order far rule). Both take 46 iterations. This sensitivity does not explain the approximately 1% medium/fine change. Near/self integration, source projection and discretization remain unresolved.

Runtime verified: Python 3.12.10 at C:/temp/hdiv-farrule-20260924/venv/Scripts/python.exe; isolated Radia 5.0.0 in this venv; NGSolve 6.2.2606 from C:/Program Files/Python312/Lib/site-packages/ngsolve. Same archived CI wheel as baseline. Foreground SSH python -u compare_far_rule_hdiv.py; OMP/MKL/OPENBLAS each 8, driver threads 8. No other compute process at preflight.

Archive including wheel, scripts, meshes, baseline and run.log: S:/Radia/validation_artifacts/hdiv_farrule_20260924/recovery.zip. Remote and LAB SHA256 matched: 5979896a751fa8e11b3a4e14967b8b1715a76e7b56be811e64570447e5674a16. Remove only this job's remote scratch after evidence commit; preserve shared mesh family.
