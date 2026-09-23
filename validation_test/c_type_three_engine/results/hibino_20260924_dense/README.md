# HDiv Gram backend diagnostic

HDiv-only diagnostic, not a three-engine validation pass. Same checked fine mesh, coil, material, BDM1, direct field evaluator and eight threads as the preceding campaign. Both backends use solve tolerance 1e-10. The driver adapter is reused and exits through a task-owned exception after recording HDiv, before other engines; no missing-engine success is fabricated.

H-matrix versus exact-dense field relative vector difference: 1.1015600186142044e-15. Native exact_dense_normalized_gram is false/true respectively. Reported final solve residuals: 9.601969136684058e-11 and 9.601945432922869e-11; both use 46 iterations. Compression does not explain the approximately 1% medium/fine field increment in this case. Shared integration/discretization or implementation effects remain unresolved; dense is not analytic continuum truth.

Runtime verified: Python 3.12.10 at C:/temp/hdiv-dense-20260924/venv/Scripts/python.exe; Radia 5.0.0 from its isolated venv; NGSolve 6.2.2606 from C:/Program Files/Python312/Lib/site-packages/ngsolve. Same CI wheel as baseline. Foreground SSH python -u compare_dense_hdiv.py with OMP/MKL/OPENBLAS=8 and driver threads=8. Dense memory guard 512 MiB. No other compute at preflight.

Source, wheel, checked inputs, baseline and log recovered in S:/Radia/validation_artifacts/hdiv_dense_20260924/recovery.zip. Remote/LAB SHA256 match: 0326f72bf06737d3fe26e8b115a20ea2e77f05b97f94e7c1f10953609262565f. Task-owned remote scratch cleanup follows evidence commit; shared mesh family remains.
