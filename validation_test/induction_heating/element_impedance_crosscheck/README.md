# Prescribed element-Zs comparison driver

Implementation preparation; no numerical pass or release qualification is
claimed yet. Review `PLAN.md`, obtain compute admission and select the final
exact-source IH artifact before running the solvers. No private model or
measured material data is used by this independently authored case.

`runBoundaryLawGPU.m` solves the complete GPU reaction block twice, with the
same four passive patch impedances on 16/64 flat triangles. It preserves the
case's row-to-patch mapping and exports complex H and q on common quadrature
points. This EM-only experiment fixes thermal state; it is not the nonlinear
temperature-BH sample and does not validate electromagnetic-to-thermal transfer.
Its output JSON is raw operational evidence with the actual loaded MEX path;
keep it outside tracked source and prepare a privacy-checked public record only
after verification.

`boundary_fem.py` implements an independent scattered-A HCurl solve. It does
not use the BEM geometry assembler, BEM matrices or BEM fields in its equations.
The four boundary coefficients remain distinct. The prescribed filament is
integrated independently; reaction work is integrated with composite segment
quadrature of maximum length 2 mm. Read the sign and peak convention in PLAN.

The default command is plan-only and performs no mesh generation or solve:

```powershell
python boundary_fem.py --work <new-job>/plan
```

After admission, run the three declared FEM mesh/domain cases against the
fine BEM sample's `result.json` (use separate fresh job directories):

```powershell
python boundary_fem.py --work <new-job>/fem-coarse --run-index 0 --samples <bem-fine>/result.json --solve
python boundary_fem.py --work <new-job>/fem-fine --run-index 1 --samples <bem-fine>/result.json --solve
python boundary_fem.py --work <new-job>/fem-exterior --run-index 2 --samples <bem-fine>/result.json --solve
python boundary_fem.py --work <new-job>/fem-quadrature --run-index 2 --samples <bem-fine>/result.json --reuse-result <new-job>/fem-exterior/result.json --quadrature-order 24 --solve
python boundary_fem.py --work <new-job>/fem-offset --run-index 2 --samples <bem-fine>/result.json --reuse-result <new-job>/fem-exterior/result.json --probe-offset 5e-10 --solve
```

Quadrature/trace controls reuse the hash-checked saved curved mesh; they never
re-Curve it. An independent mesh must pass check-vol with four named body
patches before a solve. NGSolve owns FE assembly and mapping, sparsecholesky
owns the direct subsolve, and the command owns its TaskManager region.

`compare.py` requires both BEM records, all five FEM records and the exact
native build receipt. It checks artifact/input identities and identical
physical sample positions/weights, then enforces every declared convergence,
phase-preserving field, work and total-power gate. A successful aggregation
is **same-boundary-law** accuracy only; it never marks the release, nonlinear
volume physics or experimental agreement accepted.

Run `python -m pytest test_checks.py -q -o addopts=` for the small phase,
physical-weight, integer-index and material-patch safety contracts. These tests
and static MATLAB analysis are preparation, not execution evidence.
