# Repaired energy-Newton candidate: Hibino native checks

Status: **PASS for nominal ESRF6 BDM1 three-method acceptance**. Both FEM routes
were re-solved using the same installed candidate as HDiv. This is not a
CI-artifact release approval, BDM2 qualification, or an absolute-error certificate.
Small native checks and the actual ESRF6 BDM1/mass-Riesz residual gate passed.
The finite mesh/CAD/gap/coil
audit is recorded in `../surface_overlap_20260913/`. The execution task owns
the re-solves and scratch cleanup; management owns integration and release.

## Same-wheel three-method result

`case6_fem_repaired_acceptance.json` completed successfully on hibino with
eight threads. The physical source, BH table, observations and installed Python
and native bytes match the repaired HDiv record; the FEM mesh matches the
audited identity. Reduced-A and mixed Omega use order 2. HDiv uses BDM1,
mass-Riesz, the full iron model and no IMA.

| Pair | Core 27-point relative vector RMS | All 45-point relative vector RMS |
| --- | ---: | ---: |
| HDiv / reduced-A | 0.334891% | 0.674624% |
| HDiv / mixed Omega | 0.669467% | 1.196424% |
| reduced-A / mixed Omega | 0.781827% | 1.407877% |

Each relative RMS uses the left field as denominator. The predeclared gate is
3% on the core 27 points, not the raw 45-point set or a pointwise relative
maximum. The replay regression independently recomputes the core values.

Reduced-A converged in 31 iterations with relative change 1.859221e-5;
mixed Omega converged in 14 iterations with relative B change 7.565560e-6.
Both requested 2e-5. These Picard stopping measures are not true-error bounds
and are distinct from HDiv's assembled nonlinear residual (3.278340e-7).
The runtimes were 1355.54 s and 5609.86 s respectively. This run was not an
isolated performance benchmark and must not be advertised as one.

The FEM result SHA-256 is
`66280ab5bda2e3a5a4c0b410b3daa051d011b7446b8949e65330f4ac6d9ac3b8`.
Both linked result JSON files preserve their original bytes through explicit
`-text` attributes; the HDiv evidence bytes are unchanged numerically.
The recovery archive is
`S:/Radia/validation_artifacts/esrf6_mesh_audit_20260913/esrf6-newton-4d85e72cc-recovery.zip`,
SHA-256 `3b198c0a9aba3b6aad02eeecde9014a9eea57f10bef8296bca131a9c8e006808`.
All 34 manifest entries were independently hashed after recovery, and the
installed Python/native identity was rechecked unchanged after the solve.

After evidence commit `c68a2a769`, the job-owned hibino venv, outputs, temporary
profiler, recovery ZIP/helper, and the four recovered input files were deleted.
All named targets were verified absent. An unrelated IH process was detected
and preserved. The durable cleanup report is
`S:/Radia/validation_artifacts/esrf6_mesh_audit_20260913/esrf6_acceptance_cleanup.json`,
SHA-256 `83ca8552d6564c146a1ef4fbbefa14e3d886652e22c3e502431f018b110ed0e8`.
Future runs must explicitly restage the verified inputs from the recovery
archive; do not silently substitute another mesh when old scratch paths vanish.

The follow-up driver requires `--iron-mesh` and `--wheel` in addition to its
HDiv result, FEM mesh/report and output arguments. These explicit inputs must
match the original recorded hashes. Native and all Python source identities
are checked at the current installed package, not at historical absolute paths.
`implementation` retains the original HDiv provenance; `runtime_identity`
records the paths and hashes actually checked in a new run. The input gate
requires literal boolean success, a finite nonnegative residual and a finite
positive tolerance. This runner-only hardening does not modify the completed
result JSON or require another numerical solve.

The final focused selection passed 150 tests, including the numerical-record
replay. The separate installed-wheel deep-saturation regression passed on
hibino (one test, 85.48 s). Neither result waives CI-built release requirements.

This closes the nominal BDM1 field-agreement item. Separate work remains for
iron-sensitive nonlinear observables, BDM2/IMA application qualification,
source-load quadrature convergence and CI-built release acceptance. The later
Hodge quadrature-wiring source change is not included in this candidate wheel.

## Actual ESRF6 BDM1

`case6_bdm1_mass_riesz.json` records the installed candidate on hibino:
Newton 8 iterations, one backtrack, residual 3.278340250386088e-7 at target
2e-5, 62,192 HDiv DoFs, full iron model, no IMA, 8 threads. Gram build took
267.55 s, solve 383.99 s, total including field observations 655.35 s.
These are diagnostic timings, not a controlled performance comparison: another
unrelated `freq8000Hz_mur1.0_toymodel_maxTemp.py` job (PID 12152) existed at the
initial resource probe; later it disappeared and another Python job appeared.
Neither was stopped or modified. BDM1 is not final BDM2 acceptance.

Recovery SHA-256 is `18a5b1c97d28f738191d4c40d1aee1a5d19997ac22824526b4fc3286d2d553ce`.
The durable copy is in `S:/Radia/validation_artifacts/esrf6_mesh_audit_20260913/`.
Cell-average M is retained for iron-side comparisons. Global average M is
almost zero by quadrupole symmetry and is not an adequate nonlinear observable.

The nominal core-field agreement does not certify the nonlinear material over
its full excitation range. A separate iron-sensitive acceptance should compare
section fluxes and material energy using the same constitutive law and physical
integration domain. For a reversible nonlinear B-H law, use the integral of H
with respect to B, not one half of B dot H (the latter is the linear-law
expression). Cell-average M alone cannot reconstruct this energy or pointwise
peak fields. Do not change the running nominal acceptance to add these checks.

## Build and identity

- Source: `4d85e72cce70e02417abcd7972b509dfc2bac632`.
- Branch: `codex/hdiv-energy-newton-consistency-20260913` (pushed).
- Wheel: `radia-4.95.91-cp312-cp312-win_amd64.whl`.
- Wheel SHA256: `e9e6b9f09d105552330ac9b2b6d870187d825d11bf7a78f407553c8c7a866e25`.
- Native SHA256: `932e49626c05c1e3d459f20649632d4bfa0adf9b37f325dace4093d74aeeb7b5`.
- NGSolve/Netgen 6.2.2606; NumPy 2.5.3; SciPy 1.18.1; MKL 2026.1.0.

At candidate manufacture, GitHub workflow dispatch failed twice with HTTP 500
and PR creation failed with HTTP 502. PR #231 was subsequently created and its
source CI passed, but this does not make this wheel a CI-built artifact. It was built in an
isolated LAB venv from the clean committed source with
`Build.ps1 -RequireNativeProvenance`, then `Build_Wheel.ps1 -DryRun`.
The wheel verification passed. The builder's final interactive key prompt was
stopped after verification; `-DryRun -SkipBuild` with `CI=true` independently
repeated verification and exited zero. No package was uploaded to PyPI.

The full build also regenerated the separate Cubit plugin binary. That local
generated change was restored to HEAD; it is not part of this solver repair.
No global/editable environment was redirected and no native binary was manually
copied into an installation. Hibino installed the wheel in the dedicated venv
`C:/temp/esrf6-newton-4d85e72cc/venv`. `preflight.json` checks installed package
bytes against the wheel, and records the explicit non-CI build identity.
`pip check` passed. The version number matches the previous candidate, so
source and wheel hashes, not the version string alone, distinguish this run.

## Actual native tests

`native_smoke.json` and `native_smoke_strong.json` are the original complete
outputs. Their corresponding scripts are retained alongside them. All six
TET/HEX/WEDGE x order-1/order-2 cases passed in each run, using the installed
Radia C++ demag and linear-solve kernels, not the unit tests' dense test double.

The weak 1000 A/m run already met the nonlinear residual target after the
linear warm start, so it did not exercise tangent assembly. The 300000 A/m
run explicitly required at least one tangent assembly and tested actual Newton
iterations. Both runs compared a linear BH table with the linear material solve;
the largest relative difference in average magnetization was `4.12e-7`.

| Strong-drive case | Newton iterations | Final relative residual | Backtracks |
| --- | ---: | ---: | ---: |
| TET order 1 | 13 | 6.90e-6 | 0 |
| TET order 2 | 12 | 2.42e-6 | 2 |
| HEX order 1 | 7 | 7.24e-8 | 0 |
| HEX order 2 | 12 | 1.71e-7 | 6 |
| WEDGE order 1 | 11 | 3.63e-6 | 1 |
| WEDGE order 2 | 19 | 1.58e-5 | 12 |

Every final residual is below the requested `2e-5`. These are small affine
polyhedra, not the ESRF6 magnet or curved-CAD validation. Separate local tests
checked the real NGSolve material weak-form derivatives on quadratic
deformations, distinct materials, and zero-field Hessians. The focused source,
material and tier-policy selection passed 51 tests before wheel manufacture.
The WEDGE order-2 stopping residual is not an estimate of its attainable
accuracy: a tighter tolerance requires another solve, not a prediction that the
case must fail because it stopped near the current threshold.

## Mesh audit handoff

Fresh `fem_check.json`, `iron_check.json` and `input_identity.json` passed
against the same selected iron/FEM bytes as the preceding candidate. They
establish structural, label and sampled-map checks, not all possible surface
overlaps or relative coil/CAD placement. The iron SHA remains `fdd13872...`
and the FEM SHA `dfc12b84...`; the full values are in the identity report.

The actual Hibino iron mesh was recovered to
`C:/temp/esrf6-mesh-audit-20260913/iron_conforming.vol`. The LAB asset mesh has
the same SHA; its adjacent STEP files, manifest and journal were recovered
separately for the management task's geometric audit. That finite audit is now
complete in `../surface_overlap_20260913/`: the final detector found no overlap
in its sampled surface checks, and the selected CAD, coil and gap checks passed.
It is not a proof covering every point of the continuous curved geometry.
The same-wheel mixed Omega rerun has now passed the nominal BDM1 acceptance
above. No production release claim is derived from the native smoke tests or
this finite geometry audit alone.
