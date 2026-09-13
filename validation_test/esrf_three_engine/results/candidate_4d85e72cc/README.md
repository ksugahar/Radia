# Repaired energy-Newton candidate: Hibino native checks

Status: **HOLD for ESRF6 acceptance**. The repaired candidate passed the small
native tests below. ESRF6 itself was **not re-solved**: the latest user request
prioritizes a fresh mesh/CAD/gap/coil-placement audit, owned by the management
task. Existing mesh checks do not replace that audit.

## Build and identity

- Source: `4d85e72cce70e02417abcd7972b509dfc2bac632`.
- Branch: `codex/hdiv-energy-newton-consistency-20260913` (pushed).
- Wheel: `radia-4.95.91-cp312-cp312-win_amd64.whl`.
- Wheel SHA256: `e9e6b9f09d105552330ac9b2b6d870187d825d11bf7a78f407553c8c7a866e25`.
- Native SHA256: `932e49626c05c1e3d459f20649632d4bfa0adf9b37f325dace4093d74aeeb7b5`.
- NGSolve/Netgen 6.2.2606; NumPy 2.5.3; SciPy 1.18.1; MKL 2026.1.0.

GitHub workflow dispatch failed twice with HTTP 500 and PR creation failed
with HTTP 502. This is not a CI-built or CI-approved wheel. It was built in an
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

## Mesh audit handoff

Fresh `fem_check.json`, `iron_check.json` and `input_identity.json` passed
against the same selected iron/FEM bytes as the preceding candidate. They
establish structural, label and sampled-map checks, not all possible surface
overlaps or relative coil/CAD placement. The iron SHA remains `fdd13872...`
and the FEM SHA `dfc12b84...`; the full values are in the identity report.

The actual Hibino iron mesh was recovered to
`C:/temp/esrf6-mesh-audit-20260913/iron_conforming.vol`. The LAB asset mesh has
the same SHA; its adjacent STEP files, manifest and journal were recovered
separately for the management task's geometric audit. The full ESRF6 rerun
must wait for that audit. No successful field-comparison notebook or production
release claim is derived from these smoke tests.
