# Stepped bored workpiece: linear SIBC validation

This example uses a synthetic, axisymmetric workpiece with a stepped axial
through-hole. It compares total dissipated power from the production P1 scalar
BIE with one circulation unknown and Faraday closure against an independent
HCurl air-domain FEM with SIBC. The FEM assembles curl-curl and Robin forms and
fixed filament loads; it does not use the BEM operators.

## Reproduce

Use a Radia build with its native extensions available in `src/radia` and the
matching NGSolve environment, plus `cubit-mesh-export` for the mesh checker
(no Cubit process is used). From the repository root, choose a disposable
scratch directory outside the checkout:

```powershell
$env:OPENBLAS_NUM_THREADS = '2'
$env:OMP_NUM_THREADS = '2'
python validation_test/induction_heating/bored_workpiece/validate.py --work ../bored-workpiece-scratch
```

The default run writes `results.json` beside the script. Meshes remain in the
scratch directory and are not committed. `--out` selects another evidence file.
`--sizes` and `--no-carriers` are diagnostic options; a reduced run is not a
complete validation. The script records source and native hashes, the base
checkout commit, package versions, host, thread count, conditions and gates.
Source hashes cover the raw executed checkout bytes; Git newline normalization
can change those hashes between CRLF and LF checkouts without changing code.
No installed packages are modified by the script.

## Conditions and interpretation

- Workpiece: outer radii 30/24 mm and bore radii 10/14 mm, with both steps at
  z = 0; height 50 mm, centered at z = 0. All circular corners have 1.8 mm
  fillets. OCC geometry is authored directly in metres.
- Excitation: three closed 128-segment polygons approximating radius-40 mm
  circular filaments at z = -18, 0, +18 mm; each carries a prescribed 100 A
  peak phasor. There is no circuit-current solve.
- Linear material: relative permeability 10, conductivity 5e6 S/m, 100 kHz;
  peak phasors use exp(+i omega t). Skin depth is 0.2251 mm, or 0.1251 of
  the minimum fillet radius. The two methods share the same linear SIBC
  approximation; this does not validate that approximation against a resolved
  conducting volume.
- BEM: flat P1 surface, dense in-tree operators, singular quadrature 6 and
  regular degree 7; maximum 7,000 nodes. Surface-Poisson reconstruction of the
  incident scalar potential uses the fixed filament field.
- FEM: air sphere of radius 180 mm, A = 0 at the outer boundary, HCurl order 2,
  sparse Cholesky and real gauge regularization 1e-6. The live OCC mesh is
  curved to order 2 before saving; the loaded mesh is never re-curved.
  FEM curvature-safety factor is 1.5. Area and volume are recorded before
  and after curving. A versioned label contract, curved-element mapping
  quality and domain-ownership checks must pass before the FEM solve.
- Carrier study on the finer BEM mesh: automatic carrier and explicit
  (radius, z) = (20,-13), (25,-10), (19,13) mm. Explicit rings are checked
  for material containment and segment clearance. The historical point trace
  uses area-weighted vertex normals only within a temporary comparison patch;
  production source files are unchanged.

BEM power balance compares the surface loss with the complete incident-field
reaction, including the multivalued contribution. FEM power balance compares
the surface loss with `(omega/2) Im(a^H f)`, using the independently assembled
load and solution vectors. The observer reads these vectors at solver return
and does not alter assembly or the solution. The real gauge term contributes
no dissipative power.

The gates are fixed in the script: true residual and Faraday residual 1e-6;
unit-jump error 0.005; BEM power balance 10%; FEM power balance 1e-5; total-power
BEM/FEM difference 5%; successive-mesh power change 5%; L2 carrier power width
2% and smaller than the historical point trace; frozen-circulation error at
least 20%; skin-depth/fillet ratio 0.15; full runtime 1,200 seconds.
Carrier width means `(maximum - minimum) / mean` over the four carriers.

## Recorded result (2026-10-07)

All default-run gates passed on the recorded runtime. Total-power agreement
is within 5% on both meshes; freezing circulation exceeds the 20% error gate.

| maxh (mm) | BEM nodes | FEM DOFs | BEM power (W) | FEM power (W) | Difference | Frozen power (W) | Frozen error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 6 | 2470 | 473515 | 73.612298 | 72.616185 | 1.372% | 90.981879 | 25.291% |
| 4.5 | 2602 | 462099 | 73.635849 | 72.617438 | 1.402% | 90.981935 | 25.289% |

| Carrier (r, z), mm | L2 power (W) | Point-trace power (W) | Clearance (mm) |
| --- | ---: | ---: | ---: |
| automatic | 73.635849 | 73.117080 | 9.826 |
| 20, -13 | 73.641555 | 73.116722 | 9.920 |
| 25, -10 | 73.624483 | 72.659361 | 4.930 |
| 19, 13 | 73.536148 | 73.278573 | 4.907 |

The carrier power width decreases from 0.8477%
(point trace) to 0.1432% (L2 projection).
Both FEM meshes have zero invalid mapping samples and pass label and ownership
checks. Maximum FEM true residual is 4.714e-14 and its maximum
power-balance discrepancy is 1.058e-06. Across both BEM meshes
and every carrier/method comparison, the maximum true residual is
8.210e-16, Faraday residual 6.998e-16, unit-jump error
1.129e-06, and power-balance discrepancy 2.531%.
Successive-mesh total-power changes are 0.0320% (BEM)
and 0.0017% (FEM). The complete run took
932.2 seconds (15.54 minutes).

## Limits

Agreement refers to total dissipated power within the stated gates, not to
pointwise field equality. The two-mesh comparison is a sensitivity check,
not an asymptotic error estimate. Meshing is curvature dominated: the
smaller nominal surface size increases the BEM vertex count, but does not
imply a monotonic increase in FEM degrees of freedom. Outer-boundary, filament-segmentation and
quadrature refinement are not included in this bounded run. Nonlinear
materials, off-axis holes, multiple handles, curved BEM panels and thermal
transients are outside its scope.
