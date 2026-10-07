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
$script = 'validation_test/induction_heating/bored_workpiece/validate.py'
python $script --work ../bored-workpiece-scratch/coarse --sizes .006 --fem-sizes .0018 --no-carriers --out ../coarse.json
python $script --work ../bored-workpiece-scratch/middle --sizes .0025 --fem-sizes .0025 --no-carriers --out ../middle.json
python $script --work ../bored-workpiece-scratch/fine --sizes .0017 --fem-sizes .003 --out ../fine.json
python $script --merge ../coarse.json ../middle.json ../fine.json
```

Run the three computations sequentially with a process timeout below 1,200
seconds each; the total budget is 3,600 seconds. A failed numerical gate returns
exit code 1 **after saving the evidence**. Inspect the JSON before proceeding;
a timeout or exception can leave incomplete evidence, which must not be merged.
A single-size run is deliberately incomplete and fails the refinement gate.
The merge reevaluates all gates without changing their thresholds. It rejects
incompatible conditions, limits, source/native hashes, versions and duplicate
BEM sizes. The merged `results.json` includes each run's provenance and time.

Meshes remain in the scratch directory and are not committed. `--out` selects
another evidence file. `--sizes` controls BEM sizes; `--fem-sizes` independently
controls the FEM surface sizes (defaults to the BEM sizes). Only the finest
BEM run needs the carrier comparison. The FEM sizes are 1.8, 2.5 and 3 mm; the finest FEM is paired with the
coarse BEM. Every BEM is also compared against that finest FEM in the
reported table. Refinement is assessed from actual node/DOF counts rather
than assuming that both routes refine in the same row order.

The script records source and native hashes, the base checkout commit, package
versions, host, thread count, conditions and gates. Source hashes cover the raw
executed checkout bytes; Git newline normalization can change those hashes
between CRLF and LF checkouts without changing code. No installed packages are
modified by the script.

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

The operating condition was selected to make the circulation effect visible.
An initial exploratory calculation at 10 kHz and relative permeability 100
showed a smaller effect, so the geometry was retained and the condition changed
to 100 kHz and relative permeability 10. The roughly 25% frozen-circulation
error in this example is **condition-specific**, not a general magnitude.
The effect depends on geometry, skin depth and magnetic flux through the bore.
Those two conditions have the same product of frequency and permeability and
thus the same skin depth for this conductivity; skin depth alone does not
explain their different circulation sensitivity. The initial exploratory value
is not presented as a validated FEM comparison.

BEM power balance compares the surface loss with the complete incident-field
reaction, including the multivalued contribution. FEM power balance compares
the surface loss with `(omega/2) Im(a^H f)`, using the independently assembled
load and solution vectors. The observer reads these vectors at solver return
and does not alter assembly or the solution. The real gauge term contributes
no dissipative power.

The gates are fixed in the script: true residual and Faraday residual 1e-6;
unit-jump error 0.005; BEM power balance 2%; FEM power balance 1e-5; total-power
BEM/FEM difference 2%; successive-mesh power change 5%; L2 carrier power width
2% and smaller than the historical point trace; frozen-circulation error at
least 20%; skin-depth/fillet ratio 0.15; per-run runtime 1,200 seconds and summed runtime 3,600 seconds.
A complete study also requires three BEM levels with at least 25% more nodes
per level and at least a 10% increase between the smallest and largest FEM
DOF counts. These coverage checks do not prove convergence.
Carrier width means `(maximum - minimum) / mean` over the four carriers.

## Recorded result (2026-10-07, stricter review)

Overall validation: **FAIL**. The thresholds remain 2% for BEM power
balance and BEM/FEM total-power difference. Failed coarse levels remain in the
record and are not excluded to obtain a pass. `complete_validation` indicates
study coverage, whereas `validation_pass` additionally requires every gate.

| BEM maxh (mm) | BEM nodes | BEM power (W) | BEM balance | FEM maxh (mm) | FEM DOFs | FEM power (W) | Paired difference | Difference to finest FEM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 6 | 2470 | 73.612298 | 2.3363% | 1.8 | 601029 | 72.620293 | 1.3660% | 1.3660% |
| 2.5 | 3778 | 73.630548 | 0.6408% | 2.5 | 486921 | 72.615487 | 1.3979% | 1.3911% |
| 1.7 | 6232 | 73.635503 | 0.0296% | 3 | 517256 | 72.613378 | 1.4076% | 1.3980% |

Differences use FEM power as denominator; balance uses BEM surface power.
The last column uses the same finest FEM reference for all BEM levels.

| FEM maxh (mm) | FEM DOFs | Volume elements | FEM power (W) | Difference to finest FEM |
| --- | ---: | ---: | ---: | ---: |
| 2.5 | 486921 | 90097 | 72.615487 | 0.0066% |
| 3 | 517256 | 96071 | 72.613378 | 0.0095% |
| 1.8 | 601029 | 110509 | 72.620293 | 0.0000% |

The FEM power range divided by finest FEM power is 0.0095%.
BEM balance improves with refinement, but BEM/FEM total-power discrepancy
does not decrease: it stays near 1.4% and slightly increases against the common
finest FEM reference. Refinement therefore supports the balance improvement,
not disappearance of the remaining cross-method discrepancy. This is not an
extrapolated error bound or proof of pointwise field convergence.

| Carrier (r, z), mm | L2 power (W) | Point-trace power (W) | L2 balance | Point-trace balance |
| --- | ---: | ---: | ---: | ---: |
| automatic | 73.635503 | 73.193760 | 0.0296% | 0.2168% |
| 20, -13 | 73.632260 | 73.284052 | 0.0373% | 0.1542% |
| 25, -10 | 73.612248 | 72.776613 | 0.0327% | 0.4383% |
| 19, 13 | 73.543768 | 73.354736 | 0.0448% | 0.0286% |

The carrier comparison uses the finest BEM (6,232 nodes). The power width
`(maximum - minimum) / mean` is 0.7903%
for the point trace and 0.1246% for L2.
The finest BEM frozen-circulation power is 91.100654 W,
which differs from the finest FEM by
25.4479%.

All FEM meshes passed label, curved-mapping and boundary-ownership checks.
Maximum FEM true residual is 4.812e-14;
maximum FEM power-balance discrepancy is 1.128e-06.
Split run times are 627.9, 509.6, 767.3 seconds
(in merge input order); summed solver-run time is 1904.7 seconds.
Each run is below 1,200 seconds and the sum is below 3,600 seconds.
The source base commit in the evidence is the pre-review commit; executed
validation-source hashes identify the revised script before its new commit.

## Limits

Agreement refers to total dissipated power within the stated gates, not to
pointwise field equality. The three-level BEM and independently refined FEM study is a sensitivity check,
not an asymptotic error estimate. BEM and FEM use independently specified
surface mesh sizes. Outer-boundary, filament-segmentation and
quadrature refinement are not included in this bounded run. Nonlinear
materials, off-axis holes, multiple handles, curved BEM panels and thermal
transients are outside its scope.
