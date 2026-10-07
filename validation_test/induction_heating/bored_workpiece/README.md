# Stepped bored workpiece: linear SIBC validation

Synthetic geometry and closed prescribed filament currents; independent HCurl
air-domain FEM and P1 scalar BIE with one circulation unknown. The production
loop closure now uses distributed Galerkin electric and magnetic work. The
legacy residual key `faraday_residual_rel` measures that augmented row.

## Reproduce

Use matching native extensions in `src/radia`, NGSolve and the mesh checker.
Run from the repository root. No installed package is modified.

```powershell
$env:OPENBLAS_NUM_THREADS = '2'
$env:OMP_NUM_THREADS = '2'
$script = 'validation_test/induction_heating/bored_workpiece/validate.py'
python $script --work ../bored-scratch/coarse --sizes .006 --fem-sizes .0018 --no-carriers --out ../coarse.json
python $script --work ../bored-scratch/middle --sizes .0025 --fem-sizes .0025 --no-carriers --out ../middle.json
python $script --work ../bored-scratch/fine --sizes .0017 --fem-sizes .003 --out ../fine.json
python $script --merge ../coarse.json ../middle.json ../fine.json
python validation_test/induction_heating/bored_workpiece/check_gates.py
```

Each computation has a 1,200 s budget; summed run time must be at most 3,600 s.
A single-level result fails coverage by design; it is not a complete study.
To reproduce the full diagnostic data on this runtime, allow execution beyond
1,200 seconds while retaining the time gates and their failure verdicts.
Incomplete or timed-out files must not be merged. Meshes and logs remain scratch
artifacts. Results record executed source/native hashes, versions, host and time.

## Fixed conditions and acceptance

- Metre-authored stepped cylinder: outer radii 30/24 mm, bore radii 10/14 mm,
  height 50 mm centred at zero, step at zero, circular fillets 1.8 mm.
- Three closed radius-40 mm, 128-segment filaments at z = -18, 0, 18 mm,
  each 100 A peak; exp(+i omega t), 100 kHz, permeability 10, conductivity 5e6 S/m.
  Skin depth 0.2251 mm; depth/fillet ratio 0.1251.
- BEM: undeformed flat P1, dense in-tree operators, singular quadrature 6,
  regular degree 7; distributed P0 work uses production quadrature bonus 4.
- FEM: HCurl order 2, sphere radius 180 mm, outer A=0, real gauge 1e-6,
  sparse Cholesky. The live SI CAD mesh is curved to order 2 before saving.
  Label, curved-mapping and ownership checks precede every solve.

Thresholds are unchanged: BEM balance <=2%, BEM/FEM difference <=2%,
linear/augmented-row residual <=1e-6, unit jump error <=0.005, FEM balance
<=1e-5, successive power change <=5%, L2 carrier width <=2% and smaller than
point width, frozen-circulation error >=20%, skin-depth/fillet <=0.15.
At least three BEM levels (node growth >=25% per level), FEM DOF range >=10%,
and BEM nodes <=7,000 are required. BEM balance must strictly decrease with
refinement. The absolute 2% balance and cross-method gates apply only to the
finest BEM and independently finest FEM; all-level residual, FEM and quality
gates remain required. Carrier checks apply to the finest BEM.

## Results (2026-10-08)

Overall validation: **FAIL** under the unchanged gates.

The finest-solution and carrier gates pass. The failed gates are
`balance_monotonic`, `runtime`, and `total_runtime`. Balance is
0.065102%, 0.069803%, 0.054971%: the middle level increases, so strict
monotonicity fails. This does not change the 2% finest-solution thresholds.
The middle computation was stopped after 1,200 seconds; its completed BEM
was retained and the FEM retried separately. The diagnostic finest run was
allowed to finish beyond the execution budget. Both exceed the unchanged
time gates; retry time is included. `execution_events` records these events.
No failed attempt is relabelled as a passing run.

| BEM maxh mm | Nodes | BEM W | Balance % | FEM maxh mm | FEM DOFs | FEM W | Difference to finest FEM % |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 6 | 2470 | 73.482534 | 0.065102 | 1.8 | 601029 | 72.620293 | 1.187328 |
| 2.5 | 3778 | 73.576531 | 0.069803 | 2.5 | 486921 | 72.615487 | 1.316763 |
| 1.7 | 6232 | 73.634332 | 0.054971 | 3 | 517256 | 72.613378 | 1.396357 |

Finest-solution comparison with the previous recorded cut-line closure:

| Quantity | Previous | Distributed work | Change |
| --- | ---: | ---: | ---: |
| BEM power W | 73.635503 | 73.634332 | -0.001171 |
| Balance % | 0.029599 | 0.054971 | +0.025373 |
| L2 carrier width % | 0.124630 | 0.118411 | -0.006219 |
| Point carrier width % | 0.790301 | 1.340039 | +0.549738 |
| Difference to finest FEM % | 1.397970 | 1.396357 | -0.001612 |
| Frozen circulation power W | 91.100654 | 91.100654 | 0.000000 |
| Frozen circulation error % | 25.447929 | 25.447929 | 0.000000 |

Percentage changes above are percentage points. The previous independent
synthetic results are identified by their validation commit in `results.json`.

| Carrier radius,z mm | L2 power W | Point power W | L2 balance % | Point balance % |
| --- | ---: | ---: | ---: | ---: |
| automatic | 73.634332 | 73.134926 | 0.054971 | 0.075635 |
| 20, -13 | 73.633926 | 73.180752 | 0.057377 | 0.062754 |
| 25, -10 | 73.609551 | 72.741825 | 0.045069 | 0.197023 |
| 19, 13 | 73.547174 | 73.722667 | 0.038760 | 0.090746 |

The point-trace comparison still has a distinct meaning: the patched
area-weighted vertex-normal trace enters the BIE lift and the new distributed
normal-flux work. Both routes use the new closure; this is not a reproduction
of the previous closure. Explicit rings at (20,-13), (25,-10), (19,13) mm
and the automatic ring pass containment/clearance checks. Production files
are unchanged by the temporary patch.

FEM balance compares surface loss with (omega/2) Im(a^H f). BEM balance
compares loss with the complete incident-field reaction including circulation.
A small internal balance does not bound cross-method or pointwise field error.
Refinement is a sensitivity study, not an extrapolated accuracy guarantee.
Both routes share linear SIBC; no skin-resolved conducting-volume validation
is claimed. The frozen-circulation diagnostic fixes alpha=0; its error is condition-specific.
This operating point was selected after a preliminary 10 kHz, permeability-100
study showed a smaller circulation effect. That exploratory run is not a
validated FEM comparison. Equal frequency-permeability products give equal
skin depths, but do not imply equal circulation sensitivity. Outer-boundary,
filament and quadrature refinement, nonlinear material, off-axis/multiple holes
and thermal evolution are outside this study.

Solver-run times: 1029.7 s, 1654.0 s, 3497.8 s; sum 6181.4 s.
Source hashes identify the exact executed script and loop-work implementation.

`source_dirty` concerns the library `src` directory. The recorded validator
hash identifies the executed, uncommitted diagnostic additions before this
evidence commit. Native kernels came from the installed 5.2.3 distribution;
the copied native binaries and executed Python source hashes are recorded.
