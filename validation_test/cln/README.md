# CLN / EVRS Validation

This directory holds validation-class research checks for Cauer ladder network
and Eddy-Visible Response Space work.  These are not fast CI tests and their
desktop runtimes are not benchmark claims.

## EVRS p-by-n convergence smoke

`evrs_pn_convergence.py` builds a unit-box high-order `HCurl(p)` parent space in
NGSolve, compresses it to an EVRS basis with Krylov depth `n`, maps the retained
vectors to sampled `curl(T)` current modes, assembles the reduced VIM system,
and evaluates the external-port admittance

```text
Y_r(s) = B_r^* Z_r(s)^-1 B_r.
```

The generated JSON records DoF compression and the port-admittance error
relative to the highest `(p, n)` case in the run.  Each row also includes
`system_diagnostics`, which records the interaction backend name and the
Hermitian/passive matrix checks used when replacing the dense sampled kernel,
and `topology_diagnostics`, which records exterior conductor faces versus
conductor-conductor loop bridges.  The companion `dof_policy_diagnostics`
block records how those face roles split HCurl DoFs into SIBC-surface,
loop-bridge-protected, and ordinary EVRS candidate sets.  The
`conductor_graph_diagnostics` block records the graph cycle rank, and
`reduction_plan_diagnostics` records the class-wise production plan.  The
`bridge_cycle_basis_diagnostics` and `bridge_cycle_vim_diagnostics` blocks
verify that the graph-cycle bridge class has been converted to an actual
VIM-compatible current basis.

```powershell
python validation_test/cln/evrs_pn_convergence.py
```

## p=6 depth smoke

`evrs_p6_depth12_smoke.json` is a focused p=6 sweep on the same unit-box smoke
problem, with Krylov depth `n=1..12` and the `n=12` port admittance used as the
run-local reference.  It is still a desktop smoke, not a benchmark or final
accuracy claim.

| n | rank / active DoF | kept | removed | rel. error to n=12 |
|---:|---:|---:|---:|---:|
| 1 | 2 / 1226 | 0.163% | 99.837% | 7.602e-01 |
| 2 | 4 / 1226 | 0.326% | 99.674% | 6.376e-01 |
| 3 | 6 / 1226 | 0.489% | 99.511% | 5.269e-01 |
| 4 | 8 / 1226 | 0.653% | 99.347% | 2.826e-01 |
| 5 | 10 / 1226 | 0.816% | 99.184% | 3.123e-01 |
| 6 | 12 / 1226 | 0.979% | 99.021% | 2.620e-01 |
| 7 | 14 / 1226 | 1.142% | 98.858% | 1.919e-01 |
| 8 | 16 / 1226 | 1.305% | 98.695% | 9.016e-02 |
| 9 | 18 / 1226 | 1.468% | 98.532% | 5.652e-02 |
| 10 | 20 / 1226 | 1.631% | 98.369% | 5.255e-02 |
| 11 | 22 / 1226 | 1.794% | 98.206% | 7.964e-04 |
| 12 | 24 / 1226 | 1.958% | 98.042% | 0.000e+00 |

Current interpretation: for this two-port smoke, p=6 needs roughly 20--24 EVRS
coordinates to converge the reduced port admittance, so the useful reduction is
about 98%.  Smaller ranks are valuable as low-order CLN/Cauer fits but should
not be called converged.

## EVRS + bridge-cycle + surface-Omega/SIBC mixed Schur smoke

`evrs_sibc_mixed_schur.py` uses `NgsolveTopologyAwareHybridVIM` to add
graph-cycle bridge modes and three surface-Omega modes on the same unit box,
then assembles the mixed EVRS/bridge/surface VIM matrix

```text
Z(s) = R + sL + Z_s(s) M_Gamma,
```

and checks the IGTE mixed-Galerkin Schur complement by eliminating the bulk EVRS
block while keeping the bridge-cycle and surface blocks.  The validation
records both the port-admittance error relative to the highest Krylov depth in
the run and the residual of the Schur-reconstructed solution in the original
mixed system.  The saved `system_diagnostics` block is the backend replacement
gate for the mixed SIBC case, while
`topology_diagnostics` records whether the surface/SIBC basis is being placed
only on exterior conductor faces.  `dof_policy_diagnostics` records the HCurl
DoFs touched by those SIBC faces and by conductor-conductor loop bridges.
`conductor_graph_diagnostics` records the bridge cycle basis, and
`reduction_plan_diagnostics` additionally records the mixed-mode estimate from
cycle-basis bridge modes, EVRS rank, and surface-mode count.  The bridge cycle
basis diagnostics record its mass matrix and standalone VIM passivity check.
The `topology_aware_hybrid_vim_diagnostics["eddy_bubbling"]` block is the
production summary of the same split: SIBC surface, non-SIBC trace,
bridge/cycle, and ordinary bulk eddy-bubble candidates.
The `topology_aware_hybrid_vim_diagnostics` block is the production-builder
summary tying the three retained classes to one reduced VIM system.  It also
contains a `parent_order_ledger` block when the run supplies the symbolic
bulk/bridge/surface degree requirements, so `p=6` can be identified as required
or merely conservative.

```powershell
python validation_test/cln/evrs_sibc_mixed_schur.py
python validation_test/cln/evrs_sibc_mixed_schur.py --orders 4,5,6 --steps 8,12
```

Current p=6 desktop smoke (`evrs_sibc_mixed_schur_smoke.json`) with
cycle-basis bridge reduction:

| n | bulk rank / active DoF | bridge DoF -> cycle modes | class-wise estimate | max rel. port error vs n=12 | max Schur residual |
|---:|---:|---:|---:|---:|---:|
| 8 | 16 / 1226 | 386 -> 7 | 26 | 4.759e-02 | 1.501e-14 |
| 11 | 22 / 1226 | 386 -> 7 | 32 | 1.085e-03 | 2.465e-13 |
| 12 | 24 / 1226 | 386 -> 7 | 34 | 0.000e+00 | 3.651e-13 |

The direct mixed impedance can be extremely ill-conditioned in this smoke
problem, so the Schur validation uses equation residual rather than direct
solution-vector difference as the pass/fail observable.

## HCurl-VIM / HDiv-MMM end-to-end smoke

`hcurl_vim_hdiv_mmm_end_to_end.py` runs the production-facing one-call API,
including the final native reduced solve and field reconstruction:

```powershell
python validation_test/cln/hcurl_vim_hdiv_mmm_end_to_end.py
```

The saved p=6 notched-box result
`hcurl_vim_hdiv_mmm_end_to_end_smoke.json` reduces 3557 parent HCurl DoFs to
33 hybrid eddy coordinates: 16 bulk EVRS modes, 14 conductor-cycle bridge
modes, and 3 exterior-only SIBC modes.  Two ports are solved together at each
frequency through `_HybridVIMSolve`.

| frequency | mixed residual | largest average loss | solver |
|---:|---:|---:|---|
| 100 Hz | 2.550e-14 | 1.219e2 | `radia-cpp-dense` |
| 10 kHz | 7.242e-15 | 1.320e-1 | `radia-cpp-dense` |

The default validation checks reconstruction of all 3557 parent-T coordinates,
all 267 parent BDM1-HDiv coordinates, sampled magnetization, and the
bulk/bridge/SIBC current fields.  The HDiv block solves two physical stator
fields and eight regular-solid-harmonic training fields through degree two.
The default is NGSolve BDM (`HDiv` without `RT=True`).  Actual
Raviart--Thomas is selected only by `--hdiv-family rt`.
Energy POD protects the two physical-response modes and retains six independent
training-complement modes, reducing 267 parent DoFs to 8.  The maximum snapshot
residual is `9.907e-11`; all training responses are reproduced to `3.540e-11`
in magnetic energy.

BDM2 is retained as a production option.  The companion
`hcurl_vim_hdiv_mmm_bdm2_smoke.json`, generated with
`--hdiv-order 2 --multipole-degree 3`, reduces 738 parent DoFs to 15 modes and
reproduces all degree-one/two/three training responses to `3.789e-11` in
magnetic energy.  `hcurl_vim_hdiv_mmm_rt1_smoke.json` and
`hcurl_vim_hdiv_mmm_rt2_smoke.json` are explicit family comparisons.

## Planar HDiv-MMM BDM/RT corner smoke

`planar_hdiv_mmm_response_smoke.py` validates the 2-D response-reduction path
on an L-shaped body, comparing independently solved parent coefficients,
element-average magnetization, and the re-entrant-corner neighborhood.

```powershell
python validation_test/cln/planar_hdiv_mmm_response_smoke.py
python validation_test/cln/planar_hdiv_mmm_response_smoke.py `
  --order 2 --harmonic-degree 3 `
  --output validation_test/cln/planar_hdiv_mmm_bdm2_response_smoke.json
```

| HDiv parent | parent DoF | response modes | max energy error | max corner error |
|---|---:|---:|---:|---:|
| BDM1, degree 2 | 42 | 4 | 4.906e-12 | 4.186e-12 |
| BDM2, degree 3 | 93 | 6 | 3.299e-11 | 2.298e-11 |
| RT1, degree 2 | 62 | 4 | 2.893e-13 | 1.062e-14 |
| RT2, degree 3 | 123 | 6 | 7.064e-12 | 3.935e-11 |

These response-span checks establish that BDM and explicit RT use the same
protected-POD path; they do not claim that a fixed rank covers every motor
geometry or excitation.

## p=4/5/6 parent-order smoke

`evrs_sibc_p456_depth20_smoke.json` is the first numerical attack on the
question "is p=6 really needed?"  The run compares `p=4,5,6` against the
highest case in the sweep, `p=6,n=20`.  This is still a smooth unit-box desktop
smoke, not a motor-grade benchmark, but it is enough to reject "p=6 is always
required" as a blanket statement for this case.

| p | n | active DoF | reduced modes | max rel. error vs p=6,n=20 | max Schur residual |
|---:|---:|---:|---:|---:|---:|
| 4 | 14 | 428 | 38 | 3.467e-02 | 2.068e-11 |
| 4 | 18 | 428 | 46 | 2.468e-03 | 2.670e-10 |
| 4 | 20 | 428 | 50 | 4.661e-03 | 2.644e-10 |
| 5 | 14 | 758 | 38 | 4.918e-02 | 3.350e-12 |
| 5 | 18 | 758 | 46 | 3.290e-03 | 7.489e-11 |
| 5 | 20 | 758 | 50 | 6.289e-03 | 9.899e-11 |
| 6 | 14 | 1226 | 38 | 6.683e-02 | 2.074e-12 |
| 6 | 18 | 1226 | 46 | 7.118e-03 | 4.790e-11 |
| 6 | 20 | 1226 | 50 | 0.000e+00 | 3.006e-10 |

Current interpretation: for this smooth smoke, `p=4,n=18` is already within
0.25% of the `p=6,n=20` reference while using about one third of the active
parent DoFs.  The p-ledger for the default bulk/bridge/surface requirements has
`required_parent_order=4`; therefore `p=6` is admissible but conservative here.
The non-monotone rows warn that the reduced basis and reference depth must be
treated together: this table does not prove p=4 is universally sufficient for
motors, corners, thin skin depth, or higher-order SIBC traces.

The companion `evrs_sibc_p6_depth20_smoke.json` checks reference stability on
the same p=6 parent space: `n=18` is still about 0.71% away from `n=20`, while
`n=19` is about 0.27% away.  Thus the p-comparison should use a sufficiently
deep EVRS/CLN reference; the older `p=6,n=12` smoke is not a final p-reference.

## notched-box corner stress smoke

`evrs_sibc_mixed_schur.py` also supports `--geometry notched-box`, which removes
a corner from the unit box and creates a re-entrant exterior corner while
keeping the same conductor/SIBC classification rule.  The smoke below uses one
frequency point and compares against `p=6,n=22`.

```powershell
python validation_test/cln/evrs_sibc_mixed_schur.py `
  --geometry notched-box --orders 4,6 --steps 20,22 `
  --frequencies 100 --output validation_test/cln/evrs_sibc_notched_p46_depth22_smoke.json
python validation_test/cln/evrs_sibc_mixed_schur.py `
  --geometry notched-box --orders 5,6 --steps 20,22 `
  --frequencies 100 --output validation_test/cln/evrs_sibc_notched_p56_depth22_smoke.json
```

| geometry | p | n | active DoF | reduced modes | max rel. error vs p=6,n=22 | max Schur residual |
|---|---:|---:|---:|---:|---:|---:|
| notched-box | 4 | 22 | 1252 | 61 | 2.332e-02 | 1.053e-09 |
| notched-box | 5 | 22 | 2207 | 61 | 3.273e-02 | 1.671e-09 |
| notched-box | 6 | 22 | 3557 | 61 | 0.000e+00 | 3.911e-07 |

Current interpretation: unlike the smooth box, the notched-box stress leaves a
2--3% p-error at p=4/5 even after the EVRS depth is raised to 22.  This does
not prove p=6 is universally necessary, but it gives the first numerical case
where the conservative p=6 parent space is doing real work.  The p=5 row is not
monotone in this smoke, so parent-order choice must still be judged by retained
observables, not by p alone.

## current-field comparison against full p=6

`evrs_current_field_compare.py` checks the current field itself, not only port
admittance.  It solves the dimensionless parent HCurl system
`(K + shift M) T = b`, samples `J = curl(T)`, and compares reduced EVRS
currents against the full p=6 parent solve on the same mesh and quadrature
points.  The reported energy norm is the weighted current L2 norm; Joule loss
uses `int |J|^2 / sigma dV`.  For `notched-box`, the corner metric is measured
near the re-entrant line at `(x,y)=(0.45,0.45)`.

```powershell
python validation_test/cln/evrs_current_field_compare.py `
  --geometry notched-box --orders 4,6 --steps 22 --shifts 1 `
  --output validation_test/cln/evrs_current_field_notched_p46_n22_smoke.json
python validation_test/cln/evrs_current_field_compare.py `
  --geometry notched-box --orders 6 --steps 2,4,8,12,22 --shifts 1 `
  --output validation_test/cln/evrs_current_field_notched_p6_depth_smoke.json
```

Against the full p=6 solve:

| case | active DoF | rank | max current L2 error | max Joule-loss error | max corner L2 error | max corner peak error |
|---|---:|---:|---:|---:|---:|---:|
| p=4, n=22 | 1252 | 44 | 2.101e-01 | 1.545e-01 | 1.567e-01 | 1.387e-01 |
| p=6, n=22 | 3557 | 44 | 1.605e-06 | 8.435e-08 | 1.112e-06 | 9.454e-07 |

p=6 depth sweep:

| p | n | active DoF | rank | max current L2 error | max Joule-loss error | max corner L2 error |
|---:|---:|---:|---:|---:|---:|---:|
| 6 | 2 | 3557 | 4 | 4.578e-02 | 6.011e-02 | 3.590e-02 |
| 6 | 4 | 3557 | 8 | 1.896e-03 | 1.100e-03 | 1.326e-03 |
| 6 | 8 | 3557 | 16 | 1.041e-05 | 3.123e-06 | 7.871e-06 |
| 6 | 12 | 3557 | 24 | 1.791e-06 | 1.161e-07 | 1.172e-06 |
| 6 | 22 | 3557 | 44 | 1.605e-06 | 8.435e-08 | 1.112e-06 |

Current interpretation: this is the first direct support for the eddy-bubbling
claim.  On the notched-box stress, the p=6 parent space can be compressed from
3557 active DoFs to 44 EVRS coordinates while preserving the full p=6 current
field, Joule loss, and corner-local current to about 1e-6.  The same rank in a
p=4 parent space misses the full p=6 current by 10--20%, so the high-order
parent space is not cosmetic; it supplies field content that the reduced basis
can keep after eddy bubbles are removed.
