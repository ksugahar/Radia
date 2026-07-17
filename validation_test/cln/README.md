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
25 hybrid eddy coordinates.  The response construction first finds 10
HCurl-Krylov directions, then removes 2 directions in the relative null space
of the sampled `curl(T)` current Gram matrix.  The resulting 8 current-Gram
orthonormal bulk modes are combined with 14 conductor-cycle bridge modes and 3
exterior-only SIBC modes.  Two ports are solved together at each frequency
through `_HybridVIMSolve`.

| frequency | RHS-scaled residual | backward error | largest average loss | solver |
|---:|---:|---:|---:|---|
| 100 Hz | 4.693e-16 | 2.577e-20 | 1.218e2 | `radia-cpp-dense-mixed-galerkin` |
| 10 kHz | 4.599e-16 | 2.497e-18 | 1.324e-1 | `radia-cpp-dense-mixed-galerkin` |

The default validation checks reconstruction of all 3557 parent-T coordinates,
all 267 parent BDM1-HDiv coordinates, sampled magnetization, and the
bulk/bridge/SIBC current fields.  The HDiv block solves two physical stator
fields and eight regular-solid-harmonic training fields through degree two.
The default is NGSolve BDM (`HDiv` without `RT=True`).  Actual
Raviart--Thomas is selected only by `--hdiv-family rt`.
Energy POD protects the two physical-response modes and retains six independent
training-complement modes, reducing 267 parent DoFs to 8.  The maximum snapshot
residual is `9.907e-11`; all training responses are reproduced to `3.540e-11`
in magnetic energy.  The bulk current Gram differs from the identity by only
`1.49e-14`.  Adjacency roles retain the 14 conductor-conductor cycle modes and 3
conductor-air SIBC modes, while the 8 ordinary bulk modes are eliminated by the
full-coupled mixed Galerkin step.  The complete BDM1/HCurl system therefore
reduces from 33 to 25 coupled modes.  Its maximum full-coupled Schur error is
`1.110e-16`; port response and Joule loss agree with the direct coupled solve
to `6.982e-16` and `8.140e-16`, respectively.

BDM2 is retained as a production option.  The companion
`hcurl_vim_hdiv_mmm_bdm2_smoke.json`, generated with
`--hdiv-order 2 --multipole-degree 3`, reduces 738 parent DoFs to 15 modes and
reproduces all degree-one/two/three training responses to `3.789e-11` in
magnetic energy.  `hcurl_vim_hdiv_mmm_rt1_smoke.json` and
`hcurl_vim_hdiv_mmm_rt2_smoke.json` are explicit family comparisons.  All four
v8 result files use the same current-Gram and adjacency-driven full-coupled
reduction and report
the active HACApK `_ChargeGramHMatrix` demagnetizing backend:

| HDiv parent | parent DoF -> modes | coupled modes -> retained | max backward error | max Schur error | max port error | max loss error |
|---|---:|---:|---:|---:|---:|---:|
| BDM1 | 267 -> 8 | 33 -> 25 | 2.497e-18 | 1.110e-16 | 6.982e-16 | 8.140e-16 |
| BDM2 | 738 -> 15 | 40 -> 32 | 1.759e-18 | 8.600e-17 | 4.238e-16 | 1.252e-15 |
| RT1 | 369 -> 8 | 33 -> 25 | 9.277e-18 | 7.850e-17 | 9.090e-16 | 2.075e-15 |
| RT2 | 942 -> 15 | 40 -> 32 | 1.575e-18 | 1.147e-16 | 6.657e-16 | 1.647e-15 |

Before current-Gram compression, the bulk block reached condition numbers near
`1e18`.  Removing the `curl(T)` null space lowers the complete coupled
condition number to `2.52e6` at 100 Hz and `2.53e4` at 10 kHz.  Backward error,
port, loss, and field reconstruction checks then all measure the same stable
physical response.

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
  --geometry notched-box --orders 3,4,5,6 `
  --steps 2,4,8,12,16,22 --shifts 0.01,0.1,1,10,100 `
  --output validation_test/cln/evrs_fem_notched_p36.json
```

Against the full p=6 solve:

| case | active DoF | rank | max current L2 error | max Joule-loss error | max corner L2 error | max corner peak error |
|---|---:|---:|---:|---:|---:|---:|
| p=4, n=22 | 1252 | 44 | 3.047e-01 | 7.260e-02 | 2.457e-01 | 4.184e-01 |
| p=6, n=22 | 3557 | 44 | 2.554e-06 | 8.379e-11 | 1.826e-06 | 2.252e-06 |

p=6 depth sweep:

| p | n | active DoF | rank | max current L2 error | max Joule-loss error | max corner L2 error |
|---:|---:|---:|---:|---:|---:|---:|
| 6 | 2 | 3557 | 4 | 5.917e-02 | 1.342e-02 | 4.921e-02 |
| 6 | 4 | 3557 | 8 | 2.571e-03 | 5.220e-05 | 2.159e-03 |
| 6 | 8 | 3557 | 16 | 1.372e-05 | 2.530e-09 | 1.040e-05 |
| 6 | 12 | 3557 | 24 | 2.793e-06 | 9.506e-11 | 1.934e-06 |
| 6 | 22 | 3557 | 44 | 2.554e-06 | 8.379e-11 | 1.826e-06 |

Current interpretation: this is the first direct support for the eddy-bubbling
claim.  On the notched-box stress, the p=6 parent space can be compressed from
3557 active DoFs to 44 EVRS coordinates while preserving the full p=6 current
field and corner-local current to a few parts in one million.  The same rank in
a p=4 parent space misses the full p=6 current by 30.5%, so the high-order
parent space is not cosmetic; it supplies field content that the reduced basis
can keep after eddy bubbles are removed.

## Compute-host curved and broadband validation

The result-bearing sweeps below were run on `mdx`; the JSON files record the
host, runtime, NGSolve version, and wall time.  They supersede the desktop smoke
numbers for paper claims.

```powershell
python validation_test/cln/curved_sphere_geometry_benchmark.py `
  --output validation_test/cln/curved_sphere_geometry_benchmark.json
python validation_test/cln/evrs_current_field_compare.py `
  --geometry sphere --curve-order 4 --orders 3,4,5,6 `
  --steps 2,4,8,12,16,22 --shifts 0.01,0.1,1,10,100 `
  --output validation_test/cln/evrs_fem_sphere_curved_p36.json
```

| curve order | sphere area error | sphere volume error | VIM/FEM area mismatch | tangent defect |
|---:|---:|---:|---:|---:|
| 1 | 5.266e-02 | 9.376e-02 | 5.969e-16 | 5.551e-17 |
| 2 | 1.223e-03 | 1.867e-03 | 1.415e-16 | 1.110e-16 |
| 3 | 1.157e-04 | 1.631e-04 | 2.827e-16 | 1.110e-16 |
| 4 | 5.245e-06 | 7.838e-06 | 0 | 1.665e-16 |

The curved p=6 sphere has 10970 parent DoFs.  At dimensionless shift 1,
ranks 8, 16, and 24 give current errors `1.296e-3`, `1.506e-5`, and
`1.484e-7`.  At shift 100, even rank 44 leaves `5.620e-2`.  This is a
bandwidth limitation of a single low-frequency expansion point, not evidence
that the curved parent space needs still higher order.

## Mixed Galerkin orthogonalization and DtN-SIBC

`HybridVIMSystem.mixed_galerkin_orthogonalization()` constructs separate trial
and test transformations for retained modes `s` and eliminated modes `b`:

```text
P_trial = [-Z_bb^-1 Z_bs; I]
P_test  = [-Z_bb^-T Z_sb^T; I]
P_test^T Z P_trial = Z_ss - Z_sb Z_bb^-1 Z_bs.
```

Thus the projected operator is the exact block Schur complement.  Reciprocal
complex-symmetric VIM systems satisfy `P_trial == P_test`, reducing the method
to operator-Gram orthogonalization.  The curved-sphere rank-8 hybrid case has
63 modes before bulk elimination and 55 after it.  At 100 Hz its condition
number drops from `3.117e9` to `3.155e3`; trial/test orthogonality defects are
below `1.4e-13`, and the Schur identity error is `7.50e-10`.

```powershell
python validation_test/cln/evrs_sibc_mixed_schur.py `
  --geometry sphere --curve-order 4 --order 6 --steps 4,8,12,22 `
  --frequencies 100,1000,10000,100000,1000000 `
  --output validation_test/cln/evrs_dtn_sibc_mixed_sphere_p6.json
```

With graph-cycle and exterior-only SIBC modes included, the rank-8 port error
relative to the run-local rank-44 hybrid reference is `0.2097%` at 100 Hz and
`0.2683%` at 1 MHz.  The nearly flat four-decade response is the practical
reason the DtN-SIBC branch is mandatory for the broadband formulation.

For the coupled BDM-MMM/HCurl formulation the elimination must be performed on
the complete block matrix, not on `Z_e` alone.  Radia records the adjacency
roles `bulk`, `bridge`, and `sibc`, then
`CoupledHDivHybridVIMSystem.solve_frequency_eddy_bubbled()` always retains the
HDiv modes plus conductor-conductor cycle and conductor-air SIBC modes.  Only
the ordinary bulk eddy-bubble block is eliminated.  The reduced matrix is the
exact full coupled Schur complement; a nonzero eliminated RHS is recovered by
the affine lift `A_bb^-1 f_b`.  The production BDM demag block remains backed
by the HACApK `_ChargeGramHMatrix`.

## Corner mesh and Gmsh views

`hcurl_corner_gmsh_visualization.py` exports a Gmsh v4.1 mesh with nine views:
FEM and EVRS current vectors/magnitudes, current error, local EVRS rank, basis
energy density, element size, and corner distance.  OCC re-entrant-edge
refinement with `corner_edge_maxh=0.2` produces 85 tetrahedra and an 8574-DoF
p=6 parent, then reduces it to rank 44.  At shift 10 the global and corner
current errors are `4.855e-3` and `3.666e-3`; corner basis-energy density is
1.426 times the far-field value.

```powershell
python validation_test/cln/hcurl_corner_gmsh_visualization.py `
  --corner-edge-maxh 0.2 --order 6 --steps 22 --shift 10 `
  --output C:\temp\hcurl_gmsh\corner_edge_h02_p6_n44_shift10.msh
```
