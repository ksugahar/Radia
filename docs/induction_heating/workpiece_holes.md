# Workpieces with a through-hole

Use the actual workpiece mesh, including the bore. For a linear SIBC analysis,
Radia detects surface topology and adds the circulating-current cohomology
unknown automatically. No cut, section anchor or auxiliary ring needs to be
specified for a supported workpiece.

```powershell
python -m radia.panels.calc_inductance --coil-solver peec --coil-step coil.step --vol work.vol --wp-label work --frequency 7000 --current 100 --sigma 6e6 --mu-r 200 --output result.json
```

Current is peak amperes. The defaults `--wp-loop-dof auto` and
`--wp-bem-backend auto` select HACApK for a simply connected workpiece and the
P1 dense loop-extended BEM for a workpiece with one handle. Explicit backend
requests are preserved; an incompatible request raises instead of dropping
the circulating current. Explicit `--wp-bem-backend hacapk` also supports
one handle, using the same distributed loop closure in a bordered GMRES
operator. The log reports the topology and backend selection.

The separate `--wp-loop-work-backend auto` option uses native NGSolve
Galerkin FMM for the P0 magnetic loop work at 512 or more surface faces,
with native direct products below that count. `dense` and `fmm` select
explicit routes. The Python equivalents are `loop_work_backend` on
`ScalarBIESIBCSolver` and `wp_loop_work_backend` on the two coupled solvers.
FMM compresses separated interactions of the Galerkin quadrature; direct
near/singular integration, the P1 BIE, normal-flux energy and electric work
are unchanged. It does not replace each triangle by a centroid source.
Results include `wp_loop_work_diagnostics` with the selected route,
quadrature bonus, FMM controls and NGSolve version. FMM uses minimum order
20, direct leaf capacity 100 and separation factors 2/3. Construction/apply
failure raises with the `dense` option and never changes routes silently.

The self-authored ring compression study is
`validation_test/induction_heating/run_loop_work_compression.py`. Run the
dense and FMM cases with identical meshes and thread controls to compare
heat, fields, lift sensitivity, reaction, residuals, time and process peak.
For example, run `--level 64 --backend dense --output dense.json`, then
`--level 64 --backend fmm --reference dense.json --output fmm.json`.
Use `--compact-cases` for large cost runs and `--repeats 2` to check
repeatability. Set the BLAS/OpenMP environment thread counts to match
`--threads` before starting each process.
The compression acceptance budget is at most 0.1% change in heat at the
same mesh; it does not bound the discretization or SIBC model error.

## Why the hole needs its own unknown

The scalar boundary integral equation represents the surface current as
`n × (−grad_s φ)` with a single-valued potential `φ`. Such a current carries
no net current around the hole, so the net circulating current that the
source flux drives through the bore cannot be represented. Filling the bore,
or ignoring the circulation, changes the answer. Radia keeps the hole and
adds one unknown, the net circulating current `α`, represented by the
multivalued potential `α Θ` of a unit-current carrier ring inside the wall.
The extra equation is the loop's work balance in Galerkin form: the electric
work of the surface current and the exterior magnetic work, including the
normal-flux term, are integrated over the whole surface rather than along a
single cut line. The ordinary simply connected BIE is unchanged.

How large the omitted circulation is depends on the geometry and on how much
source flux links the bore. In the public synthetic example below, fixing
`α = 0` overestimated the total loss by about 25%; that figure belongs to the
chosen example and is not a general bound.

## Automatic carrier

The carrier ring is found from the closed surface mesh alone:

- Rays from the z axis locate material intervals at several azimuths and
  heights, including the mid-planes between horizontal shelves of a stepped
  wall. Rays that touch an edge or vertex, run tangent to, or lie in a face
  are treated as ambiguous and retried, never read as a parity change.
- A candidate circle is accepted only if the axis point has winding number 0
  (it lies in the bore), the ring start has winding number ±1 (it lies in the
  material), and every ring segment keeps a positive distance from every
  surface triangle. Together these certify the whole polygon inside the
  material, not just sampled points.
- `solve_loop_extended` returns `carrier_diagnostics`: candidate count,
  rejected rays and candidates with their reasons, and the minimum wall
  clearance of the accepted ring. A failed search raises with the same record.

Supported automatic geometry: one through-hole surrounding the z axis,
including tubes, rings, stepped shafts and non-axisymmetric walls around such
a hole, provided a horizontal circle centred on the z axis fits inside the
material and is found by the candidate search. The certification holds for
the closed flat surface mesh within numerical tolerances, not for the CAD
surface. The unit potential jump, true linear residual and Faraday closure are
checked after construction and solution. Heating and reaction use the full
field, including the circulating component.

## Neumann trace of the carrier field

The carrier enters through its Neumann trace `q_Θ = −H_ring · n`. The trace is
assembled weakly, `b_i = ∫ N_i (−H_ring · n_T) dS` per
triangle, and projected with the surface mass matrix, `q_Θ = M⁻¹ b`, the same
P1 space as the rest of the operator. The former point evaluation with
averaged vertex normals is not a consistent trace at creases and made the loss
depend on where the carrier was placed. In the public example the spread of
the loss over four carriers dropped from 1.34% (point evaluation) to
0.12% (projection). The exact solution of the continuous problem does not
depend on the carrier; the remaining spread is the sensitivity of the
discrete problem (P1 trace and operators, quadrature, path integration). It is
reported, not attributed to a single cause.

## Accuracy checks and the validation example

Every solve checks the true linear residual, the Faraday residual and the
unit jump of `Θ`, and reports two complete powers: the surface loss and the
reaction power of the incident field (including the multivalued part). A large
difference between them flags an under-resolved surface, but a small one does
not bound the error: in the validation example the difference stayed
between 0.055% and 0.070% on all three meshes while the loss differed from
the independent FEM by 1.4%. Read it together with a mesh study and an
independent reference.

`validation_test/induction_heating/bored_workpiece/` compares the loss of the
loop-extended BEM with an independent air-domain HCurl FEM (order 2, SIBC
Robin boundary, no BEM operators) for a synthetic stepped bored workpiece
excited by closed filament loops. On the finest meshes (6,232 BEM vertices;
601,029 FEM unknowns) the two losses differ by 1.40%, within the 2% gate.
Across the three mesh levels of each method the loss varied by 0.21% (BEM)
and 0.010% (FEM) of the value on that method's finest mesh, and the BEM
power-balance difference stayed between 0.055% and 0.070%. Both methods share the
linear SIBC approximation; the example does not validate SIBC against a
resolved conducting volume.

The remaining difference of about 1.4% did not shrink under mesh
refinement of either method and its cause was not isolated. For the total
loss of an induction-heating workpiece this level is treated as agreement,
because for that purpose the material data are expected to dominate the
uncertainty (an engineering judgement, not a quantity measured here). Applications that need
field quality at the percent level, such as accelerator magnets, should not
inherit that tolerance and need their own convergence and reference study.

## Limits

Linear SIBC (uniform, or per-face values on the weak loop route); flat,
undeformed P1 surface triangles (a curved or deformed surface raises; small
fillet radii therefore carry an `O(h²)` geometric error); at most 14,000
surface faces for loop work. Native direct P0 assembly remains quadratic;
the dense-equivalent matrix is about 1 GB at 11,000 faces, and construction
can need several times more. The FMM route retains no NumPy dense P0
matrix; its complete native storage inventory is unavailable and is reported
as unknown, separately from measured process peaks. The default P1 body route retains dense BIE operators and a dense mixed
solve. Explicit `--wp-bem-backend hacapk` uses compressed SL/DL products
and a GMRES operator bordered by the loop and mean-gauge equations. The
mass matrix is factored once with `sparsecholesky`; tagged face impedances
and per-panel ESIM change only the source stiffness. Strong coupling uses
the same complete body operator. Every solve certifies the unscaled true
residual and Faraday row at `1e-6`; failure raises without fallback.

HACApK construction still assembles dense Galerkin entry tables before
compression. Those source tables are released after a successful build,
and the mixed complex matrix is never assembled on the compressed route.
This reduces retained storage, not quadratic construction memory. The
14,000-face and 7,000-vertex guards therefore remain unchanged. Timing and
compression depend on geometry and near-field storage; no general speed
or linear-memory construction claim follows. Record ACA parameters,
H-matrix statistics and GMRES iterations from the returned diagnostics.

Multiple handles, scalar nonlinear ESIM, and geometry with no certified
interior z-axis carrier raise with an explanation. Per-face ESIM impedances
(`--esim-per-panel`) can be passed to the genus-1 loop route; the validation
on this page covers linear SIBC only and says nothing about nonlinear
materials. The solver does not fill the hole, remesh, change material
physics, or discard its cohomology. The field basis order of the thermal FEM
is independent of this surface-BEM restriction.

For an advanced nonstandard geometry, `solve_loop_extended` accepts a
validated `section_anchor` and `carrier_ring`. This bypasses automatic
placement, not the unit-jump or residual checks.

## Meshes and geometry order

Author the CAD in metres (scale the shape before `GenerateMesh`, never the
mesh afterwards), load a checked `.vol` as saved, and do not call
`Mesh.Curve()` on a loaded mesh or while restoring a saved solution.
Cubit-exported curved meshes keep their exported mapping; a higher geometry
order means a new export. See
[Geometry order in Radia workflows](../ngsolve_integration/curve_order.md).
