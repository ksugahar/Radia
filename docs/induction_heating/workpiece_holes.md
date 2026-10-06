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
the circulating current. The log reports the topology and backend selection.

Supported automatic geometry: one through-hole surrounding the z axis,
including tubes, rings and stepped shafts. The carrier is located from the
closed surface mesh, and its vertices and edge midpoints are checked against
material intervals. The potential jump, true linear residual and Faraday
closure are checked after construction/solution. Heating and reaction use
the full field, including the circulating component.

Limits: linear, spatially uniform SIBC; P1 surface basis; at most 7000 surface
vertices for the current dense loop operator. Multiple handles, nonlinear
ESIM, and geometry with no verified interior z-axis carrier raise with an
explanation. The solver does not fill the hole, remesh, change material
physics, or discard its cohomology. The field basis order of the thermal FEM
is independent of this surface-BEM restriction.

For an advanced nonstandard geometry, `solve_loop_extended` accepts a
validated `section_anchor` and `carrier_ring`. This bypasses automatic
placement, not the unit-jump or residual checks.

Load a checked `.vol` as saved; do not call `Mesh.Curve()` while restoring a
saved solution. Cubit-exported curved meshes retain their exported mapping.
