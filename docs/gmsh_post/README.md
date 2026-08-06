# gmsh post-processing showcase

Canonical, executed examples of the standard EM post-processing lane:
`GmshPostExport` output interrogated and rendered through the
`radia-mcp` gmsh verbs (`flux_lines`, `streamlines_2d`, `streamlines`,
adaptive isosurfaces, profiles, renders).

| File | Role |
|---|---|
| `em_fieldlines.ipynb` | Executed showcase: Case A = circular coil (exact equal-flux field lines from psi = r A_theta, cross-checked against the evenly spaced tracer); Case B = opposed-PM gap (mid-plane streamlines, gap profile, nested adaptive \|B\| isosurfaces with cut-away, orbit GIF). |
| `em_fieldlines_results.json` | Synchronized result sidecar (key numbers, versions, `notebook_sha256`). |
| `output/` | Regenerated artifacts (`.msh`, `.pos`, PNG/GIF) -- not tracked; re-running the notebook rebuilds them. |

The recipe reference lives in the `radia-mcp` gmsh knowledge:
`gmsh_usage(topic="paraview")` (tool selection matrix, isosurface
beautification levers, measured pitfalls) and
`gmsh_usage(topic="policy")` (axis-equal figure policy).

The same battery applies verbatim to any solver output that goes
through `GmshPostExport` (HDiv-VIM, panels, PEEC post): export the
fields, then point the verbs at the `.msh`.
