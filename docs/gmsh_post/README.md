# gmsh post-processing showcase

Canonical, executed examples of the standard EM post-processing lane:
`GmshPostExport` output interrogated and rendered through the
`radia-mcp` gmsh verbs (`flux_lines`, `streamlines_2d`, `streamlines`,
adaptive isosurfaces, profiles, renders).

| File | Role |
|---|---|
| `em_fieldlines.ipynb` | Executed showcase: Case A = circular coil (exact equal-flux field lines from psi = r A_theta, cross-checked against the evenly spaced tracer); Case B = opposed-PM gap (mid-plane streamlines, gap profile, nested adaptive \|B\| isosurfaces with cut-away, orbit GIF). |
| `em_fieldlines_results.json` | Synchronized result sidecar (key numbers, versions, `notebook_sha256`). |
| `em_post_gallery.ipynb` | Executed gallery of the 2026-08 verb generations on ONE problem -- a `saddle_coil` dipole with and without an HDiv-VIM soft-iron flux return: figure controls (camera presets, pinned colour range, labelled axes), `render_panels` shared-scale air-vs-iron comparison, `select` compound queries with carried values, TRUE ray-cast volume + LIC with the coil STEP depth-composited, and `time_series` over a moving-magnet file series. This is also the standard-post-set-applied-to-iron example the earlier handoff deferred until `radia.vim` settled. |
| `em_post_gallery_results.json` | Synchronized result sidecar for the gallery. |
| `output/` | Regenerated artifacts (`.msh`, `.pos`, STEP, PNG/GIF) -- not tracked; re-running the notebooks rebuilds them. |

The recipe reference lives in the `radia-mcp` gmsh knowledge:
`gmsh_usage(topic="paraview")` (tool selection matrix, isosurface
beautification levers, measured pitfalls) and
`gmsh_usage(topic="policy")` (axis-equal figure policy).

The same battery applies verbatim to any solver output that goes
through `GmshPostExport` (HDiv-VIM, panels, PEEC post): export the
fields, then point the verbs at the `.msh`.
