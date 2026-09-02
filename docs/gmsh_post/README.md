# gmsh post-processing showcase

Canonical, executed examples of the standard EM post-processing lane:
`GmshPostExport` output interrogated and rendered through the
`radia-mcp` gmsh verbs (`flux_lines`, `streamlines_2d`, `streamlines`,
adaptive isosurfaces, profiles, renders).

| File | Role |
|---|---|
| `em_fieldlines.ipynb` | Executed showcase: Case A = circular coil (exact equal-flux field lines from psi = r A_theta, cross-checked against the evenly spaced tracer); Case B = opposed-PM gap (mid-plane streamlines, gap profile, nested adaptive \|B\| isosurfaces with cut-away, orbit GIF). |
| `em_fieldlines_results.json` | Optional historical result record; the executed notebook owns the public demonstration. |
| `em_post_gallery.ipynb` | Executed gallery of the 2026-08 verb generations on ONE problem -- a `saddle_coil` dipole with and without an HDiv-VIM soft-iron flux return: figure controls (camera presets, pinned colour range, labelled axes), `render_panels` shared-scale air-vs-iron comparison, `select` compound queries with carried values, TRUE ray-cast volume + LIC with the coil STEP depth-composited, and `time_series` over a moving-magnet file series. This is also the standard-post-set-applied-to-iron example the earlier handoff deferred until `radia.vim` settled. |
| `em_post_gallery_results.json` | Optional historical result record for the gallery. |
| `em_particle_orbits.ipynb` | Executed showcase of CHARGED-PARTICLE ORBITS (`track_lorentz_ivp` + `export_particle_tracks_msh`, and the `gmsh_particle_trace` verb): a `saddle_coil` dipole sorting 15/20/25 MeV electrons into a dispersion fan, the same tracks as a flying-beam GIF, a permanent-magnet quadrupole focusing in one plane within 1.1% of the thick-lens focal length, and edge focusing on the shipped Maxwellian tilted fringe -- where flying both charge signs splits the measurement into an edge term matching the hard-edge law to 0.07% and a fringe term matching the Enge form to 0.01%. |
| `em_particle_orbits_results.json` | Optional historical result record for the orbits notebook. |
| `output/` | Regenerated artifacts (`.msh`, `.pos`, STEP, PNG/GIF) -- not tracked; re-running the notebooks rebuilds them. |

The recipe reference lives in the `radia-mcp` gmsh knowledge:
`gmsh_usage(topic="paraview")` (tool selection matrix, isosurface
beautification levers, measured pitfalls) and
`gmsh_usage(topic="policy")` (axis-equal figure policy -- a beamline
envelope is the documented exception, and the orbits notebook states
its magnification).

The same battery applies verbatim to any solver output that goes
through `GmshPostExport` (HDiv-VIM, panels, PEEC post): export the
fields, then point the verbs at the `.msh`.
