# GMSH Animation

Displacement animation using GMSH NodeData views with STEP geometry overlay.

The saved `gmsh_animation.ipynb` is a historical calculation, not a current
format acceptance example. Its original assets are not tracked; the archived
record retains its original calculation date, Radia version and format under
[`validation_test/gmsh_animation`](../../validation_test/gmsh_animation/).
It is marked stale and pins the unchanged notebook bytes. Do not infer today's
output contract from those saved results. Current exports accept MSH 4.2 only.
Regenerating this example awaits recovery of the original rotor assets.

`moving_body_demo.py` is the notebook-coupled generator/helper that creates the
rotor `.vol`, GMSH `.msh`, and `.geo` companion.

`gmsh_animation_export.ipynb` is the result-saved export lesson. It opens the
existing `animation.geo`, confirms the `animation.geo.opt` / `.msh.opt`
sidecar roles, synchronizes Gmsh time steps, exports PNG frames, and writes
`gmsh_animation_export.gif` plus `gmsh_animation_export.mp4`. The executed
notebook owns the public result;
[`gmsh_animation_export_results.json`](../../validation_test/gmsh_animation/gmsh_animation_export_results.json)
is the checked validation record outside the docs tree.
