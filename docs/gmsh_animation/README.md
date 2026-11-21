# GMSH Animation

Displacement animation using GMSH NodeData views with STEP geometry overlay.

`gmsh_animation.ipynb` was regenerated on 2026-11-21 with Radia 5.4.3,
NGSolve 6.2.2608 and the Gmsh 5 development build: the rotor `.msh` is MSH 4.2
with 21 displacement `$NodeData` frames. The checked records live under
[`validation_test/gmsh_animation`](../../validation_test/gmsh_animation/).
Exports accept MSH 4.2 only.

`moving_body_demo.py` is the notebook-coupled generator/helper that creates the
rotor `.vol`, the stator STEP and the GMSH `.msh`; the tracked `animation.geo`
(+ `.geo.opt`) merges them and is not rewritten.

`gmsh_animation_export.ipynb` is the result-saved export lesson. It opens the
existing `animation.geo`, confirms the `animation.geo.opt` / `.msh.opt`
sidecar roles, synchronizes Gmsh time steps, exports PNG frames, and writes
`gmsh_animation_export.gif` plus `gmsh_animation_export.mp4`. The executed
notebook owns the public result;
[`gmsh_animation_export_results.json`](../../validation_test/gmsh_animation/gmsh_animation_export_results.json)
is the checked validation record outside the docs tree.
