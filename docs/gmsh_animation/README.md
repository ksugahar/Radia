# GMSH Animation

Displacement animation using GMSH NodeData views with STEP geometry overlay.

The result-saved notebook is `gmsh_animation.ipynb`. It inspects the docs-local
`rotor_animation.msh` artifact and verifies the GMSH v4.1 headers and 21 vector
`$NodeData` frames. Its public result is embedded in the executed notebook;
checked evidence lives in
[`validation_test/gmsh_animation`](../../validation_test/gmsh_animation/).

`moving_body_demo.py` is the notebook-coupled generator/helper that creates the
rotor `.vol`, GMSH `.msh`, and `.geo` companion.

`gmsh_animation_export.ipynb` is the result-saved export lesson. It opens the
existing `animation.geo`, confirms the `animation.geo.opt` / `.msh.opt`
sidecar roles, synchronizes Gmsh time steps, exports PNG frames, and writes
`gmsh_animation_export.gif` plus `gmsh_animation_export.mp4`. The executed
notebook owns the public result;
[`gmsh_animation_export_results.json`](../../validation_test/gmsh_animation/gmsh_animation_export_results.json)
is the checked validation record outside the docs tree.
