# Gmsh animation validation

This directory owns the checked machine-readable evidence behind the executed
animation notebooks in `docs/gmsh_animation/`. The notebooks retain their
display outputs; the helper CLIs regenerate these JSON records here without
creating docs-local result sidecars.

```powershell
python docs/gmsh_animation/moving_body_demo.py
python docs/gmsh_animation/gmsh_animation_inspect.py
python docs/gmsh_animation/gmsh_animation_export.py
python -m pytest validation_test/gmsh_animation/test_gmsh_animation_evidence.py -q
```

The `.geo.opt` and `.msh.opt` files remain beside the launch artifacts because
they are Gmsh execution settings, not result bookkeeping.
