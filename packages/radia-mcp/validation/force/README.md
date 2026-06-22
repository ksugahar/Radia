# Force validation scripts

This directory contains the heavy electromagnetic force/energy
cross-validation suite for `radia_mcp.radia_ngsolve.force`.

Run explicitly from `packages/radia-mcp`:

```bash
python validation/force/validate_force_xval.py
```

The full force suite may take minutes because it runs real 3D FEM solves.
