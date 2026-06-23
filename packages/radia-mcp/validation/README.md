# radia-mcp validation

This directory contains long-running, physics-oriented validation scripts.
They are intentionally outside `tests/` so the default pytest gate stays
focused on fast MCP/API contracts and small numerical invariants.

Run a validation script directly:

```bash
python validation/age/validate_age_nonlinear_saturation.py
python validation/comsol_ngsolve/validate_carter.py
python validation/force/validate_force_xval.py
```

Validation scripts may run real FEM solves, compare against brute-force
references, or reproduce research examples. They should be executable and
self-contained, but they are not part of the default pytest collection.
