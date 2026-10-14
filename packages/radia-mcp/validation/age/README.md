# AGE validation scripts

These scripts validate the analytic-gap element (AGE) rotating-machine
workflow against brute meshed-gap solves, analytic motor formulas, or
longer physics sweeps. They were moved out of `tests/` because several
cases take tens of seconds and are better treated as explicit validation.

Run one case directly from `packages/radia-mcp`:

```bash
python validation/age/validate_age_pmsm_physical.py
python validation/age/validate_age_nonlinear_saturation.py
```

Keep new pytest tests small and fast. If a check needs a real FEM sweep,
brute-reference solve, or publication-style physics validation, add it
here instead of under `tests/`.
