# COMSOL/NGSolve validation scripts

These scripts hold the older cross-validation cases that compare closed-form
CAE formulas against NGSolve FE anchors. They are intentionally outside
`tests/` because they can require optional FE dependencies and real solves.

Run one case directly from `packages/radia-mcp`:

```bash
python validation/multiphysics_ngsolve/validate_carter.py
python validation/multiphysics_ngsolve/validate_waveguide_dispersion.py
```

Keep `tests/` for fast closed-form contracts. Put FE anchors and
tool-to-tool cross-checks here.
