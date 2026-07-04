# 2D planar machine modelling — MMMM / HDiv-VIM + shared eddy/PM coupling

Executable showcase for the per-unit-length (motor cross-section) stack in `radia`:

- **`radia.mmmm2d`** — collocation MMMM soft-iron demag (C++ 2D log-kernel moment core).
- **`radia.vim._vim2d`** — HDiv-VIM (RT1 charge Gram) soft-iron demag.
- **`radia.planar_charges`** — shared exterior field / A_z / complex field / Maxwell torque + force /
  permanent-magnet field (method-agnostic: either demag method's per-element M feeds it).
- **`radia.planar_eddy`** — shared staggered eddy-current coupling (maglev / induction machine /
  eddy-current brake): analytic iron demag ↔ NGSolve reduced-potential complex `A_z` eddy FEM.

## Files

| File | Role |
|------|------|
| `planar_coupling_showcase.ipynb` | **Executed** notebook (outputs + figures embedded): demag factors, reluctance-torque sweep, eddy vs analytic Bessel, PM rotor (design B), unified PM+iron+eddy rotor. |
| `planar_coupling_results.json` | Synchronized result sidecar (values + `generated_at_utc` + version/runtime metadata + `notebook_sha256`). |
| `build_showcase.py` | Reproducible builder: constructs + executes the notebook, then finalizes the JSON sidecar. Re-run after editing; commit the `.ipynb` + `.json` pair together. |

```bash
python docs/planar_coupling/build_showcase.py
```

## Validation

Every figure/number is gated against a closed form or a monolithic FEM (goldens under
`validation_test/feec/`): Clausius–Mossotti demag factors, ellipse `Dx=b/(a+b)`, conducting-cylinder
Bessel `⟨Bx⟩/B0 = 2 I₁(z)/(z I₀(z))`, 2D dipole `a²M/(2r²)`, and the staggered coupling vs a
monolithic AC(+PM) FEM. Queryable knowledge: MCP `motor_planar_coupling`.
