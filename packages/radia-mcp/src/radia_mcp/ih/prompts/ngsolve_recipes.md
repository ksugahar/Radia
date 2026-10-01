# NGSolve Implementation Recipes for IH SIBC

Production-grade NGSolve code patterns extracted from the lab's
`calc_inductance.py`, `calc_fem_kelvin.py`, `bem_sibc_solver.py`,
and the supporting `radia.sparsesolv_ngsolve` (Compact AMS / COCR).

Use these as starting points when extending the IH pipeline. All
code follows the lab's policies:
- Verify-First Policy (FES check before physics solve)
- Compact HX preconditioner (HYPRE-free, TaskManager-native)
- Shifted preconditioner for air + conductor problems
- No fallback chains (fail-fast)

## Recipe 1: Uniform scalar BIE + SIBC

Production reference: `radia.bem_sibc_solver.ScalarBIESIBCSolver`.
The implemented surface equation is
`(0.5 M - DL + gamma SL M^-1 K) phi = rhs`, with a mean-potential
Lagrange multiplier and uniform `gamma = Z_s / (1j omega mu0)`.
Use the assembled solver's `solve(phi_inc, Z_s, omega)` for dense solving,
`solve_hacapk(...)` for its assembled HACApK operator, or
`solve_iterative(...)` for its matrix-free operator. Iterative paths use GMRES,
not HCurl AMS or COCR. Caller owns NGSolve TaskManager where required.

All three paths recompute the augmented true residual against the fixed load
norm and require the shared `RELATIVE_LIMIT = 1e-6`; GMRES nonconvergence raises.
Returned `linear_residual_rel` and `linear_residual_limit` expose this check.

## Recipe 2: Variable impedance is not implemented in scalar BEM

A scalar or constant nodal array is accepted. A nonuniform `Z_s` array raises
`NotImplementedError`. Multiplying output rows by impedance is not equivalent
to source-side weighted surface stiffness `K_gamma` followed by `SL M^-1`.
Correct variable-coefficient assembly and local loss integration still need
implementation and independent validation. Curvature or ESIM cell values must
not be fed into the removed output-row scaling route.

## Recipe 3: Use the supported FEM-SIBC entry point

Run `python -m radia.panels.calc_fem_kelvin --help` for the actual model and
source arguments. The implementation owns the geometry, Kelvin mapping,
source integration and residual checks; do not reconstruct them from the
removed schematic snippet with undefined `nu_kelvin_cf`/`A`/`prec` objects.

For supported nonperiodic HCurl spaces, `--solver auto` uses AMS at p=1 and
BDDC+AMS at p=2/3. The latter explicitly selects an edge-only wirebasket and
three AMS coarse cycles. Periodic Kelvin and compound-space restrictions
remain distinct; inspect the current dispatch errors rather than assume a
preconditioner applies to every space. `sparsecholesky` is the direct option.
No universal speed claim or new validation result follows from this recipe.

## Recipe 4: Verify-First Policy on a typical IH FES setup

**MANDATORY** before any 5-minute physics solve. Each check is <1 s.

```python
from ngsolve import (H1, HCurl, Periodic, GridFunction, Integrate,
                     CoefficientFunction as CF)

mesh = Mesh("model.vol")

# Check 1: materials and boundaries spelled as expected
print("Materials:", mesh.GetMaterials())   # expect: air, kelvin (NOT wp!)
print("Boundaries:", mesh.GetBoundaries())  # expect: sibc, gnd,
                                            #         kelvin_int, kelvin_ext

# Check 2: Periodic actually constrains slave DOFs
fes_base = HCurl(mesh, order=1, dirichlet="gnd", complex=True)
fes = Periodic(fes_base)
slaved = sum(fes_base.FreeDofs()) - sum(fes.FreeDofs())
print(f"Kelvin slaved DOFs: {slaved} (should be > 0)")
assert slaved > 0

# Check 3: Functional Kelvin pair test
gfu = GridFunction(fes); gfu.vec[:] = 0
gfu.Set(CF((1.0, 0, 0)), definedon=mesh.Boundaries("kelvin_int"))
norm_int = Integrate(InnerProduct(gfu, gfu), mesh,
                      definedon=mesh.Boundaries("kelvin_int"))
norm_ext = Integrate(InnerProduct(gfu, gfu), mesh,
                      definedon=mesh.Boundaries("kelvin_ext"))
ratio = norm_ext / max(norm_int, 1e-30)
print(f"Kelvin int->ext ratio: {ratio:.6f} (should be 1.0)")
assert abs(ratio - 1.0) < 1e-3

# Check 4: SIBC boundary has nonzero area
A_sibc = Integrate(CF(1), mesh, BND,
                   definedon=mesh.Boundaries("sibc"))
print(f"SIBC area: {A_sibc:.6e} m^2 (should match coil-facing wp surf)")
assert A_sibc > 1e-9

# Check 5: GND vertex is identified
gnd_dofs = [i for i in range(fes_base.ndof) if not fes_base.FreeDofs()[i]]
print(f"GND-constrained DOFs: {len(gnd_dofs)} (should be > 0)")
assert len(gnd_dofs) > 0
```

If any check fails, FIX THE GEOMETRY/LABELS first. Do NOT iterate
solver runs to discover an FES-level bug.

Real-world cost saving (lab 2026-04-25):
- Wrong-debug path: 10 minutes physics solve x N iterations
- Right path: 5 FES checks at <1 s each + 1 physics solve

## Cross-MCP recipe references

- `radia_mcp.matrix_solvers` -- Krylov solver selection (COCR vs BiCGStab)
- `radia_mcp.bem` -- ngsolve.bem operator catalog (LaplaceDL/SL etc.)
- `radia_mcp.fem` -- FEM formulation theory (A-Omega, T-Omega, gauging)
- `radia_mcp.radia_ngsolve.analytical_formulas('validation_use_cases')`
  -- closed-form reference for each FEM/BEM analysis
