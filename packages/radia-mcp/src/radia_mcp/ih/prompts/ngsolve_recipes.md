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
Periodic Kelvin `auto` selects SparseCholesky; explicit AMS/BDDC requests fail.
Compound A-V and scalar thermal H1 do not inherit HCurl AMS support.

Before scaling, retrieve `induction_heating("memory")` or
`matrix_solvers_direct("capacity")` for the diagnostic report checklist.
Direct can be fast when the factor fits; factor fill-in and Windows commit
headroom govern capacity, not element count alone. Record outer solver,
preconditioner and coarse solver: generic BDDC can retain direct factorization.
This IH BDDC route explicitly configures AMS coarse cycles. Aggregate t_solve_s
is not isolated factor time; absent memory/factor telemetry stays unavailable.
No universal speed claim or new validation result follows from this recipe.

### Gauge fixing for a port-driven A–V coil with a SIBC workpiece

For the FEM gauge reference, see `fem_gauge_open_boundary(topic="coil_av_tree")`.

This applies to a meshed conducting coil with source and sink faces,
surrounding air, and a workpiece represented by a surface impedance. See
`radia.panels.calc_fem_coilmesh` for the solver in the installed Radia
version. The spanning-tree method below describes the reviewed solver
candidate. Check which gauge the installed solver uses before reproducing
its results.

With the solver's scaled scalar potential `phi` and `s = 1j * omega`,
the coil electric field is `E = -s * (A + grad(phi))`. The transformation
`A -> A + grad(chi)`, `phi -> phi - chi` preserves this field and
`curl(A)`. Lowest-order edge elements still contain vertex gradients with
`nograds=True`. A point constraint on an H1 potential or a choice of
sparse solver cannot remove this HCurl nullspace.

- **SIBC trace.** The surface law uses the tangential trace of `A`.
  A gauge function must be constant on each connected SIBC component so
  its tangential gradient vanishes there. Collapse each component to its
  own node in the gauge graph, accounting for periodic vertex
  identifications. Keep the physical trace degrees of freedom.
- **Ports.** Keep the source and sink as separate graph nodes. Since
  `phi` is prescribed on the ports, the homogeneous gauge must preserve
  those data. The candidate also fixes the line integral of `A` along a
  designated physical source-to-sink path. This defines the voltage
  convention across the gap. Choose the path on actual mesh edges,
  accounting for periodic copies and avoiding the SIBC surface. Reject
  the setup if no such path exists. A unit `phi` lift corresponds to a
  voltage of magnitude `abs(s)` in this convention.
- **Tree and verification.** With the path held fixed, use a spanning
  tree of the resulting graph to constrain the lowest-order edge
  degrees of freedom. Check that different admissible trees have
  different masks and reproduce current-normalized inductance and both
  losses. Check the true residual and the unnormalized terminal current.
  A small residual can coexist with an incorrect port drive. When
  comparing formulations, keep the mesh, excitation and complex
  `Z_s` fixed first, then check against an independent reference with
  compatible ports and return path.

An air-coil inductance comparison checks the port convention in that
setup. Coil losses still require mesh refinement against the copper
skin depth. When a surface-current reference reports external
inductance, account for the conductor's internal contribution before
comparing it with total volume-FEM inductance.

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
