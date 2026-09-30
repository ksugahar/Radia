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

## Recipe 3: FEM + Kelvin + Robin BC (hole approach)

Production reference: `src/radia/panels/calc_fem_kelvin.py`.

**KEY**: workpiece is SUBTRACTED from mesh (hole approach), NOT meshed.
SIBC = Robin BC on the hole boundary. Avoids the -34% systematic error
of the interface approach.

```python
from ngsolve import (Mesh, HCurl, Periodic, BilinearForm, LinearForm,
                     GridFunction, curl, dx, ds, InnerProduct,
                     specialcf, x, y, z, Preconditioner, BVP)
import radia.sparsesolv_ngsolve as ssn

mesh = Mesh("model.vol")  # hole approach + Kelvin pair

# Verify-First Policy: check FES and Kelvin pair BEFORE solving
print("Materials:", mesh.GetMaterials())
print("Boundaries:", mesh.GetBoundaries())

fes_base = HCurl(mesh, order=1, dirichlet="gnd",
                  complex=True, gradientdomains={})
fes = Periodic(fes_base)  # Kelvin pair = master/slave
slaved = sum(fes_base.FreeDofs()) - sum(fes.FreeDofs())
assert slaved > 0, "Kelvin Periodic identification failed"

# Functional test: Set 1 on kelvin_int -> ratio on kelvin_ext should be 1.0
gfu_test = GridFunction(fes); gfu_test.vec[:] = 0
gfu_test.Set(1.0, definedon=mesh.Boundaries("kelvin_int"))
ratio = Integrate(gfu_test*gfu_test, mesh,
                   definedon=mesh.Boundaries("kelvin_ext")) / \
        Integrate(gfu_test*gfu_test, mesh,
                   definedon=mesh.Boundaries("kelvin_int"))
assert abs(ratio - 1.0) < 1e-3, f"Kelvin ratio = {ratio}, expected 1.0"

# Build system
u, v = fes.TnT()
nu = 1.0 / (mu0 * mur_air)
nu_kelvin = nu_kelvin_cf(mesh)  # = nu0 * (r'/R)^2

a = BilinearForm(fes, symmetric=False)
a += nu * curl(u) * curl(v) * dx("air")
a += nu_kelvin * curl(u) * curl(v) * dx("kelvin")
a += (1j * omega / Zs) * InnerProduct(u.Trace(), v.Trace()) * \
     ds(definedon=mesh.Boundaries("sibc"))  # Robin SIBC

# Shifted preconditioner for air + conductor (eps = 1e-6 * nu)
a_shifted = BilinearForm(fes, symmetric=False)
a_shifted += nu * curl(u) * curl(v) * dx
a_shifted += 1e-6 * nu * u * v * dx   # eps * mass shift on prec only

c = Preconditioner(a_shifted, "compactams")  # Compact HX (HYPRE-free)

# RHS: Biot-Savart from PEEC filaments
f = LinearForm(fes)
f += InnerProduct(H_inc_cf, v.Trace()) * \
     ds(definedon=mesh.Boundaries("sibc"))

# Solve via BiCGStab + Compact HX
A.Assemble(); c.Update(); f.Assemble()
inv = ssn.BiCGStab(matrix=A, c=c, tol=1e-7, maxiter=200)
gfu.vec.data = inv * f.vec
```

**Validated**: 3D FEM-Kelvin L = 90.71 nH vs analytical torus 88.5 nH
(+2.5% on coarse mesh, lab 2026-04-12). 2D axisym: L < 1%, P < 2%.

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
