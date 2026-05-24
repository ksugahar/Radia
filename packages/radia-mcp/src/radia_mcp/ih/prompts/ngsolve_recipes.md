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

## Recipe 1: Scalar BIE + SIBC Robin term with COCR + ComplexCompactAMS

Production reference: `src/radia/bem_sibc_solver.py::ScalarBIESIBCSolver`.

```python
from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction,
                     grad, ds, dx, InnerProduct, specialcf)
from ngsolve.bem import LaplaceDL, LaplaceSL, SingleLayerPotential
import radia.sparsesolv_ngsolve as ssn

# Surface H1 space (BND only)
fes = H1(mesh, order=1, definedon=mesh.Boundaries("wp_bnd"))
u, v = fes.TnT()

# System: (1/2 M - DL + gamma * SL * M^-1 * K) phi = M phi_inc
# gamma = Z_s / (jw * mu0)  -- complex
gamma = Zs / (1j * omega * mu0)

M = BilinearForm(fes); M += u*v*ds; M.Assemble()
K = BilinearForm(fes); K += grad(u).Trace() * grad(v).Trace() * ds; K.Assemble()
DL = LaplaceDL(fes, fes).Assemble()
SL = LaplaceSL(fes, fes).Assemble()

# Build A = 0.5 * M - DL + gamma * SL * M^-1 * K (Lagrange mult for gauge)
# ... assembly details in bem_sibc_solver.py

# COCR solver (Sogabe-Zhang 2007) for complex symmetric system
solver = ssn.COCRSolver(
    matrix=A, preconditioner=prec,
    tol=1e-7, maxiter=400
)
phi.vec.data = solver * f.vec
```

**When to use**: BEM Scalar BIE (no FEM volume), nonlinear ESIM
workpieces via Karl outer iteration. Reference: lab 2026-04-13
Smythe sphere validation (-1.6% accuracy).

## Recipe 2: Per-node Z_s array via per-panel curvature extraction

Production reference: `src/radia/panels/calc_inductance.py::
_compute_panel_local_radii` and `_build_per_node_Zs`.

```python
import numpy as np
from ngsolve import H1, GridFunction, SurfaceL2

def per_node_curvature_h1(mesh, surface_label, percentile=10):
    # Discrete normal-angle radius extractor on a workpiece surface.
    # Returns one R_local per H1 vertex on `surface_label`. Used to
    # compute per-node Z_s for the SIBC Robin block.
    fes_p1 = H1(mesh, order=1, definedon=mesh.Boundaries(surface_label))
    ndof = fes_p1.ndof
    R_local = np.zeros(ndof)

    # Build panel adjacency from BND elements
    panels = list(mesh.Boundaries(surface_label).Elements())
    # ... per-panel normal, centroid, neighbor pairs (see calc_inductance.py)

    # For each panel: angle_ij = arccos(n_i . n_j), dist_ij = |c_j - c_i|
    # R_ij = dist_ij / angle_ij ; R_panel[i] = percentile(R_ij, 10)
    R_panel = compute_per_panel_R(panels, percentile=percentile)

    # Adaptive sliver clamp: R_floor[i] = 0.5 * panel_diameter[i]
    R_panel = np.maximum(R_panel, 0.5 * panel_diameters)

    # Vertex averaging to project panel -> H1 nodes
    for vid in range(ndof):
        neighbor_panels = vertex_to_panels[vid]
        R_local[vid] = np.mean([R_panel[p] for p in neighbor_panels])

    return R_local

# Use per-node Z_s in ScalarBIESIBCSolver
R_local = per_node_curvature_h1(mesh, "wp_bnd")
Zs_per_node = (1 + 1j) / (sigma * skin_depth(omega, mu0*mur, sigma)) * \
              (1 + (1 + 1j)/2 * skin_depth(omega, mu0*mur, sigma) / R_local)
solver.solve(phi_inc, Z_s=Zs_per_node, omega=omega)   # ndarray accepted
```

**Validation** (lab 2026-04-12, prolate spheroid `a/b=4`):
- Cu 10 Hz, scalar(mean R): -11.2% error
- Cu 10 Hz, per-node(mesh): +3.1% error (4x improvement)

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
