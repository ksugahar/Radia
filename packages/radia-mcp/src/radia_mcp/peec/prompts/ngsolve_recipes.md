# NGSolve Implementation Recipes for HOIBC

Production-grade NGSolve code patterns for the 3-Laplace-BVP cascade
(Dong-Di Rienzo 2020). Goes BEYOND the existing
`calc_fem_kelvin.py --formulation total` 1st-order Robin BC.

Status: design recipes for the upcoming `--ibc-order {0,1,2}` flag.
The 1st-order path is already production; orders 0 and 2 are new.

## Recipe 1: BVP_0 (PEC approximation)

```python
from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction,
                     grad, dx, ds, InnerProduct, specialcf, x, y, z,
                     Preconditioner)

mesh = Mesh("model.vol")  # workpiece is meshed as a HOLE (subtracted)

# Reduced magnetic scalar potential phi: H_ext = H_s - grad(phi)
fes = H1(mesh, order=p, dirichlet="kelvin_outer", complex=True)

phi_0 = GridFunction(fes)
u, v = fes.TnT()

# LHS: Laplace
a = BilinearForm(fes, symmetric=True)
a += grad(u) * grad(v) * dx
c = Preconditioner(a, "bddc")  # register before assembly
a.Assemble()

# RHS: Neumann data on conductor surface = -H_s.n
n = specialcf.normal(3)
H_s_cf = biot_savart_H_cf_from_peec(peec_topology, current=I_coil)

f = LinearForm(fes)
f += -InnerProduct(H_s_cf, n) * v * ds(
    definedon=mesh.Boundaries("conductor_bnd"))
f.Assemble()

# Real Laplace operator with a Dirichlet outer boundary: CG + BDDC
from ngsolve.krylovspace import CGSolver
inv = CGSolver(a.mat, c.mat, tol=1e-10, maxiter=1000)
phi_0.vec.data = inv * f.vec

H_ext_PEC = H_s_cf - grad(phi_0)  # PEC result
```

## Recipe 2: BVP_1 (Leontovich correction RHS)

```python
# BVP_1: same LHS as BVP_0 (Laplace), different Neumann RHS
# RHS_1 = -(1/sqrt(s)) * div_s ((H_s - grad phi_0)_tangential)
# where div_s = surface divergence (tangential gradient component)

phi_1 = GridFunction(fes)

# Use specialcf.Weingarten for tangential derivatives (NGSolve 6.2.2603)
# Or compute surface gradient via projection onto tangent plane

# Compute H_t at conductor surface (tangential of PEC field)
H_t_PEC = (H_s_cf - grad(phi_0)) - \
          InnerProduct((H_s_cf - grad(phi_0)), n) * n  # tangential proj

# Surface divergence of H_t (use H1 boundary FES + integration by parts)
# RHS_1 = -(1/sqrt(j*omega)) * div_s(H_t_PEC) . v
f1 = LinearForm(fes)
f1 += -1.0/sqrt(1j*omega) * InnerProduct(grad(v).Trace(), H_t_PEC) * \
      ds(definedon=mesh.Boundaries("conductor_bnd"))
f1.Assemble()

phi_1.vec.data = inv * f1.vec   # reuse a.mat inverse (same Laplace)

# Leontovich result: H_ext_1 = H_s - grad(phi_0 + p_tilde * phi_1)
p_tilde = mur * delta(omega, sigma, mur) / D_char  # small parameter
H_ext_Leontovich = H_s_cf - grad(phi_0) - p_tilde * grad(phi_1)
```

## Recipe 3: BVP_2 (Mitzner curvature correction)

The headline HOIBC advance — captures surface curvature explicitly.

```python
# BVP_2: same LHS (Laplace), Neumann RHS includes curvature term
# RHS_2 = -(1/sqrt(s)) * div_s(-grad phi_1)
#       - (1/(mur * sqrt(s))) * 2*H_mean * H_t_PEC_tangential

phi_2 = GridFunction(fes)

# Surface mean curvature: 2 * H_mean = 1/d_1 + 1/d_2
# NGSolve specialcf for Weingarten map (3D)
W = specialcf.Weingarten(3)    # 2x2 Weingarten map on boundary
H_mean = 0.5 * Trace(W)        # mean curvature in 1/m

# Term 1: tangential derivative of (-grad phi_1)
grad_phi1_t = grad(phi_1) - InnerProduct(grad(phi_1), n) * n
term1 = 1.0/sqrt(1j*omega) * InnerProduct(grad(v).Trace(), grad_phi1_t)

# Term 2: curvature correction
term2 = -1.0/(mur * sqrt(1j*omega)) * 2 * H_mean * \
        InnerProduct(H_t_PEC, v.Trace())

f2 = LinearForm(fes)
f2 += (term1 + term2) * ds(definedon=mesh.Boundaries("conductor_bnd"))
f2.Assemble()

phi_2.vec.data = inv * f2.vec   # reuse a.mat inverse

# Mitzner result: H_ext_2 = H_s - grad(phi_0 + p_tilde phi_1 + p_tilde^2 phi_2)
H_ext_Mitzner = H_s_cf - grad(phi_0) - p_tilde * grad(phi_1) - \
                p_tilde**2 * grad(phi_2)
```

## Recipe 4: Multi-frequency post-process (zero-cost frequency sweep)

The key computational win: BVPs 0/1/2 are **frequency-independent**
Laplace problems. Once solved, evaluate at any frequency by plugging
in `p_tilde(omega)`:

```python
# Pre-solve (frequency-independent)
solve_bvp_0()  # phi_0 -- 33 s
solve_bvp_1()  # phi_1 -- +52 s
solve_bvp_2()  # phi_2 -- +55 s
# Total: ~2 min 20 s ONE TIME for all frequencies

# Frequency sweep is O(1) per point
def H_ext_at_freq(omega, mur, sigma, D_char):
    delta = sqrt(2/(omega * mu0 * mur * sigma))
    p_tilde = mur * delta / D_char
    return H_s_cf - grad(phi_0) - \
           p_tilde * grad(phi_1) - \
           p_tilde**2 * grad(phi_2)

# 100-point frequency sweep: 100 * O(1) = milliseconds
for f in logspace(2, 6, 100):
    omega = 2*pi*f
    H_ext = H_ext_at_freq(omega, mur=100, sigma=2e6, D_char=0.01)
    # Compute L, P, etc. from H_ext via boundary integrals
```

**Compare to existing pipeline**: standard FEM-Kelvin would re-solve
the volume eddy-current problem per frequency (~2 min each x 100 =
3.3 hours). HOIBC: 2 min one-time + 100 * milliseconds = ~2 min total.

## Recipe 5: Validation against Dong-Di Rienzo 2020 prolate ellipsoid

```python
# Geometry: prolate ellipsoid x^2/b^2 + y^2/b^2 + z^2/a^2 = 1
#           a = 20 mm, b = 10 mm (a/b = 2)
#           Surrounded by circular coil at I = 10 A

# Lab equivalent geometry script: prolate_spheroid_hoibc.jou (planned)
# validation_test/mixed_galerkin/sphere/02_hoibc_gamma1.py

# Validation table (Dong-Di Rienzo 2020 IEEE Access Table II):
# | omega | mu_r | Order needed |
# | 1 kHz | 1    | Leontovich   |
# | 50 Hz | 1    | Mitzner      |
# | 1 kHz | 100  | Mitzner      |
# | 1 kHz | 1000 | Mitzner not enough; order 3 |

# Golden criteria for the new `--ibc-order 2` panel mode:
#   order=1 vs order=2 difference: > 5% at any (f, mur) combo
#                                  in the Mitzner-needed band
#   order=2 vs full-volume FE reference: < 2% error
```

## Recipe 6: Per-panel curvature for HOIBC

`radia.panels.calc_inductance` provides no per-panel curvature extractor.
Use the geometric
`specialcf.Weingarten(3)` of BVP_2 on a curved (`mesh.Curve(order)`)
surface mesh; if a per-element mean curvature is needed, project
`Trace(specialcf.Weingarten(3))/2` = (kappa_1 + kappa_2)/2 onto
`SurfaceL2(mesh, order=0, definedon=mesh.Boundaries("conductor_bnd"))`.
Check it against the analytic mean curvature before using it: 1/R on a
sphere, 1/(2R) on a cylinder side (principal curvatures 1/R and 0).  The
sign follows the boundary normal: with the conductor meshed as a hole the
normal points into the conductor and the value changes sign (a spherical
hole gives -1/R).

## Cross-MCP recipe references

- `radia_mcp.peec.hoibc('dong_di_rienzo_2020')` -- formula reference
- `radia_mcp.peec.hoibc('mathematica_derivation')` -- symbolic check
- `radia_mcp.peec.hoibc('radia_application')` -- TODO checklist for
  full Radia HOIBC integration
- `radia_mcp.fem` -- FEM gauging + Kelvin transform for outer boundary
- `radia_mcp.matrix_solvers` -- Krylov solver and preconditioner selection
