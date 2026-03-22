"""
Shifted AMS Experiment: solve curl-curl with air WITHOUT regularization.

Idea: Add epsilon*mass ONLY to the preconditioner, not to the system matrix.
The Krylov solver sees the original (singular) system.
The preconditioner sees a shifted (non-singular) system.

Compare:
  1. Standard regularized: solve (A + eps*M) x = b, precond = BDDC(A + eps*M)
  2. Shifted precond:      solve A x = b,           precond = BDDC(A + eps*M)
  3. No shift (baseline):  solve A x = b,           precond = BDDC(A)
"""

import sys
import time
import numpy as np

from ngsolve import *
from ngsolve.krylovspace import CGSolver
from netgen.csg import CSGeometry, Pnt, Cylinder, OrthoBrick, Sphere

# ---------- geometry: conducting sphere in air ----------
geo = CSGeometry()
sphere = Sphere(Pnt(0, 0, 0), 0.3).bc("interface")
box = OrthoBrick(Pnt(-1, -1, -1), Pnt(1, 1, 1)).bc("outer")
geo.Add(sphere.mat("cond"))
geo.Add((box - sphere).mat("air"))

mesh = Mesh(geo.GenerateMesh(maxh=0.3))
print(f"Mesh: nv={mesh.nv} ne={mesh.ne}")
print(f"Materials: {mesh.GetMaterials()}")

# ---------- parameters ----------
mu0 = 4e-7 * np.pi
freq = 1e3
omega = 2 * np.pi * freq
sigma_cu = 5.8e7

nu_cf = 1.0 / mu0
sigma_cf = mesh.MaterialCF({"cond": sigma_cu}, default=0)

order = 1
B0 = 0.01
A0 = CoefficientFunction((0.5 * B0 * y, -0.5 * B0 * x, 0))

eps_values = [1e-4, 1e-6, 1e-8]

print(f"\n{'='*60}")
print(f"  Shifted AMS Experiment (freq={freq} Hz)")
print(f"{'='*60}")


def solve_with_precond(label, a_sys, a_pre, fes):
    """Solve a_sys * x = b using preconditioner from a_pre."""
    u, v = fes.TnT()

    # Preconditioner must be registered BEFORE Assemble
    c = Preconditioner(a_pre, "bddc")
    a_pre.Assemble()
    if a_sys is not a_pre:
        a_sys.Assemble()

    gf_bc = GridFunction(fes)
    gf_bc.Set(A0, definedon=mesh.Boundaries("outer"))

    rhs = gf_bc.vec.CreateVector()
    rhs.data = -a_sys.mat * gf_bc.vec

    t0 = time.perf_counter()
    inv = CGSolver(a_sys.mat, c.mat, maxiter=2000, tol=1e-10,
                   printrates=False, conjugate=False)
    gf_sol = GridFunction(fes)
    gf_sol.vec.data = gf_bc.vec
    gf_sol.vec.data += inv * rhs
    t_solve = time.perf_counter() - t0

    B_sol = curl(gf_sol)
    B2 = Integrate(B_sol * Conj(B_sol), mesh).real
    print(f"  {label}: {inv.iterations} iter, {t_solve*1000:.0f} ms, ||B||^2 = {B2:.6e}")
    return B2, inv.iterations


# ============================================================
# Method 1: Standard regularized (eps in BOTH system and precond)
# ============================================================
print("\n--- Method 1: Regularized (eps in system + precond) ---")
for eps in eps_values:
    fes = HCurl(mesh, order=order, dirichlet="outer", complex=True)
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu_cf * curl(u) * curl(v) * dx
    a += 1j * omega * sigma_cf * u * v * dx("cond")
    a += eps * nu_cf * u * v * dx
    try:
        solve_with_precond(f"eps={eps:.0e}", a, a, fes)
    except Exception as e:
        print(f"  eps={eps:.0e}: FAILED - {e}")

# ============================================================
# Method 2: Shifted precond (eps in precond ONLY)
# ============================================================
print("\n--- Method 2: Shifted precond (eps in precond only) ---")
for eps in eps_values:
    fes = HCurl(mesh, order=order, dirichlet="outer", complex=True)
    u, v = fes.TnT()

    # Original system (no regularization)
    a_orig = BilinearForm(fes)
    a_orig += nu_cf * curl(u) * curl(v) * dx
    a_orig += 1j * omega * sigma_cf * u * v * dx("cond")

    # Shifted system (for precond only)
    a_shift = BilinearForm(fes)
    a_shift += nu_cf * curl(u) * curl(v) * dx
    a_shift += 1j * omega * sigma_cf * u * v * dx("cond")
    a_shift += eps * nu_cf * u * v * dx

    try:
        solve_with_precond(f"eps={eps:.0e}", a_orig, a_shift, fes)
    except Exception as e:
        print(f"  eps={eps:.0e}: FAILED - {e}")

# ============================================================
# Method 3: No regularization, no shift (baseline)
# ============================================================
print("\n--- Method 3: No regularization (baseline) ---")
fes = HCurl(mesh, order=order, dirichlet="outer", complex=True)
u, v = fes.TnT()
a_none = BilinearForm(fes)
a_none += nu_cf * curl(u) * curl(v) * dx
a_none += 1j * omega * sigma_cf * u * v * dx("cond")
try:
    solve_with_precond("no shift", a_none, a_none, fes)
except Exception as e:
    print(f"  no shift: FAILED - {e}")

print(f"\n{'='*60}")
print("  If Method 2 matches Method 1 ||B||^2 -> Shifted AMS works")
print(f"{'='*60}")
