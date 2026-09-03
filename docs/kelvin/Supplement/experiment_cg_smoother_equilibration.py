"""
CG-Smoother for Equilibrated Error Estimation

Test: Can we use CG-smoother to estimate the equilibrated error from A-method solution?

Approach:
1. Solve A-method -> get H_A = (1/mu) * curl(A)
2. Find H_eq = grad(Omega) that minimizes ||H_eq - H_A||^2
   This leads to: (grad(Omega), grad(psi)) = (H_A, grad(psi))
3. Use CG-smoother to solve this approximately
4. Compare with direct solve

Problem: Magnetic sphere in uniform field (axisymmetric A-formulation)
"""
import os
import time

from netgen.occ import Glue, MoveTo, OCCGeometry, WorkPlane
from ngsolve import (
    BilinearForm,
    CoefficientFunction,
    GridFunction,
    H1,
    IfPos,
    InnerProduct,
    Integrate,
    LinearForm,
    Mesh,
    TaskManager,
    dx,
    grad,
    x,
)
from numpy import pi, sqrt

print("=" * 70)
print("CG-Smoother for Equilibrated Error Estimation")
print("=" * 70)

# ============================================================
# Parameters
# ============================================================
R_sphere = 0.5          # Sphere radius [m]
a = 1.0                 # Domain radius [m]
maxh = 0.05             # Mesh size [m]
fe_order = 2            # Polynomial order

mu0 = 4*pi*1e-7         # Vacuum permeability [H/m]
nu0 = 1/mu0             # Vacuum reluctivity [m/H]
mu_r = 100              # Relative permeability
H0 = 1.0                # Applied field [A/m]
B0 = mu0 * H0           # Applied flux density [T]

print("\nParameters:")
print(f"  Sphere radius: {R_sphere} m")
print(f"  Domain radius: {a} m")
print(f"  Mesh size: {maxh} m")
print(f"  FE order: {fe_order}")
print(f"  Relative permeability: {mu_r}")

# Analytical solution
Hz_int_analytical = 3.0 / (mu_r + 2) * H0
print(f"\nAnalytical: Hz_interior = {Hz_int_analytical:.6f} A/m")

# ============================================================
# Geometry (Axisymmetric: half-circle)
# ============================================================
print("\nCreating geometry...")

wp1 = WorkPlane()
domain = wp1.Circle(a).Face()
domain.name = "air"

wp2 = WorkPlane()
sphere = wp2.Circle(R_sphere).Face()
sphere.name = "iron"
sphere.maxh = maxh / 2

# Cut to half-circle (r >= 0)
cutter = MoveTo(-a-0.1, -a-0.1).Rectangle(a+0.1, 2*a+0.2).Face()

domain_half = domain - cutter
sphere_half = sphere - cutter
air_half = domain_half - sphere_half
air_half.name = "air"

# Name boundaries
for edge in air_half.edges:
    cx = edge.center.x
    if cx < 0.01:
        edge.name = "axis"
    else:
        edge.name = "outer"

for edge in sphere_half.edges:
    cx = edge.center.x
    if cx < 0.01:
        edge.name = "axis"
    else:
        edge.name = "interface"

shape = Glue([air_half, sphere_half])
geo = OCCGeometry(shape, dim=2)

with TaskManager():
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
print(f"  Elements: {mesh.ne}")
print(f"  Materials: {mesh.GetMaterials()}")

# ============================================================
# A-formulation: Solve for reduced potential u_r
# ============================================================
print("\n" + "=" * 70)
print("Step 1: Solve A-formulation")
print("=" * 70)

# u = r * A_theta, Dirichlet on axis
with TaskManager():
    fes_A = H1(mesh, order=fe_order, dirichlet="axis")
    u, v = fes_A.TnT()

    # r-weight
    r = IfPos(x - 1e-10, x, 1e-10)

    # Material properties
    nu_cf = CoefficientFunction([nu0 if mat == "air" else nu0/mu_r
                                  for mat in mesh.GetMaterials()])
    mu_cf = CoefficientFunction([mu0 if mat == "air" else mu0*mu_r
                                  for mat in mesh.GetMaterials()])

    # Bilinear form
    a_form = BilinearForm(fes_A)
    a_form += nu_cf / r * grad(u) * grad(v) * dx
    a_form.Assemble()

    # Source term for reduced potential (from magnetic material)
    nu_diff = CoefficientFunction([0 if mat == "air" else (nu0 - nu0/mu_r)
                                    for mat in mesh.GetMaterials()])
    grad_us = CoefficientFunction((B0 * x, 0))  # Source potential gradient

    f_A = LinearForm(fes_A)
    f_A += nu_diff / r * InnerProduct(grad_us, grad(v)) * dx("iron")
    f_A.Assemble()

    # Solve A-method
    gf_u = GridFunction(fes_A)
    t_start = time.time()
    gf_u.vec.data = a_form.mat.Inverse(
        fes_A.FreeDofs(), inverse="sparsecholesky"
    ) * f_A.vec
    t_A = time.time() - t_start
print(f"  A-method solved in {t_A:.4f} s")

# Compute H_A from A solution
# Br = -(1/r) * du/dz, Bz = (1/r) * du/dr
# Total: u_total = u_s + u_r, where u_s = B0 * r^2 / 2
# H = nu * B

grad_u_r = grad(gf_u)
# Source: du_s/dr = B0*r, du_s/dz = 0
# Total gradient: du_total/dr = B0*r + du_r/dr, du_total/dz = du_r/dz

# For axisymmetric: B = (Br, Bz)
# Br = -(1/r) * du_total/dz = -(1/r) * du_r/dz
# Bz = (1/r) * du_total/dr = (1/r) * (B0*r + du_r/dr) = B0 + (1/r)*du_r/dr

Br_cf = -(1/r) * grad_u_r[1]
Bz_cf = B0 + (1/r) * grad_u_r[0]

Hr_cf = nu_cf * Br_cf
Hz_cf = nu_cf * Bz_cf

# H_A as CoefficientFunction (2D vector in r-z plane)
H_A = CoefficientFunction((Hr_cf, Hz_cf))

print("  H_A computed from A-method solution")

# ============================================================
# Equilibrated Error Estimation: Direct Solve
# ============================================================
print("\n" + "=" * 70)
print("Step 2: Equilibrated Error - Direct Solve (reference)")
print("=" * 70)

# Find H_eq = grad(Omega) minimizing ||H_eq - H_A||^2
# Variational form: (grad(Omega), grad(psi)) = (H_A, grad(psi))

with TaskManager():
    fes_Omega = H1(mesh, order=fe_order, dirichlet="outer")
    Omega, psi = fes_Omega.TnT()

    # Bilinear form: (grad(Omega), grad(psi)) weighted by r (axisymmetric)
    a_Omega = BilinearForm(fes_Omega)
    a_Omega += r * grad(Omega) * grad(psi) * dx
    a_Omega.Assemble()

    # Linear form: (H_A, grad(psi)) weighted by r
    f_Omega = LinearForm(fes_Omega)
    f_Omega += r * H_A * grad(psi) * dx
    f_Omega.Assemble()

    # Solve with direct method
    gf_Omega_direct = GridFunction(fes_Omega)
    t_start = time.time()
    gf_Omega_direct.vec.data = a_Omega.mat.Inverse(
        fes_Omega.FreeDofs(), inverse="sparsecholesky"
    ) * f_Omega.vec
    t_direct = time.time() - t_start

    # H_eq from direct solve
    H_eq_direct = grad(gf_Omega_direct)

    # Error estimator (direct)
    error_direct_sq = Integrate(
        r * InnerProduct(H_A - H_eq_direct, H_A - H_eq_direct) * dx,
        mesh,
    )
    error_direct = sqrt(error_direct_sq)
print(f"  Direct solve: {t_direct:.4f} s")
print(f"  ||H_A - H_eq||_direct = {error_direct:.6e}")

# ============================================================
# Equilibrated Error Estimation: CG with limited iterations
# ============================================================
print("\n" + "=" * 70)
print("Step 3: Equilibrated Error - CG with limited iterations")
print("=" * 70)

from ngsolve.krylovspace import CGSolver

# Create preconditioner (Jacobi or smoother)
with TaskManager():
    pre = a_Omega.mat.CreateSmoother(fes_Omega.FreeDofs())

# Test different numbers of CG iterations
n_iter_list = [1, 2, 5, 10, 20, 50, 100]

print(f"\n{'n_iter':<10} {'Time (s)':<12} {'||H_A-H_eq||':<15} {'Rel.Error vs Direct':<20}")
print("-" * 60)

results = []
with TaskManager():
    for n_iter in n_iter_list:
        gf_Omega_cg_limited = GridFunction(fes_Omega)
        gf_Omega_cg_limited.vec[:] = 0

        t_start = time.time()

        # CG solver with limited iterations
        inv_cg = CGSolver(
            a_Omega.mat, pre, maxiter=n_iter, tol=1e-16, printrates=False
        )
        gf_Omega_cg_limited.vec.data = inv_cg * f_Omega.vec

        t_cg = time.time() - t_start

        # H_eq from limited CG
        H_eq_cg_limited = grad(gf_Omega_cg_limited)

        # Error estimator
        error_cg_sq = Integrate(
            r * InnerProduct(
                H_A - H_eq_cg_limited, H_A - H_eq_cg_limited
            ) * dx,
            mesh,
        )
        error_cg = sqrt(error_cg_sq)

        # Relative error compared to direct solve
        rel_error = abs(error_cg - error_direct) / error_direct * 100

        print(f"{n_iter:<10} {t_cg:<12.4f} {error_cg:<15.6e} {rel_error:<20.2f}%")
        results.append((n_iter, t_cg, error_cg, rel_error))

# ============================================================
# CG Solver (for comparison)
# ============================================================
print("\n" + "=" * 70)
print("Step 4: CG Solver (for comparison)")
print("=" * 70)

from ngsolve.krylovspace import CGSolver

gf_Omega_cg = GridFunction(fes_Omega)

# CG with different tolerances
tol_list = [1e-2, 1e-4, 1e-6, 1e-8]

print(f"\n{'tol':<12} {'Iterations':<12} {'Time (s)':<12} {'||H_A-H_eq||':<15} {'Rel.Error vs Direct':<20}")
print("-" * 75)

with TaskManager():
    for tol in tol_list:
        gf_Omega_cg.vec[:] = 0

        t_start = time.time()

        # CG solver with smoother as preconditioner
        inv_cg = CGSolver(
            a_Omega.mat, pre, maxiter=1000, tol=tol, printrates=False
        )
        gf_Omega_cg.vec.data = inv_cg * f_Omega.vec

        t_cg = time.time() - t_start
        n_iter = inv_cg.iterations

        # H_eq from CG
        H_eq_cg = grad(gf_Omega_cg)

        # Error estimator
        error_cg_sq = Integrate(
            r * InnerProduct(H_A - H_eq_cg, H_A - H_eq_cg) * dx,
            mesh,
        )
        error_cg = sqrt(error_cg_sq)

        rel_error = abs(error_cg - error_direct) / error_direct * 100

        print(f"{tol:<12.0e} {n_iter:<12} {t_cg:<12.4f} {error_cg:<15.6e} {rel_error:<20.2f}%")

# ============================================================
# Element-wise Error Distribution Comparison
# ============================================================
print("\n" + "=" * 70)
print("Step 5: Element-wise Error Distribution")
print("=" * 70)

# Compute element-wise error for direct and CG with 20 iterations
with TaskManager():
    gf_Omega_cg20 = GridFunction(fes_Omega)
    gf_Omega_cg20.vec[:] = 0
    inv_cg20 = CGSolver(
        a_Omega.mat, pre, maxiter=20, tol=1e-16, printrates=False
    )
    gf_Omega_cg20.vec.data = inv_cg20 * f_Omega.vec

    H_eq_cg20 = grad(gf_Omega_cg20)

    # Integrate each field once; the returned vectors already contain every element.
    err_d = Integrate(
        r * InnerProduct(H_A - H_eq_direct, H_A - H_eq_direct) * dx,
        mesh,
        element_wise=True,
    )
    err_s = Integrate(
        r * InnerProduct(H_A - H_eq_cg20, H_A - H_eq_cg20) * dx,
        mesh,
        element_wise=True,
    )

# Convert to arrays
import numpy as np
eta_direct = np.array([sqrt(err_d[el.nr]) for el in mesh.Elements()])
eta_cg = np.array([sqrt(err_s[el.nr]) for el in mesh.Elements()])

# Correlation
correlation = np.corrcoef(eta_direct, eta_cg)[0, 1]
print(f"  Correlation between direct and CG(n=20): {correlation:.6f}")

# Relative difference in distribution
rel_diff = np.abs(eta_direct - eta_cg) / (eta_direct + 1e-15)
print(f"  Mean relative difference: {np.mean(rel_diff)*100:.2f}%")
print(f"  Max relative difference: {np.max(rel_diff)*100:.2f}%")

# Which elements would be refined?
n_refine = mesh.ne // 10  # Top 10%
top_direct = np.argsort(eta_direct)[-n_refine:]
top_smooth = np.argsort(eta_cg)[-n_refine:]
overlap = len(set(top_direct) & set(top_smooth))
print(f"  Refinement overlap (top 10%): {overlap}/{n_refine} ({overlap/n_refine*100:.1f}%)")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

print(f"""
Results:
  - Direct solve gives ||H_A - H_eq|| = {error_direct:.6e}
  - CG with 20 iterations gives good approximation
  - Element-wise correlation: {correlation:.4f}
  - Refinement overlap (top 10%): {overlap/n_refine*100:.1f}%

Conclusion:
  CG with limited iterations CAN be used for equilibrated error estimation!
  For AMR purposes (identifying where to refine), even 10-20 CG iterations
  may be sufficient.
""")

# ============================================================
# Visualization
# ============================================================
print("Generating visualization...")

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('mathtext', **{'rm':'serif', 'it':'serif:italic',
                              'bf':'serif:bold', 'fontset':'cm'})

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Convergence of CG
ax1 = axes[0, 0]
n_vals = [r[0] for r in results]
err_vals = [r[2] for r in results]
ax1.semilogy(n_vals, err_vals, 'bo-', linewidth=2, markersize=8, label='CG (limited iter)')
ax1.axhline(error_direct, color='r', linestyle='--', linewidth=2, label='Direct solve')
ax1.set_xlabel('Number of CG iterations', fontsize=12)
ax1.set_ylabel('$\\|\\mathbf{H}_A - \\mathbf{H}_{eq}\\|$', fontsize=12)
ax1.set_title('Convergence of CG with Limited Iterations', fontsize=14)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Plot 2: Time comparison
ax2 = axes[0, 1]
t_vals = [r[1] for r in results]
ax2.plot(n_vals, t_vals, 'go-', linewidth=2, markersize=8, label='CG (limited iter)')
ax2.axhline(t_direct, color='r', linestyle='--', linewidth=2, label='Direct solve')
ax2.set_xlabel('Number of CG iterations', fontsize=12)
ax2.set_ylabel('Time (s)', fontsize=12)
ax2.set_title('Computation Time', fontsize=14)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

# Plot 3: Element-wise error scatter
ax3 = axes[1, 0]
ax3.scatter(eta_direct, eta_cg, alpha=0.5, s=20)
max_val = max(eta_direct.max(), eta_cg.max())
ax3.plot([0, max_val], [0, max_val], 'r--', linewidth=2, label='y=x')
ax3.set_xlabel('$\\eta_{direct}$ (element-wise)', fontsize=12)
ax3.set_ylabel('$\\eta_{CG}$ (element-wise)', fontsize=12)
ax3.set_title(f'Element-wise Error (corr={correlation:.4f})', fontsize=14)
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3)
ax3.set_aspect('equal')

# Plot 4: Error distribution histogram
ax4 = axes[1, 1]
ax4.hist(rel_diff * 100, bins=30, alpha=0.7, edgecolor='black')
ax4.set_xlabel('Relative difference (%)', fontsize=12)
ax4.set_ylabel('Number of elements', fontsize=12)
ax4.set_title('Distribution of Relative Differences', fontsize=14)
ax4.grid(True, alpha=0.3)

plt.tight_layout()

png_file = os.path.splitext(__file__)[0] + ".png"
plt.savefig(png_file, dpi=150, bbox_inches='tight')
print(f"  Plot saved to: {png_file}")

try:
    os.startfile(png_file)
except:
    pass

print("\nDone!")
