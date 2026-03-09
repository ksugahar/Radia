"""
Nodal Force (Eggshell) Method for 2D Electromagnetic Force Computation
======================================================================
Compares 3 force calculation methods on a wire + iron cylinder problem:
  1. Analytical (multipole series, exact for linear material)
  2. Maxwell Stress Tensor (MST, point-wise integration on circle)
  3. Nodal Force / Virtual Work (Eggshell + mesh.SetDeformation)

Benchmark: straight wire (I=1000A) at origin, iron cylinder (mu_r=1000,
R=0.05m) at distance d=0.3m. Force is attractive (negative Fx).

IMPORTANT: mesh.Curve(3) is required for accurate force computation
with curved geometries. Without curving, the polygon approximation
of circles loses ~2% area, causing ~9% force error.

References:
  - Coulomb (1983): Virtual work principle in FEM
  - Henrotte & Hameyer (2004): Maxwell stress vs virtual work
  - Henrotte, Deliege & Hameyer (2003): Eggshell approach

Updated: 2026-03-09
"""
import os
import sys
import numpy as np
from math import pi

from ngsolve import *
from netgen.occ import *

# ============================================================
# Parameters
# ============================================================
mu0 = 4 * pi * 1e-7
mu_r = 1000.0           # Iron relative permeability
I_wire = 1000.0         # Wire current [A]
R_wire = 0.01           # Wire radius [m]
R_cyl = 0.05            # Iron cylinder radius [m]
d = 0.3                 # Distance wire-to-iron center [m]
t_egg = 0.008           # Eggshell thickness [m]
R_outer = 3.0           # Outer boundary radius [m]
J_z = I_wire / (pi * R_wire**2)  # Current density [A/m^2]

print("=" * 60)
print("Nodal Force (Eggshell) Method - 2D Benchmark")
print("=" * 60)
print(f"  Wire current:    I = {I_wire} A")
print(f"  Iron cylinder:   R = {R_cyl} m, mu_r = {mu_r}")
print(f"  Distance:        d = {d} m")
print(f"  Current density: J = {J_z:.2e} A/m^2")


# ============================================================
# Method 1: Analytical (Multipole Series, Exact)
# ============================================================
# Line current near a permeable cylinder (2D), closed-form:
# F_x = -(mu_0*I^2*alpha*R^2) / (2*pi*d*(d^2-R^2))
# where alpha = (mu_r - 1) / (mu_r + 1)
alpha = (mu_r - 1) / (mu_r + 1)
F_analytical = -mu0 * I_wire**2 * alpha * R_cyl**2 / (2 * pi * d * (d**2 - R_cyl**2))
F_dipole = -mu0 * alpha * I_wire**2 * R_cyl**2 / (2 * pi * d**3)

print(f"\n--- Method 1: Analytical ---")
print(f"  F_dipole (n=1 only) = {F_dipole:.6e} N/m")
print(f"  F_exact  (closed)   = {F_analytical:.6e} N/m")


def compute_force_fem(maxh_iron, order=3, curve_order=3):
    """Compute force using MST and nodal force for a given mesh size.

    Parameters
    ----------
    maxh_iron : float
        Maximum element size in iron and eggshell regions.
    order : int
        Polynomial order for the H1 finite element space.
    curve_order : int
        Order of mesh curving. Use >= 2 for curved boundaries.
    """

    # ============================================================
    # Geometry (OCC 2D)
    # ============================================================
    wp = WorkPlane()

    # Wire at origin
    wire_shape = wp.MoveTo(0, 0).Circle(R_wire).Face()
    wire_shape.name = "wire"
    wire_shape.maxh = R_wire / 2

    # Iron cylinder at (d, 0)
    iron_shape = wp.MoveTo(d, 0).Circle(R_cyl).Face()
    iron_shape.name = "iron"
    iron_shape.maxh = maxh_iron

    # Eggshell: thin air ring around iron
    R_egg_outer = R_cyl + t_egg
    egg_outer_shape = wp.MoveTo(d, 0).Circle(R_egg_outer).Face()
    eggshell_shape = egg_outer_shape - iron_shape
    eggshell_shape.name = "eggshell"
    eggshell_shape.maxh = maxh_iron

    # Outer air domain
    outer_shape = wp.MoveTo(0, 0).Circle(R_outer).Face()
    for edge in outer_shape.edges:
        edge.name = "outer"

    # Air = outer - wire - eggshell outer circle
    air_shape = outer_shape - wire_shape - egg_outer_shape
    air_shape.name = "air"

    # Glue all together and generate mesh
    geo = Glue([wire_shape, iron_shape, eggshell_shape, air_shape])
    mesh = Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=R_outer / 5))

    # Curve mesh for accurate representation of circular boundaries
    if curve_order > 1:
        mesh.Curve(curve_order)

    # ============================================================
    # A-formulation Magnetostatic Solve
    # ============================================================
    fes = H1(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()

    mu_cf = mesh.MaterialCF({"iron": mu0 * mu_r}, default=mu0)
    nu_cf = 1.0 / mu_cf

    a = BilinearForm(nu_cf * grad(u) * grad(v) * dx)
    a.Assemble()

    f = LinearForm(J_z * v * dx("wire"))
    f.Assemble()

    gf_A = GridFunction(fes)
    gf_A.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    # B field from A: B = curl(A_z) = (-dA/dy, dA/dx) -- sign irrelevant for MST
    B = CF((-grad(gf_A)[1], grad(gf_A)[0]))

    # ============================================================
    # Method 2: Maxwell Stress Tensor (point-wise on circle)
    # ============================================================
    R_mst = R_cyl + t_egg / 2  # Integration circle in eggshell (air)
    n_quad = 360
    theta = np.linspace(0, 2 * pi, n_quad + 1)[:-1]
    dtheta = 2 * pi / n_quad

    Fx_MST = 0.0
    for k in range(n_quad):
        x_q = d + R_mst * np.cos(theta[k])
        y_q = R_mst * np.sin(theta[k])
        nx = np.cos(theta[k])
        ny = np.sin(theta[k])

        mip = mesh(x_q, y_q)
        Bx_val = B[0](mip)
        By_val = B[1](mip)
        Bn = Bx_val * nx + By_val * ny
        B2 = Bx_val**2 + By_val**2

        Fx_MST += (1.0 / mu0) * (Bx_val * Bn - 0.5 * B2 * nx) * R_mst * dtheta

    # ============================================================
    # Method 3: Nodal Force (Eggshell + Virtual Work)
    # ============================================================
    # Blending function: 1 inside iron, smooth transition to 0 in eggshell
    fes_blend = H1(mesh, order=1)
    gf_blend = GridFunction(fes_blend)
    gf_blend.Set(1.0, definedon=mesh.Materials("iron"))

    # Virtual displacement field: V = (blend * delta, 0)
    fes_def = VectorH1(mesh, order=1)
    gf_def = GridFunction(fes_def)
    delta = 1e-7

    # Magnetic energy density: W = 0.5 * nu * |grad(A)|^2
    energy_density = 0.5 * nu_cf * (grad(gf_A)[0]**2 + grad(gf_A)[1]**2)

    def compute_energy():
        return Integrate(energy_density, mesh, order=2 * order + 4)

    # Forward displacement (+delta in x)
    gf_def.components[0].vec.data = delta * gf_blend.vec
    gf_def.components[1].vec[:] = 0.0
    mesh.SetDeformation(gf_def)
    W_plus = compute_energy()
    mesh.UnsetDeformation()

    # Backward displacement (-delta in x)
    gf_def.components[0].vec.data = (-delta) * gf_blend.vec
    gf_def.components[1].vec[:] = 0.0
    mesh.SetDeformation(gf_def)
    W_minus = compute_energy()
    mesh.UnsetDeformation()

    # Central difference: F_x = -dW/dx
    Fx_nodal = -(W_plus - W_minus) / (2 * delta)

    return mesh.ne, Fx_MST, Fx_nodal


# ============================================================
# Convergence Study
# ============================================================
print("\n--- Convergence Study (mesh.Curve(3) for accurate geometry) ---")
print(f"{'maxh':>8s}  {'NE':>6s}  {'F_MST':>14s}  {'F_nodal':>14s}  "
      f"{'err_MST%':>9s}  {'err_nodal%':>11s}")
print("-" * 75)

maxh_list = [0.02, 0.012, 0.008, 0.005, 0.003]
results = []

for maxh in maxh_list:
    ne, Fx_MST, Fx_nodal = compute_force_fem(maxh, order=3, curve_order=3)
    err_MST = (Fx_MST - F_analytical) / abs(F_analytical) * 100
    err_nodal = (Fx_nodal - F_analytical) / abs(F_analytical) * 100
    results.append((maxh, ne, Fx_MST, Fx_nodal, err_MST, err_nodal))
    print(f"{maxh:8.4f}  {ne:6d}  {Fx_MST:14.6e}  {Fx_nodal:14.6e}  "
          f"{err_MST:+9.2f}  {err_nodal:+11.2f}")

print(f"\n  Analytical (exact): F_x = {F_analytical:.6e} N/m")

# MST vs nodal agreement at finest mesh
_, _, Fx_MST_fine, Fx_nodal_fine, _, _ = results[-1]
mst_nodal_diff = abs(Fx_MST_fine - Fx_nodal_fine) / abs(F_analytical) * 100
print(f"  MST vs Nodal difference (finest): {mst_nodal_diff:.4f}%")

# ============================================================
# Visualization
# ============================================================
print("\nGenerating convergence plot...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm': 'serif', 'it': 'serif:italic',
                             'bf': 'serif:bold', 'fontset': 'cm'})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=150)

ne_arr = [r[1] for r in results]
err_MST_arr = [abs(r[4]) for r in results]
err_nodal_arr = [abs(r[5]) for r in results]
F_MST_arr = [r[2] for r in results]
F_nodal_arr = [r[3] for r in results]

# Left: Force values vs number of elements
ax1.axhline(y=F_analytical, color='gray', linestyle='--', linewidth=1,
            label='Analytical (exact)')
ax1.plot(ne_arr, F_MST_arr, 'bo-', linewidth=1.5, markersize=6, label='MST')
ax1.plot(ne_arr, F_nodal_arr, 'rs-', linewidth=1.5, markersize=6, label='Nodal force')
ax1.set_xlabel('Number of elements', fontname='Times New Roman', fontsize=11)
ax1.set_ylabel('$F_x$ (N/m)', fontname='Times New Roman', fontsize=11)
ax1.set_title('Force convergence', fontname='Times New Roman', fontsize=12)
ax1.legend(fontsize=9, frameon=False)
ax1.minorticks_on()
ax1.tick_params(which='major', direction='in', top=True, right=True)
ax1.tick_params(which='minor', direction='in', top=True, right=True)
ax1.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)

# Right: Relative error vs number of elements (log scale)
ax2.loglog(ne_arr, err_MST_arr, 'bo-', linewidth=1.5, markersize=6, label='MST')
ax2.loglog(ne_arr, err_nodal_arr, 'rs-', linewidth=1.5, markersize=6, label='Nodal force')
ax2.set_xlabel('Number of elements', fontname='Times New Roman', fontsize=11)
ax2.set_ylabel('Relative error vs analytical (%)', fontname='Times New Roman', fontsize=11)
ax2.set_title('Error convergence', fontname='Times New Roman', fontsize=12)
ax2.legend(fontsize=9, frameon=False)
ax2.minorticks_on()
ax2.tick_params(which='major', direction='in', top=True, right=True)
ax2.tick_params(which='minor', direction='in', top=True, right=True)
ax2.grid(axis='both', which='major', c='gainsboro', linestyle=':', linewidth=0.3)

plt.tight_layout()

png_file = os.path.join(os.path.dirname(__file__), "demo_nodal_force_2d.png")
plt.savefig(png_file, dpi=150, bbox_inches='tight')
print(f"  Plot saved to: {png_file}")

print("\n" + "=" * 60)
print("Computation completed successfully")
print("=" * 60)
