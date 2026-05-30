"""
Magnetostatic Problem: Magnetic Sphere in Uniform External Field

A sphere of relative permeability mu_r placed in a uniform external magnetic field H0.

Analytical Solution:
  Inside sphere (r < R):
    H_in = 3 / (mu_r + 2) * H0  (uniform, parallel to H0)
    B_in = mu0 * mu_r * H_in

  Outside sphere (r > R):
    Perturbed field with dipole moment m = 4*pi*R^3 * (mu_r-1)/(mu_r+2) * H0

Reference: Jackson, Classical Electrodynamics, Section 5.11

This script compares:
1. A method (vector potential formulation): curl(1/mu * curl(A)) = 0
2. Omega method (scalar potential formulation): div(mu * grad(Omega)) = 0

Both methods should give identical results within discretization error.

Equilibrated Error Estimation (Prager-Synge):
  ||H - H_A|| <= ||H_A - H_eq||
  where H_eq satisfies curl(H_eq) = 0.

  For the A method, H_A = B/mu = curl(A)/mu does NOT satisfy curl(H_A) = 0 in general.
  The best equilibrated field H_eq is the Omega method solution: H_Omega = grad(Omega).

  Therefore:
  ||H_A - H_Omega||_mu gives the Prager-Synge error bound for the A method.
"""
import os
import sys

os.environ['NETGENDIR'] = r'S:\NGSolve\01_GitHub\install_ngsolve\bin'
os.environ['PATH'] = r'S:\NGSolve\01_GitHub\install_ngsolve\bin;' + os.environ.get('PATH', '')
sys.path.insert(0, r'S:\NGSolve\01_GitHub\install_ngsolve\Lib\site-packages')

from numpy import pi, sqrt
from ngsolve import *
from netgen.occ import *

print("=" * 70)
print("Magnetic Sphere in Uniform External Field")
print("=" * 70)

# ============================================================
# Parameters
# ============================================================
R_sphere = 0.25    # Sphere radius [m]
R_outer = 1.0      # Outer boundary radius [m]
mu_r = 100         # Relative permeability of sphere
mu0 = 4 * pi * 1e-7  # Vacuum permeability [H/m]
H0 = 1000.0        # External field strength [A/m]

# Analytical solution
H_in_exact = 3.0 / (mu_r + 2) * H0
B_in_exact = mu0 * mu_r * H_in_exact

print(f"\nParameters:")
print(f"  Sphere radius: {R_sphere} m")
print(f"  Outer radius: {R_outer} m")
print(f"  Relative permeability: {mu_r}")
print(f"  External field: H0 = {H0} A/m")

print(f"\nAnalytical Solution (inside sphere):")
print(f"  H_in = {H_in_exact:.4f} A/m")
print(f"  B_in = {B_in_exact*1e3:.4f} mT")

# ============================================================
# Geometry
# ============================================================
print("\nCreating geometry...")

iron = Sphere(Pnt(0, 0, 0), R_sphere).mat("iron")
outer = Sphere(Pnt(0, 0, 0), R_outer)
outer.faces.name = "outer"
air = (outer - iron).mat("air")

geo = Glue([iron, air])
mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=0.1))
# Important: Curve order must match FES polynomial order for good convergence
with TaskManager():
    mesh.Curve(2)
    print(f"  Mesh elements: {mesh.ne}")

    # Material properties
    mu_cf = mesh.MaterialCF({"iron": mu_r * mu0, "air": mu0}, default=mu0)
    inv_mu_cf = mesh.MaterialCF({"iron": 1/(mu_r * mu0), "air": 1/mu0}, default=1/mu0)

    # Test point inside sphere
    test_point = mesh(0.1, 0.1, 0.1)

    # ============================================================
    # Method 1: A-method (Vector Potential)
    # ============================================================
    print("\n" + "=" * 60)
    print("A Method (Vector Potential)")
    print("=" * 60)

    # External field B0 = mu0 * H0 in z-direction
    # Vector potential: A = (B0/2)(-y, x, 0)
    B0 = mu0 * H0
    Av = CF((B0/2 * (-y), B0/2 * x, 0))

    fes_A = HCurl(mesh, order=2, dirichlet="outer", nograds=True)
    A, v_A = fes_A.TnT()

    a_A = BilinearForm(fes_A)
    a_A += inv_mu_cf * InnerProduct(curl(A), curl(v_A)) * dx
    a_A.Assemble()

    gf_A = GridFunction(fes_A)
    gf_A.Set(Av, BND)

    f_A = LinearForm(fes_A)
    f_A.Assemble()

    res_A = gf_A.vec.CreateVector()
    res_A.data = f_A.vec - a_A.mat * gf_A.vec

    from ngsolve import CGSolver
    from ngsolve import TaskManager
    pre_A = a_A.mat.CreateSmoother(fes_A.FreeDofs())
    inv_A = CGSolver(a_A.mat, pre_A, maxsteps=500, printrates=False, precision=1e-10)
    update_A = gf_A.vec.CreateVector()
    update_A.data = inv_A * res_A
    gf_A.vec.data += update_A

    B_from_A = curl(gf_A)
    H_from_A = inv_mu_cf * B_from_A

    H_A_val = H_from_A(test_point)
    print(f"\nA method at (0.1, 0.1, 0.1):")
    print(f"  H = ({H_A_val[0]:.4f}, {H_A_val[1]:.4f}, {H_A_val[2]:.4f}) A/m")
    print(f"  H_z error: {(H_A_val[2] - H_in_exact)/H_in_exact*100:.2f}%")

    # ============================================================
    # Method 2: Omega-method (Scalar Potential)
    # ============================================================
    print("\n" + "=" * 60)
    print("Omega Method (Scalar Potential)")
    print("=" * 60)

    # H = -grad(Omega) where div(mu * grad(Omega)) = 0
    # BC: Omega = -H0 * z on outer (gives H -> H0 * z_hat at infinity)

    fes_O = H1(mesh, order=2, dirichlet="outer")
    Omega, psi = fes_O.TnT()

    a_O = BilinearForm(fes_O)
    a_O += mu_cf * grad(Omega) * grad(psi) * dx
    a_O.Assemble()

    gf_Omega = GridFunction(fes_O)
    gf_Omega.Set(H0 * z, BND)  # H = grad(Omega) -> H0*z at boundary

    f_O = LinearForm(fes_O)
    f_O.Assemble()

    res_O = gf_Omega.vec.CreateVector()
    res_O.data = f_O.vec - a_O.mat * gf_Omega.vec
    gf_Omega.vec.data += a_O.mat.Inverse(fes_O.FreeDofs()) * res_O

    H_from_Omega = grad(gf_Omega)

    H_O_val = H_from_Omega(test_point)
    print(f"\nOmega method at (0.1, 0.1, 0.1):")
    print(f"  H = ({H_O_val[0]:.4f}, {H_O_val[1]:.4f}, {H_O_val[2]:.4f}) A/m")
    print(f"  H_z error: {(H_O_val[2] - H_in_exact)/H_in_exact*100:.2f}%")

    # ============================================================
    # H1 Gradient Recovery from A Method
    # ============================================================
    print("\n" + "=" * 60)
    print("H1 Gradient Recovery from A Method")
    print("=" * 60)

    # Recover curl-free H field from B_from_A
    # Solve: <mu * grad(Omega_rec), grad(psi)> = <B, grad(psi)>

    fes_rec = H1(mesh, order=3, dirichlet="outer")
    Omega_rec, psi_rec = fes_rec.TnT()

    gf_Omega_rec = GridFunction(fes_rec)
    gf_Omega_rec.Set(H0 * z, BND)

    a_rec = BilinearForm(fes_rec)
    a_rec += mu_cf * grad(Omega_rec) * grad(psi_rec) * dx
    a_rec.Assemble()

    f_rec = LinearForm(fes_rec)
    f_rec += InnerProduct(B_from_A, grad(psi_rec)) * dx
    f_rec.Assemble()

    res_rec = gf_Omega_rec.vec.CreateVector()
    res_rec.data = f_rec.vec - a_rec.mat * gf_Omega_rec.vec
    gf_Omega_rec.vec.data += a_rec.mat.Inverse(fes_rec.FreeDofs()) * res_rec

    H_rec = grad(gf_Omega_rec)

    H_rec_val = H_rec(test_point)
    print(f"\nH1 recovery at (0.1, 0.1, 0.1):")
    print(f"  H = ({H_rec_val[0]:.4f}, {H_rec_val[1]:.4f}, {H_rec_val[2]:.4f}) A/m")
    print(f"  H_z error: {(H_rec_val[2] - H_in_exact)/H_in_exact*100:.2f}%")

    # ============================================================
    # Error Estimation (Prager-Synge)
    # ============================================================
    print("\n" + "=" * 70)
    print("PRAGER-SYNGE ERROR ESTIMATION")
    print("=" * 70)

    # Full Prager-Synge bound: ||H_A - H_Omega||_mu
    H_diff_AO = H_from_A - H_from_Omega
    eta_AO = sqrt(Integrate(mu_cf * InnerProduct(H_diff_AO, H_diff_AO) * dx, mesh))
    print(f"\nFull Prager-Synge bound:")
    print(f"  ||H_A - H_Omega||_mu = {eta_AO:.6e}")

    # Recovery-based estimator: ||H_A - H_rec||_mu
    H_diff_AR = H_from_A - H_rec
    eta_AR = sqrt(Integrate(mu_cf * InnerProduct(H_diff_AR, H_diff_AR) * dx, mesh))
    print(f"\nRecovery-based estimator (from A method only):")
    print(f"  ||H_A - H_rec||_mu = {eta_AR:.6e}")
    print(f"  Ratio (rec/full): {eta_AR/eta_AO:.4f}")

    # Difference between recovery and Omega
    H_diff_RO = H_rec - H_from_Omega
    eta_RO = sqrt(Integrate(mu_cf * InnerProduct(H_diff_RO, H_diff_RO) * dx, mesh))
    print(f"\nRecovery quality:")
    print(f"  ||H_rec - H_Omega||_mu = {eta_RO:.6e}")

    # ============================================================
    # Summary
    # ============================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"""
    Analytical: H_in = {H_in_exact:.4f} A/m

    Method                    H_z [A/m]    Error [%]
    ------------------------------------------------
    A method (B/mu)           {H_A_val[2]:.4f}       {(H_A_val[2]-H_in_exact)/H_in_exact*100:.2f}
    Omega method              {H_O_val[2]:.4f}       {(H_O_val[2]-H_in_exact)/H_in_exact*100:.2f}
    H1 recovery from A        {H_rec_val[2]:.4f}       {(H_rec_val[2]-H_in_exact)/H_in_exact*100:.2f}

    Prager-Synge Error Bounds (mu-weighted norm):
      Full bound ||H_A - H_Omega||_mu = {eta_AO:.6e}
      Recovery   ||H_A - H_rec||_mu   = {eta_AR:.6e} ({eta_AR/eta_AO*100:.1f}%)

    Key insight:
      The A method gives curl(H_A) != 0 in general.
      The Omega method gives curl(H_Omega) = 0 exactly.
      For Prager-Synge equilibrated estimation, we need curl(H_eq) = 0.

      Therefore:
      - If only A method is solved: use H1 recovery (ratio {eta_AR/eta_AO:.2f})
      - If both methods are solved: use ||H_A - H_Omega|| (full bound)
    """)

    print("\nDone.")
