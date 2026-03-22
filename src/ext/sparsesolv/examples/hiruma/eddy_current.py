"""
Eddy current analysis on Hiruma mesh.

Frequency-domain electromagnetic analysis of a copper coil around iron core.

Approach 1 (A-Phi + BDDC):
    Product space HCurl x H1(cond), exact coupling, NGSolve BDDC.
    Converges in ~3 CG iterations.

Approach 2 (Two-step + ICCG):
    Step a: Solve DC current distribution (Laplace for phi in conductor)
    Step b: HCurl A-only with pre-computed source, SparseSolv ICCG.
    Demonstrates that IC is ineffective for HCurl curl-curl with
    high-contrast materials (mur=1000 core).

Regions:
    core:  mur=1000, sigma=0      (ferromagnetic)
    cond:  mur=1,    sigma=5.96e7 (copper coil)
    air:   mur=1,    sigma=0
"""

import sys
import time
import numpy as np
from netgen.read_gmsh import ReadGmsh
from ngsolve import *
from ngsolve.krylovspace import CGSolver

# ---------- mesh ----------
mesh_name = "mesh1_2.5T"
m = ReadGmsh(mesh_name)
mesh = Mesh(m)
print(f"Mesh: {mesh_name}  nv={mesh.nv}  ne={mesh.ne}")
print(f"Materials: {mesh.GetMaterials()}")
print(f"Boundaries: {mesh.GetBoundaries()}")
sys.stdout.flush()

# ---------- physical parameters ----------
mu0 = 4e-7 * np.pi
freq = 50.0  # Hz
s = 1j * 2 * np.pi * freq

nu_cf = mesh.MaterialCF({
    "core": 1 / (1000 * mu0),
    "cond": 1 / mu0,
    "air":  1 / mu0,
})
sigma_cf = mesh.MaterialCF({"cond": 5.96e7}, default=0)

order = 1

# ============================================================
# Approach 1: A-Phi product space + NGSolve BDDC
# ============================================================
print("\n=== A-Phi + NGSolve BDDC ===")
sys.stdout.flush()

fesA = HCurl(mesh, order=order, nograds=True,
             dirichlet="dirichlet", complex=True)
fesPhi = H1(mesh, order=order, definedon="cond",
            dirichlet="gamma_in|gamma_out", complex=True)
fes = fesA * fesPhi
(A, phi), (N, psi) = fes.TnT()
print(f"DOFs: {fes.ndof:,}  (HCurl {fesA.ndof:,} + H1 {fesPhi.ndof:,})")

# Dirichlet: phi=1 on gamma_in, phi=0 on gamma_out
gf = GridFunction(fes)
gfA, gfPhi = gf.components
gfPhi.Set(1, definedon=mesh.Boundaries("gamma_in"))

a = BilinearForm(fes)
a += nu_cf * curl(A) * curl(N) * dx
a += 1e-6 * nu_cf * A * N * dx                                       # stabilization
a += s * sigma_cf * (A + grad(phi)) * (N + grad(psi)) * dx("cond")

c = Preconditioner(a, "bddc")
t0 = time.perf_counter()
a.Assemble()
c.Update()
t_asm = time.perf_counter() - t0

rhs = gf.vec.CreateVector()
rhs.data = -a.mat * gf.vec

t0 = time.perf_counter()
inv = CGSolver(a.mat, c.mat, maxiter=2000, tol=1e-8,
               printrates=False, conjugate=False)
gf.vec.data += inv * rhs
t_solve = time.perf_counter() - t0

gfA_sol, gfPhi_sol = gf.components
B1 = curl(gfA_sol)
B2_bddc = Integrate(B1 * Conj(B1), mesh).real
print(f"  assembly: {t_asm*1000:.0f} ms")
print(f"  solve:    {t_solve*1000:.0f} ms  ({inv.iterations} iters)")
print(f"  ||B||^2 = {B2_bddc:.6e}")
sys.stdout.flush()

# ============================================================
# Approach 2: Two-step A-only + SparseSolv ICCG
# ============================================================
print("\n=== Two-step A-only + SparseSolv ICCG ===")
sys.stdout.flush()

from sparsesolv_ngsolve import SparseSolvSolver

# Step a: DC current distribution (Laplace for phi in conductor)
print("  Step a: DC current distribution...")
fes_phi = H1(mesh, order=order, definedon="cond",
             dirichlet="gamma_in|gamma_out")
u_phi, v_phi = fes_phi.TnT()

a_phi = BilinearForm(fes_phi)
a_phi += sigma_cf * grad(u_phi) * grad(v_phi) * dx("cond")
a_phi.Assemble()

gf_phi = GridFunction(fes_phi)
gf_phi.Set(1, definedon=mesh.Boundaries("gamma_in"))

r_phi = gf_phi.vec.CreateVector()
r_phi.data = -a_phi.mat * gf_phi.vec
gf_phi.vec.data += a_phi.mat.Inverse(fes_phi.FreeDofs(),
                                      inverse="sparsecholesky") * r_phi

J_dc = sigma_cf * grad(gf_phi)
print(f"  ||J_dc||^2 = {Integrate(J_dc * J_dc, mesh):.4e}")

# Step b: HCurl A-only with pre-computed source
# curl(nu*curl A) + s*sigma*A = -s*sigma*grad(phi_dc)  in cond
# curl(nu*curl A) = 0                                   elsewhere
# No stabilization needed: IC auto_shift handles semi-definiteness
print("  Step b: HCurl solve with ICCG...")
fesA2 = HCurl(mesh, order=order, nograds=True,
              dirichlet="dirichlet", complex=True)
u2, v2 = fesA2.TnT()

a2 = BilinearForm(fesA2)
a2 += nu_cf * curl(u2) * curl(v2) * dx
a2 += s * sigma_cf * u2 * v2 * dx("cond")
a2.Assemble()

f2 = LinearForm(fesA2)
f2 += (-s * sigma_cf) * grad(gf_phi) * v2 * dx("cond")
f2.Assemble()

solver = SparseSolvSolver(a2.mat, method="ICCG",
                          freedofs=fesA2.FreeDofs(),
                          tol=1e-6, maxiter=5000,
                          save_best_result=True,
                          use_abmc=True)
solver.auto_shift = True

t0 = time.perf_counter()
gfu2 = GridFunction(fesA2)
gfu2.vec.data = solver * f2.vec
t_ic = time.perf_counter() - t0

result = solver.last_result
B2 = curl(gfu2)
B2_ic = Integrate(B2 * Conj(B2), mesh).real
print(f"  ICCG:   {t_ic*1000:.0f} ms  ({result.iterations} iters, "
      f"converged={result.converged})")
print(f"  ||B||^2 = {B2_ic:.6e}")
sys.stdout.flush()

# ============================================================
# Summary
# ============================================================
rel_diff = abs(B2_bddc - B2_ic) / max(abs(B2_bddc), 1e-30)
print(f"\n--- Summary ---")
print(f"BDDC:   {inv.iterations:5d} iters, ~{(t_asm + t_solve)*1000:.0f} ms total")
print(f"ICCG:   {result.iterations:5d} iters, ~{t_ic*1000:.0f} ms total")
print(f"||B||^2 relative diff: {rel_diff:.2e}")
print(f"\nBDDC is essential for HCurl problems with high-contrast materials.")

# ============================================================
# VTK Output (BDDC solution)
# ============================================================
print("\n=== VTK Output ===")

# Physical fields from BDDC A-Phi solution
B_re = curl(gfA_sol).real
B_im = curl(gfA_sol).imag
B_abs = sqrt(curl(gfA_sol).real * curl(gfA_sol).real
             + curl(gfA_sol).imag * curl(gfA_sol).imag)

# Current density in conductor: J = s*sigma*(A + grad(phi))
J_cond = s * sigma_cf * (gfA_sol + grad(gfPhi_sol))
J_re = J_cond.real
J_im = J_cond.imag

# Joule heat: Q = 0.5 * sigma * |E|^2 = 0.5 * Re(J . conj(E))
# E = -(s*A + grad(phi)*s) ... simplified: Q = 0.5/sigma * |J|^2
# where |J|^2 = J_re.J_re + J_im.J_im
J2 = J_re * J_re + J_im * J_im
Q_joule = 0.5 * IfPos(sigma_cf - 1, J2 / sigma_cf, CF(0))

# Material ID for visualization
mat_names = mesh.GetMaterials()
mat_id = mesh.MaterialCF({mat: i + 1 for i, mat in enumerate(mat_names)})

vtk = VTKOutput(mesh,
                coefs=[B_re, B_im, B_abs, J_re, J_im, Q_joule, mat_id],
                names=["B_re", "B_im", "B_abs", "J_re", "J_im", "Q_joule", "MaterialID"],
                filename=f"eddy_current_{mesh_name}",
                subdivision=0, legacy=False)
vtk.Do()
print(f"  Written: eddy_current_{mesh_name}.vtu")

# Total Joule power in conductor [W]
P_cond = Integrate(Q_joule, mesh, definedon=mesh.Materials("cond"))
print(f"  Joule power in conductor: {P_cond:.4e} W")

# Magnetic energy: W = 0.5 * Re(int nu * B . conj(B) dV)
W_mag = 0.5 * Integrate(nu_cf * B1 * Conj(B1), mesh).real
print(f"  Magnetic energy: {W_mag:.4e} J")
print(f"  (values are unnormalized; scale by (I_target/I_computed)^2)")
