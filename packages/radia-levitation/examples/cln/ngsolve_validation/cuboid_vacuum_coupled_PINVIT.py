"""5x2x1 mm Cu cuboid IN VACUUM - full FEM eigenvalue (CASE_VAC).

Independent verification of the "isolated cuboid in vacuum" eddy current
leading τ. Distinct from cuboid_521_PINVIT_caseA.py which uses
Dirichlet on conductor surface (= interior-only / PEC walls problem).

Setup:
- Mesh: conductor (5x2x1 mm) + surrounding vacuum box (e.g., 50x20x10 mm = 10x)
- HCurl space across BOTH regions (continuous A field)
- sigma field: sigma_Cu in conductor, 0 in vacuum
- Stiffness: (1/mu0) curl-curl over entire domain
- Mass: sigma × inner product (zero contribution in vacuum)
- BC: n × A = 0 on outer vacuum box boundary (truncates infinity)
- Eigenvalue: S A = λ M A, τ = 1/λ

Compare leading τ to:
  - 25.5 μs (interior Dirichlet analytical)
  - 19.77 μs (Phase 7-11 BEM Newton kernel)
  - To answer: which is the true vacuum-coupled τ?

Date: 2026-05-03
"""
from netgen.occ import Box, Pnt, OCCGeometry, Glue
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction, Preconditioner,
    CoefficientFunction, curl, dx, x, y, z, BitArray, BND,
    Integrate, TaskManager, ngsglobals,
)
from ngsolve.solvers import PINVIT
from math import pi
import time
import numpy as np

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7

# Conductor dimensions
ax_cond, ay_cond, az_cond = 5e-3, 2e-3, 1e-3

# Vacuum box dimensions (10x conductor)
VAC_FACTOR = 5  # vacuum box side is VAC_FACTOR times the largest conductor dim
ax_vac = VAC_FACTOR * ax_cond
ay_vac = VAC_FACTOR * ax_cond  # use largest dim for symmetry
az_vac = VAC_FACTOR * ax_cond

H_COND = 0.4e-3   # mesh size in conductor
H_VAC = 5.0e-3    # mesh size in vacuum (coarser)
ORDER = 2
N_EIGEN = 10
PINVIT_MAXIT = 600


def build_geometry():
    # Conductor centered at origin
    conductor = Box(
        Pnt(-ax_cond/2, -ay_cond/2, -az_cond/2),
        Pnt(+ax_cond/2, +ay_cond/2, +az_cond/2)
    )
    conductor.mat("conductor").bc("interface")
    conductor.maxh = H_COND

    # Vacuum box (centered at origin)
    vac_outer = Box(
        Pnt(-ax_vac/2, -ay_vac/2, -az_vac/2),
        Pnt(+ax_vac/2, +ay_vac/2, +az_vac/2)
    )
    vac_outer.mat("vacuum").bc("outer")
    vac_outer.maxh = H_VAC

    # Glue conductor and vacuum (subtract conductor from vacuum first)
    vacuum = vac_outer - conductor
    vacuum.mat("vacuum").bc("outer")
    vacuum.maxh = H_VAC

    geo = Glue([conductor, vacuum])
    return OCCGeometry(geo)


def main():
    ngsglobals.msg_level = 0
    print("=" * 64)
    print("Cuboid 5x2x1 mm Cu IN VACUUM - full FEM eigenvalue")
    print("=" * 64)
    print(f"  Conductor: {ax_cond*1000} x {ay_cond*1000} x {az_cond*1000} mm")
    print(f"  Vacuum box: {ax_vac*1000} x {ay_vac*1000} x {az_vac*1000} mm "
          f"({VAC_FACTOR}x conductor max dim)")
    print(f"  H_cond={H_COND*1000} mm, H_vac={H_VAC*1000} mm, order={ORDER}")

    print(f"\nGenerating mesh...")
    t0 = time.time()
    geo = build_geometry()
    mesh = Mesh(geo.GenerateMesh(maxh=H_VAC))
    print(f"  ne={mesh.ne}, nv={mesh.nv}, time={time.time()-t0:.1f}s")
    print(f"  materials: {mesh.GetMaterials()}")
    print(f"  boundaries: {mesh.GetBoundaries()}")

    # HCurl space, Dirichlet on outer boundary
    # nograds=True removes the gradient (curl=0) null space, improves conditioning
    fes = HCurl(mesh, order=ORDER, dirichlet="outer", nograds=True)
    print(f"\nHCurl ndof = {fes.ndof}")

    u, v = fes.TnT()

    # Stiffness: (1/mu0) curl-curl over entire mesh
    a = BilinearForm(fes)
    a += (1.0/mu0) * curl(u) * curl(v) * dx

    # Mass: sigma * inner product
    # sigma_vac is set small but nonzero to regularize M (avoid singular matrix)
    # For physical consistency, sigma_vac << sigma_Cu (e.g., 10^-6 ratio)
    sigma_vac = 1.0   # regularization, 10^7x smaller than sigma_Cu
    sigma_cf = CoefficientFunction(
        [sigma_Cu if mat == "conductor" else sigma_vac for mat in mesh.GetMaterials()]
    )
    m = BilinearForm(fes)
    m += sigma_cf * u * v * dx

    print("Assembling...")
    pre = Preconditioner(a, type="local")
    with TaskManager():
        a.Assemble()
        m.Assemble()
    print(f"  done")

    # PINVIT for smallest eigenvalues (= longest τ)
    print(f"\nRunning PINVIT for {N_EIGEN} smallest eigenpairs...")
    t0 = time.time()
    gf_modes = GridFunction(fes, multidim=N_EIGEN, name="modes")

    with TaskManager():
        evals = PINVIT(a.mat, m.mat, pre=pre.mat, num=N_EIGEN,
                       maxit=PINVIT_MAXIT, printrates=False)
    print(f"  PINVIT done in {time.time()-t0:.1f}s")

    # Parse eigenvalues
    if isinstance(evals, tuple):
        evals_vec = evals[0]
        evec_multi = evals[1]
    else:
        evals_vec = evals
        evec_multi = gf_modes.vecs

    if hasattr(evals_vec, "NumPy"):
        evals_arr = np.array(evals_vec.NumPy())
    elif hasattr(evals_vec, "FV"):
        evals_arr = np.array(list(evals_vec.FV()))
    else:
        s = str(evals_vec).strip()
        evals_arr = np.array([float(xv.strip()) for xv in s.split("\n") if xv.strip()])

    print(f"\nGeneralized eigenvalues (S v = λ M v): {evals_arr}")
    print(f"\nDiagnostic: τ_n raw, conductor energy fraction, β_n (H_z coupling):")
    A_ext = CoefficientFunction((-y/2, x/2, 0))
    sigma_cu_only = CoefficientFunction(
        [sigma_Cu if mat == "conductor" else 0.0 for mat in mesh.GetMaterials()]
    )
    sigma_vac_only = CoefficientFunction(
        [0.0 if mat == "conductor" else sigma_vac for mat in mesh.GetMaterials()]
    )
    for k in range(len(evals_arr)):
        ev = float(evals_arr[k])
        gf_k = GridFunction(fes)
        gf_k.vec.data = evec_multi[k]
        m_cu = float(Integrate(sigma_cu_only * gf_k * gf_k * dx, mesh))
        m_vac = float(Integrate(sigma_vac_only * gf_k * gf_k * dx, mesh))
        m_tot = m_cu + m_vac
        cu_frac = m_cu / m_tot if m_tot > 0 else 0.0
        # normalize
        if m_tot > 1e-30:
            gf_k.vec.data = (1.0/m_tot**0.5) * gf_k.vec
            beta = float(Integrate(sigma_cf * gf_k * A_ext * dx, mesh))
            beta_sq = beta**2
        else:
            beta_sq = 0.0
        tau_us = 1.0 / ev * 1e6 if abs(ev) > 1e-30 else float('inf')
        print(f"  mode {k:2d}: lam={ev:+.4e}, tau={tau_us:+.3e} us, "
              f"cu_frac={cu_frac:.3f}, beta2={beta_sq:.3e}")

    print()
    print("=" * 64)
    print("References:")
    print("  Interior Dirichlet TE_z(1,1,0) analytical: 25.5 μs")
    print("  Phase 7-11 BEM (parity-restricted):        19.77 μs")
    print("  HDiv div-free (our buggy):                 43.84 μs")
    print("=" * 64)


if __name__ == "__main__":
    main()
