"""5x2x1 mm Cu cuboid IN VACUUM - scipy eigsh with shift-invert.

Robust solver for the rank-deficient M problem (sigma_vac->0).
Use shift-invert at target lambda = 4e4 (corresponding to tau ~ 25 us)
so we hone in on conductor modes, not vacuum garbage.

Date: 2026-05-03
"""
from netgen.occ import Box, Pnt, OCCGeometry, Glue
from ngsolve import (
    Mesh, HCurl, BilinearForm, GridFunction, CoefficientFunction,
    curl, dx, x, y, z, Integrate, TaskManager, ngsglobals,
)
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import numpy as np
from math import pi
import time

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7

ax_cond, ay_cond, az_cond = 5e-3, 2e-3, 1e-3
VAC_FACTOR = 5
ax_vac = VAC_FACTOR * ax_cond
ay_vac = VAC_FACTOR * ax_cond
az_vac = VAC_FACTOR * ax_cond

H_COND = 0.4e-3
H_VAC = 5.0e-3
ORDER = 2
SIGMA_VAC = 1.0    # vacuum modes get HUGE M (rel to S), so lambda_vac is small. We want SMALLEST=conductor.
N_TARGET = 12
SHIFT_TARGET = None  # use which='SA' (smallest algebraic) instead


def build_geometry():
    conductor = Box(
        Pnt(-ax_cond/2, -ay_cond/2, -az_cond/2),
        Pnt(+ax_cond/2, +ay_cond/2, +az_cond/2)
    )
    conductor.mat("conductor").bc("interface")
    conductor.maxh = H_COND

    vac_outer = Box(
        Pnt(-ax_vac/2, -ay_vac/2, -az_vac/2),
        Pnt(+ax_vac/2, +ay_vac/2, +az_vac/2)
    )
    vac_outer.mat("vacuum").bc("outer")
    vac_outer.maxh = H_VAC

    vacuum = vac_outer - conductor
    vacuum.mat("vacuum").bc("outer")
    vacuum.maxh = H_VAC

    return OCCGeometry(Glue([conductor, vacuum]))


def to_scipy_csr(mat):
    rows, cols, vals = mat.COO()
    return sp.csr_matrix((vals, (rows, cols)), shape=(mat.height, mat.width))


def main():
    ngsglobals.msg_level = 0
    print("=" * 64)
    if SHIFT_TARGET is None:
        print(f"Cuboid 5x2x1 Cu in vacuum - scipy eigsh which='SA' (smallest algebraic)")
    else:
        print(f"Cuboid 5x2x1 Cu in vacuum - scipy eigsh shift-invert at lambda={SHIFT_TARGET:.1e}")
    print(f"sigma_vac={SIGMA_VAC}, sigma_Cu={sigma_Cu}, ratio={sigma_Cu/SIGMA_VAC:.1e}")
    print("=" * 64)

    print("Generating mesh...")
    t0 = time.time()
    geo = build_geometry()
    mesh = Mesh(geo.GenerateMesh(maxh=H_VAC))
    print(f"  ne={mesh.ne}, nv={mesh.nv}, time={time.time()-t0:.1f}s")
    print(f"  materials: {mesh.GetMaterials()}")

    fes = HCurl(mesh, order=ORDER, dirichlet="outer", nograds=True)
    print(f"  HCurl ndof = {fes.ndof}")

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += (1.0/mu0) * curl(u) * curl(v) * dx

    sigma_cf = CoefficientFunction(
        [sigma_Cu if mat == "conductor" else SIGMA_VAC for mat in mesh.GetMaterials()]
    )
    m = BilinearForm(fes)
    m += sigma_cf * u * v * dx

    print("Assembling...")
    with TaskManager():
        a.Assemble()
        m.Assemble()

    print("Converting to scipy CSR...")
    Ssp = to_scipy_csr(a.mat)
    Msp = to_scipy_csr(m.mat)
    print(f"  S nnz={Ssp.nnz}, M nnz={Msp.nnz}, shape={Ssp.shape}")

    # Symmetrize for safety
    Ssp = (Ssp + Ssp.T) / 2.0
    Msp = (Msp + Msp.T) / 2.0

    print(f"\nRunning scipy eigsh which='SA' (smallest algebraic)...")
    t0 = time.time()
    try:
        if SHIFT_TARGET is None:
            evals, evecs = spla.eigsh(
                Ssp, k=N_TARGET, M=Msp,
                which='SA', tol=1e-6, maxiter=2000
            )
        else:
            evals, evecs = spla.eigsh(
                Ssp, k=N_TARGET, M=Msp,
                sigma=SHIFT_TARGET, which='LM',
                tol=1e-6, maxiter=2000
            )
        print(f"  eigsh done in {time.time()-t0:.1f}s")
        print(f"  Eigenvalues: {evals}")
    except Exception as e:
        print(f"  eigsh failed: {e}")
        return

    # Sort by ascending eigenvalue
    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]

    print(f"\n=== Mode analysis ===")
    sigma_cu_only = CoefficientFunction(
        [sigma_Cu if mat == "conductor" else 0.0 for mat in mesh.GetMaterials()]
    )
    sigma_vac_only = CoefficientFunction(
        [0.0 if mat == "conductor" else SIGMA_VAC for mat in mesh.GetMaterials()]
    )
    A_ext = CoefficientFunction((-y/2, x/2, 0))

    for k in range(len(evals)):
        lam = float(evals[k])
        gf = GridFunction(fes)
        gf.vec.FV().NumPy()[:] = evecs[:, k]
        m_cu = float(Integrate(sigma_cu_only * gf * gf * dx, mesh))
        m_vac = float(Integrate(sigma_vac_only * gf * gf * dx, mesh))
        m_tot = m_cu + m_vac
        cu_frac = m_cu / m_tot if m_tot > 0 else 0.0
        if m_tot > 1e-30:
            gf.vec.FV().NumPy()[:] /= m_tot**0.5
            beta = float(Integrate(sigma_cf * gf * A_ext * dx, mesh))
            beta_sq = beta**2
        else:
            beta_sq = 0.0
        tau_us = 1.0 / lam * 1e6 if abs(lam) > 1e-30 else float('inf')
        print(f"  mode {k:2d}: lam={lam:+.5e}, tau={tau_us:+.4f} us, "
              f"cu_frac={cu_frac:.4f}, beta2={beta_sq:.3e}")

    print()
    print("=" * 64)
    print("References:")
    print("  Interior Dirichlet TE_z(1,1,0) analytical: 25.5 us")
    print("  Phase 7-11 BEM (parity-restricted):        19.77 us")
    print("=" * 64)


if __name__ == "__main__":
    main()
