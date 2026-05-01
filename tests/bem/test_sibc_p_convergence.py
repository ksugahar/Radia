"""SIBC p-convergence sanity: order=1 vs order=2 on a conducting sphere.

Reference: Smythe, "Static and Dynamic Electricity" 3rd ed, Section 11.

For a conducting sphere of radius R in a uniform external field
H = H0 * z_hat with surface impedance Z_s applied on r=R, the surface
tangential current is::

    J_s(theta) = (3/2) * j*omega*B0*R / (j*omega*mu_0*R + Z_s) * sin(theta)
    H_t_rms    = max(|J_s|) * sqrt(2/3)

This test asserts that H1 order=2 reduces the error vs analytical by
at least 100x compared to order=1 on the same mesh.  Both orders run
through the in-tree path (use_intree_bem=True): order=1 uses the C++
P1 Galerkin assembler (`AssembleSLDL`), order=2 uses the Python
Lagrange-P2 assembler (`assemble_SL/DL_dense_curved_p2`) with the
solver bypassing NGSolve's H1 hierarchical FES so the matrices are
self-consistent in the Lagrange basis.

Phase 0 / Phase A baseline (2026-05-01):
    maxh/R   order    ndof    rel.err  (in-tree)   notes
    R/3       1       128     ~2.5e-3              C++ P1 assembler
    R/3       2       506     ~2.1e-6              Python Lagrange P2
                                                   (~1000x improvement;
                                                   ngsolve.bem oracle
                                                   gives 7.7e-7).
"""
import math
import numpy as np
import pytest


MU_0 = 4e-7 * math.pi


def _have_netgen_occ():
    try:
        from netgen.occ import Sphere  # noqa: F401
        return True
    except Exception:
        return False


@pytest.mark.skipif(not _have_netgen_occ(),
                     reason="netgen.occ not available")
def test_p_convergence_sphere_intree():
    """In-tree path: p=2 must give >= 100x lower error than p=1."""
    from ngsolve import Mesh
    from netgen.occ import Sphere, Pnt, OCCGeometry, Glue
    from radia.bem_sibc_solver import ScalarBIESIBCSolver

    R = 0.030
    sigma = 5.8e7
    mu_r = 1.0
    frequency = 50000
    B0 = 1e-3
    omega = 2 * math.pi * frequency
    mu_eff = MU_0 * mu_r
    Z_s = (1 + 1j) * math.sqrt(omega * mu_eff / (2.0 * sigma))
    beta = 1j * omega * MU_0 * R
    J_s_max = abs((3.0 / 2.0) * 1j * omega * B0 * R / (beta + Z_s))
    H_t_rms_ana = J_s_max * math.sqrt(2.0 / 3.0)

    sph = Sphere(Pnt(0, 0, 0), R)
    for f in sph.faces:
        f.name = "sibc"
    geo = OCCGeometry(Glue(sph.faces))
    mesh = Mesh(geo.GenerateMesh(maxh=R / 3))
    mesh.Curve(2)
    H0 = B0 / MU_0

    errs = {}

    # order = 1 : in-tree P1 path (uses C++ AssembleSLDL).  phi_inc must be
    # an ndarray of length n_v (vertex DOFs), so we sample analytical phi
    # at each surface vertex.
    solver1 = ScalarBIESIBCSolver(mesh, order=1, use_intree_bem=True,
                                    assemble_dense=True)
    phi_inc1 = np.zeros(solver1.ndof, dtype=complex)
    for i in range(mesh.nv):
        phi_inc1[i] = -H0 * mesh.vertices[i].point[2]
    sol1 = solver1.solve(phi_inc1, Z_s=Z_s, omega=omega)
    errs[1] = abs(sol1["H_t_rms"] - H_t_rms_ana) / H_t_rms_ana

    # order = 2 : in-tree Lagrange-P2 path.  Caller samples phi at each
    # Lagrange node via solver.dof_coords (vertices + edge mid-points).
    solver2 = ScalarBIESIBCSolver(mesh, order=2, use_intree_bem=True,
                                    assemble_dense=True)
    phi_inc2 = (-H0 * solver2.dof_coords[:, 2]).astype(complex)
    sol2 = solver2.solve(phi_inc2, Z_s=Z_s, omega=omega)
    errs[2] = abs(sol2["H_t_rms"] - H_t_rms_ana) / H_t_rms_ana

    print(f"\np=1 rel.err = {errs[1]:.3e}  (ndof={solver1.ndof})")
    print(f"p=2 rel.err = {errs[2]:.3e}  (ndof={solver2.ndof})")
    speedup = errs[1] / errs[2]
    print(f"p=2 vs p=1 error ratio = {errs[2]/errs[1]:.3e} ({speedup:.0f}x lower err)")

    assert errs[1] > 0 and errs[2] > 0
    assert speedup > 100, (
        f"p=2 should give at least 100x lower error than p=1 on the "
        f"same mesh (got {speedup:.1f}x). p=1 err={errs[1]:.3e}, "
        f"p=2 err={errs[2]:.3e}")


if __name__ == "__main__":
    test_p_convergence_sphere_intree()
    print("OK: p-convergence regression test passed.")
