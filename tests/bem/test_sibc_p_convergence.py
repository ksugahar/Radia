"""SIBC p-convergence sanity: order=1 vs order=2 on a conducting sphere.

Reference: Smythe, "Static and Dynamic Electricity" 3rd ed, Section 11.

For a conducting sphere of radius R in a uniform external field
H = H0 * z_hat with surface impedance Z_s applied on r=R, the surface
tangential current is::

    J_s(theta) = (3/2) * j*omega*B0*R / (j*omega*mu_0*R + Z_s) * sin(theta)
    H_t_rms    = max(|J_s|) * sqrt(2/3)

This test asserts that H1 order=2 reduces the error vs analytical by
at least 100x compared to order=1 on the same mesh.  Currently the
in-tree C++ assembler is P1-only (`AssembleSLDL` in rad_bem_galerkin.cpp
documents "P1-H1 surface mesh"), so this test runs through the
ngsolve.bem path (use_intree_bem=False).  Once the C++ assembler is
extended to P2 H1, set DEFAULT_USE_INTREE=True below to gate on the
in-tree path instead.

Phase 0 baseline (2026-05-01, ngsolve.bem path):
    maxh/R   order    ndof    rel.err
    R/3       1       128     2.5e-3
    R/3       2       506     7.7e-7   (3300x better)
    R/5       1       363     8.7e-4
    R/5       2      1446     8.7e-8   (10000x better)
"""
import math
import os
import sys

import numpy as np
import pytest


MU_0 = 4e-7 * math.pi


def _have_ngsolve_bem():
    try:
        from ngsolve.bem import LaplaceDL, LaplaceSL  # noqa: F401
        return True
    except Exception:
        return False


@pytest.mark.skipif(not _have_ngsolve_bem(),
                     reason="ngsolve.bem not available")
def test_p_convergence_sphere():
    """p=2 must give >= 100x lower error than p=1 on the same mesh."""
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
    for order in (1, 2):
        solver = ScalarBIESIBCSolver(mesh, order=order, use_intree_bem=False,
                                      assemble_dense=True)
        # phi_inc = -H0 * z, projected onto the FES at the BND
        from ngsolve import GridFunction, z
        gf = GridFunction(solver.fes)
        gf.Set(-H0 * z, definedon=mesh.Boundaries(".*"))
        phi_inc = np.asarray(gf.vec.FV().NumPy(), dtype=complex).copy()

        sol = solver.solve(phi_inc, Z_s=Z_s, omega=omega)
        H_t_rms = float(sol["H_t_rms"])
        errs[order] = abs(H_t_rms - H_t_rms_ana) / H_t_rms_ana

    assert errs[1] > 0, f"p=1 error vanished unexpectedly: {errs[1]}"
    assert errs[2] > 0, f"p=2 error vanished unexpectedly: {errs[2]}"
    speedup = errs[1] / errs[2]
    print(f"\np=1 rel.err = {errs[1]:.3e}")
    print(f"p=2 rel.err = {errs[2]:.3e}")
    print(f"p=2/p=1 ratio (smaller is better) = {errs[2]/errs[1]:.3e}")
    print(f"speedup factor = {speedup:.0f}x")
    assert speedup > 100, (
        f"p=2 should give at least 100x lower error than p=1 on the "
        f"same mesh (got {speedup:.1f}x). p=1 err={errs[1]:.3e}, "
        f"p=2 err={errs[2]:.3e}")


if __name__ == "__main__":
    test_p_convergence_sphere()
    print("OK: p-convergence regression test passed.")
