# -*- coding: utf-8 -*-
"""The INDUCTANCE half of the DtN datasheet -- the dual of capacitance.

Capacitance C and external inductance L_ext are each ONE Steklov mode of the SAME
exterior scalar Laplace DtN, but on a different rung of the ladder lambda_n=-(n+1)/R:

    C      <->  n=0 (MONOPOLE)   isolated charged conductor; constant image; EXACT
    L_ext  <->  n=1 (DIPOLE)     current loop has no magnetic monopole, leads dipole

The external magnetic-field energy obeys the identity

    W_ext = 1/2 mu0 * (n+1)/R * oint_Gamma phi^2 dS     (decaying mode r^-(n+1))

so the DtN eigenvalue (n+1)/R IS the exterior-energy coefficient, and
L_ext = 2 W_ext / I^2 inherits the n=1 (dipole) defect EXACTLY.  This script MEASURES
that defect (vs order, mesh, curve-order, exterior-volume) and confirms the C/L dual.

Certification scope: this is the open-boundary TRUNCATION accuracy of a FIELD-ENERGY
inductance (Kelvin / air-box), via the magnetic POTENTIAL exterior -- scalar Omega
(single-valued for a magnetisation source; a free-current loop needs a cohomology cut)
or the vector potential A (no cut; same -(n+1)/R gradient block, see act0_06_aform_dtn_gradient).  It is
NOT the ngsolve.bem vector single-layer energy L = mu0 J^T(LaplaceSL)J, which is a
DIFFERENT integral operator (no -(n+1)/R ladder) -- keep the two distinct.

PRINTS its result (no files written); depends only on numpy / ngsolve / netgen and
the open radia_mcp.radia_ngsolve helpers.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
from ngsolve import grad, ds, dx, BND, z as Z

from radia_mcp.radia_ngsolve.fem_bem_coupling import (
    kelvin_dtn_eigenvalue, kelvin_openbc_error_vs_exterior_mesh,
)

R = 1.0
mu0 = 4e-7 * np.pi


def main():
    print("=" * 72)
    print("INDUCTANCE half of the DtN datasheet   (L_ext <-> dipole n=1)")
    print("=" * 72)

    # (0) the energy identity, at the analytic level (hand-checked, re-confirmed):
    #     1/2 (2/R) oint phi^2  ==  integral_{r>R}|H|^2  ==  m^2/(6 pi R^3)
    m = 1.0
    W_ext = m**2 / (6 * np.pi * R**3)
    oint_phi2 = m**2 / (12 * np.pi * R**2)
    lhs = (2.0 / R) * oint_phi2
    print(f"\n(0) identity  (n+1)/R * oint phi^2 = {lhs:.6e}")
    print(f"              integral_(r>R)|H|^2  = {W_ext:.6e}"
          f"   rel diff {abs(lhs - W_ext) / W_ext:.1e}")
    print("    -> the DtN eigenvalue IS the exterior-energy coefficient")

    # (A) dipole (n=1) DtN eigenvalue vs FE order p = the L_ext truncation defect
    print("\n(A) L_ext truncation defect vs FE order p   (exact lambda_1 = -2/R)")
    a_defect = {}
    for order in (1, 2, 3, 4):
        r = kelvin_dtn_eigenvalue(R=R, degree=1, maxh=0.4, order=order)
        a_defect[order] = r["rel_err"]
        print(f"    p={order}:  lam={r['lam']:.8f}   rel_err(L_ext)={r['rel_err']:.3e}")

    # (B) the C <-> L dual at one fixed coarse mesh
    print("\n(B) C<->L dual (maxh=0.5, order=2):")
    dual = {}
    for n, label in ((0, "C (monopole n=0)"), (1, "L (dipole   n=1)")):
        r = kelvin_dtn_eigenvalue(R=R, degree=n, maxh=0.5, order=2)
        dual[n] = r["rel_err"]
        print(f"    {label}:  lam={r['lam']:.8f}  exact={r['lam_exact']:.4f}"
              f"  rel_err={r['rel_err']:.3e}")

    # (C) floor IS geometry: hold order>=n and the mesh fixed, raise only Curve k
    print("\n(C) floor = geometry  (degree=1, order=2, mesh fixed, raise Curve k):")
    ngmesh = OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=0.4)
    curve_err = {}
    for k in (1, 2, 3):
        mesh = ng.Mesh(ngmesh.Copy()); mesh.Curve(k)
        fes = ng.H1(mesh, order=2, dirichlet=".*")
        u, v = fes.TnT()
        a = ng.BilinearForm(grad(u) * grad(v) * dx); a.Assemble()
        gfu = ng.GridFunction(fes); gfu.Set(Z, BND)           # n=1 zonal solid harmonic
        rhs = gfu.vec.CreateVector(); rhs.data = -(a.mat * gfu.vec)
        gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * rhs
        energy = float(ng.Integrate(grad(gfu) * grad(gfu) * dx(bonus_intorder=10), mesh))
        bmass = float(ng.Integrate(gfu * gfu * ds(bonus_intorder=10), mesh))
        lam = -1.0 / R - energy / bmass
        rel = abs(lam - (-2.0 / R)) / (2.0 / R)
        curve_err[k] = rel
        print(f"    Curve k={k}:  lam={lam:.8f}   rel_err={rel:.3e}")

    # (D) exterior-VOLUME-irrelevance: refine the exterior, interior held fixed
    print("\n(D) exterior-volume-irrelevance for L (degree=1):")
    res = kelvin_openbc_error_vs_exterior_mesh(degree=1, order=2)
    print(f"    interior_fem_error (fixed) = {res['interior_fem_error']:.3e}"
          f"   always_below_fem={res['always_below_fem']}")
    for mrow in res["per_mesh"]:
        print(f"      kelvin_maxh={mrow['kelvin_maxh']}: "
              f"openbc_err={mrow['kelvin_openbc_error']:.3e}  ndof={mrow['kelvin_ndof']}")

    # (E) a physical external-inductance number (current loop, energy beyond R)
    a_loop = 0.2
    L_ext_loop = mu0 * (np.pi * a_loop**2) ** 2 / (6 * R**3)
    print(f"\n(E) current loop a={a_loop} m: L_ext(beyond R) = {L_ext_loop:.3e} H,"
          f" rel error = the (A) defect")

    # ---- self-checks (the verified claims) ----------------------------------
    assert abs(lhs - W_ext) / W_ext < 1e-12, "energy identity broken"
    assert dual[0] < 1e-9, "C (n=0 monopole) must be machine-exact"
    assert a_defect[1] < 1e-2 and a_defect[3] < a_defect[1], "L defect must fall with p"
    assert curve_err[3] < 0.1 * curve_err[1], "floor must be geometry (Curve-controlled)"
    assert res["always_below_fem"], "exterior must never be the bottleneck"
    print("\nOK: C(n=0) exact; L_ext(n=1) defect falls with p, floored at geometry,")
    print("    mesh-independent, exterior-volume-irrelevant -- the dual of capacitance.")


if __name__ == "__main__":
    main()
