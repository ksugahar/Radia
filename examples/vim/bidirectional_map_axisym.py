"""Axisymmetric flux-coordinate bidirectional map (Tampere -> Radia, Phase 2).

The axisymmetric (helicity-0) generalisation of the 2D potential-coordinate map
(bidirectional_map_2d.py).  In the (r, z) meridian the two natural coordinates of
a CURRENT-FREE axisymmetric field are

    psi  = r * A_phi    -- the Stokes STREAM function  (B-flux coordinate)
    Phi  -- the magnetic SCALAR potential              (B = grad Phi)

a conjugate pair under the AXISYMMETRIC (1/r-weighted) Cauchy-Riemann relations

    dPhi/dr = -(1/r) dpsi/dz  ( = B_r ),   dPhi/dz = +(1/r) dpsi/dr  ( = B_z ).

FORWARD (geometry -> potentials): given the meridian geometry + boundary data,
solve psi (Grad-Shafranov weak form, current-free) and Phi (axisymmetric Laplace)
-> the potential coordinates of every point.
INVERSE (potentials -> geometry): r(psi, Phi), z(psi, Phi) -- "redraw the geometry
in potential coordinates".

KEY (verified sympy 2026-06-11, now in FEM):
  * Jacobian  det[grad psi; grad Phi] = r |grad Phi|^2 = r|B|^2 (meridian); the map
    is a valid chart exactly where r|B|^2 != 0 (degeneracy = the axis / field nulls).
  * FINDING -- the r-metric breaks conformality: the 2D map is CONFORMAL (|grad A| =
    |grad phi| -> the inverse x,y are HARMONIC in (A, phi), a plain Laplace solve).
    The axisymmetric map is only r-CONFORMAL (|grad psi| = r|grad Phi|, angle- but
    not scale-preserving) -> the inverse is NOT a plain harmonic image solve; it is
    a pointwise Newton inversion of the (invertible, Jacobian r|B|^2) forward map.
    (For our pair psi = 2z(z^2 - Phi) is a CUBIC in z -- genuinely non-harmonic.)

Validation pair = the axisymmetric GRADIENT (l = 2 zonal solid-harmonic) field, the
accelerator-relevant field shape (B_r = -r, B_z = 2z, current-free, Laplace Phi = 0):
    psi = r^2 z ,    Phi = z^2 - r^2/2
on the meridian box r in [0.3, 1.0], z in [-0.4, 0.4] (away from the axis null).

Run:  python examples/vim/bidirectional_map_axisym.py
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../src/radia"))

from ngsolve import (                                      # noqa: E402
    Mesh, H1, GridFunction, grad, InnerProduct, dx, sqrt,
    x, y, BilinearForm, LinearForm, TaskManager, Integrate,
)
from netgen.geom2d import SplineGeometry                   # noqa: E402

# meridian box (NGSolve 2D coords: x = r, y = z), away from the axis r = 0
R0, R1, Z0, Z1 = 0.3, 1.0, -0.4, 0.4

# analytic axisymmetric gradient (l = 2 zonal) field, current-free
PSI_exact = x * x * y           # Stokes stream function  r^2 z
PHI_exact = y * y - x * x / 2   # magnetic scalar potential  z^2 - r^2/2


def _solve(maxh, order, bc_cf, weight):
    """Harmonic extension of bc_cf under the weighted operator int weight*grad.grad."""
    geo = SplineGeometry()
    geo.AddRectangle((R0, Z0), (R1, Z1), bc="bnd")
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    fes = H1(mesh, order=order, dirichlet="bnd")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += weight * InnerProduct(grad(u), grad(v)) * dx
    f = LinearForm(fes)
    a.Assemble()
    f.Assemble()
    inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    gf = GridFunction(fes)
    gf.Set(bc_cf, definedon=mesh.Boundaries("bnd"))
    r = gf.vec.CreateVector()
    r.data = f.vec - a.mat * gf.vec
    gf.vec.data += inv * r
    return mesh, gf


def solve_forward(maxh=0.04, order=3):
    # psi: Grad-Shafranov weak form  int (1/r) grad psi . grad v   (current-free)
    mesh, gfPsi = _solve(maxh, order, PSI_exact, 1.0 / x)
    # Phi: axisymmetric Laplace      int r grad Phi . grad v
    _, gfPhi = _solve(maxh, order, PHI_exact, x)
    return mesh, gfPsi, gfPhi


def forward_checks(mesh, gfPsi, gfPhi):
    errPsi = sqrt(Integrate((gfPsi - PSI_exact) ** 2, mesh) /
                  Integrate(PSI_exact ** 2, mesh))
    errPhi = sqrt(Integrate((gfPhi - PHI_exact) ** 2, mesh) /
                  Integrate(PHI_exact ** 2, mesh))
    gPsi, gPhi = grad(gfPsi), grad(gfPhi)
    gPhi2 = InnerProduct(gPhi, gPhi)
    # axisymmetric Cauchy-Riemann: dPhi/dr + (1/r)dpsi/dz = 0, dPhi/dz - (1/r)dpsi/dr = 0
    c1 = gPhi[0] + gPsi[1] / x
    c2 = gPhi[1] - gPsi[0] / x
    conj = sqrt(Integrate(c1 ** 2 + c2 ** 2, mesh) / Integrate(gPhi2, mesh))
    # Jacobian det[grad psi; grad Phi] = r |grad Phi|^2  ( = r|B|^2 )
    detJ = gPsi[0] * gPhi[1] - gPsi[1] * gPhi[0]
    jac_rel = sqrt(Integrate((detJ - x * gPhi2) ** 2, mesh) /
                   Integrate((x * gPhi2) ** 2, mesh))
    return {"errPsi": float(errPsi), "errPhi": float(errPhi),
            "conj": float(conj), "jac_rel": float(jac_rel)}


def _fwd(r, z):
    return r * r * z, z * z - r * r / 2


def _inv_newton(psi, phi, r_guess, z_guess, iters=40):
    """Pointwise inverse: Newton on F(r,z) = (psi, phi).  J = [[2rz, r^2],[-r, 2z]]."""
    r, z = r_guess, z_guess
    for _ in range(iters):
        fr, fz = _fwd(r, z)
        gr, gz = fr - psi, fz - phi
        j11, j12, j21, j22 = 2 * r * z, r * r, -r, 2 * z
        det = j11 * j22 - j12 * j21          # = r(r^2 + 4z^2) = r|B|^2 != 0
        dr = (j22 * gr - j12 * gz) / det
        dz = (-j21 * gr + j11 * gz) / det
        r, z = r - dr, z - dz
        if math.hypot(dr, dz) < 1e-14:
            break
    return r, z


def roundtrip_error(mesh, gfPsi, gfPhi, samples=None):
    """Physical p -> (psi,Phi) via FORWARD FEM -> inverse Newton -> physical; max err."""
    if samples is None:
        samples = [(0.5, 0.2), (0.7, -0.3), (0.9, 0.1), (0.4, -0.15), (0.8, 0.35)]
    rg, zg = 0.5 * (R0 + R1), 0.5 * (Z0 + Z1)   # fixed guess = domain centre
    worst = 0.0
    for (r, z) in samples:
        mp = mesh(r, z)
        psi, phi = float(gfPsi(mp)), float(gfPhi(mp))
        rr, zz = _inv_newton(psi, phi, rg, zg)
        worst = max(worst, math.hypot(rr - r, zz - z))
    return worst


def main():
    print("=== axisymmetric flux-coordinate bidirectional map (Tampere, Phase 2) ===")
    with TaskManager():
        mesh, gfPsi, gfPhi = solve_forward()
        chk = forward_checks(mesh, gfPsi, gfPhi)
        rt = roundtrip_error(mesh, gfPsi, gfPhi)

    print("\n[forward]  geometry -> (psi = r A_phi, Phi)   (GS + axisym Laplace)")
    print(f"   L2 err psi         = {chk['errPsi']:.2e}")
    print(f"   L2 err Phi         = {chk['errPhi']:.2e}")
    print(f"   axisym CR conj.    = {chk['conj']:.2e}   "
          f"(dPhi/dr+(1/r)dpsi/dz, dPhi/dz-(1/r)dpsi/dr ~ 0)")
    print(f"   Jacobian rel-err   = {chk['jac_rel']:.2e}   (det = r |grad Phi|^2)")
    print("\n[inverse]  (psi, Phi) -> geometry  (pointwise Newton; r-conformal, NOT")
    print("            a plain harmonic image solve -- the r-metric breaks conformality)")
    print(f"   round-trip max err = {rt:.2e}")

    ok = (chk["errPsi"] < 1e-3 and chk["errPhi"] < 1e-3 and
          chk["conj"] < 1e-2 and chk["jac_rel"] < 1e-2 and rt < 1e-2)
    print(f"\n=== PASS: {ok} ===")
    return chk, rt, ok


if __name__ == "__main__":
    *_, ok = main()
    raise SystemExit(0 if ok else 1)
