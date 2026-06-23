# -*- coding: utf-8 -*-
r"""Regression test: CURVED P2 (isoparametric, mesh.Curve(2)) axisymmetric
magnetostatics on the uniformly magnetized sphere.

Selected implicitly by solve_axi_magnetostatic(..., order=2) on a curved mesh: the
P2 triangle samples its node positions from the Curve(2)-d ElementTransformation, so
the sphere's circular boundary is captured to 2nd order instead of by straight
chords.  This pins the CORRECT behaviour so a silent regression (curving stops
engaging, or the geometry stops being used in the 1/r integration) is caught.

WHERE curving shows up -- evaluate the VOLUME (extensive), not the average.
The magnet half-disk has exact volume integral  D = INT r dr dz = 2/3, and the
interior B_z is exactly uniform (= B_EXACT), so the exact total axial flux is
N = INT r B_z dr dz = B_EXACT * 2/3 (proportional to the magnet's TOTAL MAGNETIC
MOMENT, which sets the entire external field).  The volume-AVERAGE is avg = N/D.

  * Volume D:    straight chords UNDER-integrate the disk by ~-3.9% at h=0.4;
                 Curve(2) recovers it to -0.025% (machine-precision circle, ~150x).
  * Total flux N (~ total moment): straight -5.07%, Curve(2) -1.07% at h=0.4 --
                 curving cuts the error ~5x, because N keeps the volume error.
  * Average N/D: straight -1.22%, Curve(2) -1.05% -- the -3.9% volume error CANCELS
                 in the ratio (avg_err ~= N_err - D_err), so the average looks almost
                 unchanged.  That ~0.18% is NOT a sign curving is impotent; it is the
                 demagnetizing factor N=1/3 being variationally stationary at the
                 sphere (first-order faceting cancels in the average).  The average is
                 simply the WRONG observable -- evaluate the volume/moment instead.

The default quadrature already integrates the curved 1/r form adequately (a
bonus_intorder 0->12 sweep moves N by <0.01%), so the gain is geometry-limited.

Requires the rebuilt axifem extension (Build.ps1 -AxiFemmOnly); loads src/radia and
the radia-mcp package source.
"""
import math
import os
import sys

_HERE = os.path.dirname(__file__)
_SRC = os.path.abspath(os.path.join(_HERE, "..", "..", "src"))
_PKG = os.path.abspath(os.path.join(_HERE, "..", "..", "packages", "radia-mcp", "src"))
for _p in (_SRC, _PKG):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from ngsolve import (CoefficientFunction as CF, Mesh, grad, Integrate, TaskManager,
                     x as r_cf)
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import solve_axi_magnetostatic

MU0, A, MUR, HC = 4e-7 * math.pi, 1.0, 2.0, 3.0e5
B_EXACT = 2 * MU0 * MUR * HC / (MUR + 2)        # uniform interior field, 0.376991 T
D_EXACT = 2.0 / 3.0                              # INT_{half disk r>=0, radius 1} r dr dz
N_EXACT = B_EXACT * D_EXACT                       # total axial flux ~ total moment
RNEAR, RFAR = 3.0, 40.0
H_COARSE = 0.4                                    # coarse: curving benefit is largest here


def _build(h):
    """Graded sphere mesh: magnet (half-disk) AND the near-air shell (radius RNEAR)
    both at maxh=h so the shared boundary is resolved on BOTH sides; far air coarse.
    (Resolving only the magnet, air left coarse, would let the coarse air side of the
    same boundary dominate and mask curving.)"""
    disk = WorkPlane().Circle(0, 0, A).Face()
    half = MoveTo(0, -RFAR).Rectangle(RFAR, 2 * RFAR).Face()
    mag = disk * half
    mag.faces.name = "magnet"; mag.maxh = h; mag.edges.Min(X).name = "axis"
    near = (WorkPlane().Circle(0, 0, RNEAR).Face() * half) - mag
    near.faces.name = "air"; near.maxh = h; near.edges.Min(X).name = "axis"
    far = (MoveTo(0, -RFAR).Rectangle(RFAR, 2 * RFAR).Face()) - WorkPlane().Circle(0, 0, RNEAR).Face() * half
    far.faces.name = "air"; far.edges.Min(X).name = "axis"
    far.edges.Max(X).name = "outer"; far.edges.Max(Y).name = "outer"; far.edges.Min(Y).name = "outer"
    return Mesh(OCCGeometry(Glue([far, near, mag]), dim=2).GenerateMesh(maxh=4.0))


def _magnet_volume(mesh):
    """D = INT r dr dz over the magnet (= magnet volume / 2pi)."""
    return Integrate(r_cf, mesh, definedon=mesh.Materials("magnet"))


def _vol_flux_avg(mesh):
    """Solve order-2 magnetostatics, return (D, N, avg) for the magnet:
    D = INT r dr dz (volume), N = INT r B_z dr dz (total flux ~ moment), avg = N/D."""
    nu = CF(1.0 / (MU0 * mesh.MaterialCF({"magnet": MUR}, default=1.0)))
    with TaskManager():
        gfu = solve_axi_magnetostatic(mesh, nu, magnets={"magnet": (HC, 90.0)}, order=2)
    Bz = grad(gfu)[0] + gfu / r_cf
    mag = mesh.Materials("magnet")
    D = Integrate(r_cf, mesh, definedon=mag)
    N = Integrate(r_cf * Bz, mesh, definedon=mag)
    return D, N, N / D


def test_curved_p2_geometry_engages():
    """mesh.Curve(2) moves the P2 nodes onto the true circle: the magnet volume
    INT r dr dz (exact 2/3) is recovered from a ~3.9% straight-chord error to
    <0.1%.  If curving silently stops engaging, straight and curved volumes coincide
    and this fails."""
    err_s = abs(_magnet_volume(_build(H_COARSE)) - D_EXACT) / D_EXACT * 100
    mesh_c = _build(H_COARSE); mesh_c.Curve(2)
    err_c = abs(_magnet_volume(mesh_c) - D_EXACT) / D_EXACT * 100
    print(f"  volume err: straight {err_s:.3f}%   Curve(2) {err_c:.4f}%")
    assert err_s > 1.0, f"straight chords should under-integrate the disk (got {err_s:.3f}%)"
    assert err_c < 0.1, f"Curve(2) must recover the true circle (got {err_c:.4f}%)"
    assert err_c < err_s / 10.0, "curving must reduce the geometry error by >=10x"
    print("[OK] Curve(2) captures the circular boundary to machine precision")


def test_curved_p2_total_flux_tracks_volume():
    """Evaluate the VOLUME, not the average.  The total axial flux N = INT r B_z dr dz
    (~ total magnetic moment) keeps the volume error, so curving cuts its error ~5x
    (straight -5.07%, Curve(2) -1.07% at h=0.4).  The volume-AVERAGE N/D cancels the
    volume error and barely moves (~0.18%) -- asserted here only as the CONTRAST that
    exposes the mechanism (avg_err ~= N_err - D_err), not as a curving failure."""
    Ds, Ns, avg_s = _vol_flux_avg(_build(H_COARSE))
    mesh_c = _build(H_COARSE); mesh_c.Curve(2)
    Dc, Nc, avg_c = _vol_flux_avg(mesh_c)

    eD_s = (Ds - D_EXACT) / D_EXACT * 100
    eD_c = (Dc - D_EXACT) / D_EXACT * 100
    eN_s = (Ns - N_EXACT) / N_EXACT * 100
    eN_c = (Nc - N_EXACT) / N_EXACT * 100
    eA_s = (avg_s - B_EXACT) / B_EXACT * 100
    eA_c = (avg_c - B_EXACT) / B_EXACT * 100
    print(f"  D(vol)  err%: straight {eD_s:+8.4f}  Curve(2) {eD_c:+8.4f}")
    print(f"  N(flux) err%: straight {eN_s:+8.4f}  Curve(2) {eN_c:+8.4f}   (ratio {eN_c/eN_s:.3f})")
    print(f"  avg     err%: straight {eA_s:+8.4f}  Curve(2) {eA_c:+8.4f}   (delta {abs(eA_c-eA_s):.4f})")

    # the dramatic effect: curving cuts the total-flux / total-moment error >=2x
    assert abs(eN_c) < 0.5 * abs(eN_s), \
        f"Curve(2) must cut total-flux error >=2x (straight {eN_s:+.3f}%, curved {eN_c:+.3f}%)"
    # mechanism: total-flux error = field-shape error (the average) + volume error
    assert abs(eN_s - (eA_s + eD_s)) < 0.15 and abs(eN_c - (eA_c + eD_c)) < 0.05, \
        "N_err ~= avg_err + D_err should hold (the volume error is what curving removes)"
    # contrast: curving moves the total flux a lot but the average almost not at all
    assert abs(eN_s - eN_c) > 2.0 and abs(eA_s - eA_c) < 0.5, \
        "the volume/flux must show curving dramatically while the average cancels it"
    print("[OK] evaluating the volume/total-flux exposes curving (the average hides it)")


if __name__ == "__main__":
    test_curved_p2_geometry_engages()
    test_curved_p2_total_flux_tracks_volume()
    print("\nCurved P2 magnetostatic tests passed.")
