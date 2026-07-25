"""Golden: A-Phi (A-V) eddy formulation + Periodic Kelvin WORKS at p=2 --
and beats the plain A-method, which hits a p-independent accuracy floor.

Problem: conducting non-magnetic sphere (radius a = 0.4 m, sigma = 5.8e7 S/m)
in a uniform harmonic field B0 z_hat at a/delta = 2 (moderate skin effect),
open boundary by the two-sphere Periodic Kelvin transformation.  Classical
analytic induced dipole moment (Smythe):

    m_z = -2 pi a^3 (B0/mu0) [ 1 - (3/x) coth x + 3/x^2 ],  x = (1+j) a/delta

(low-frequency limit -x^2/15 and high-frequency PEC limit -1 both verified).

Reduced split A = A_s + A_r with A_s = (B0/2)(-y, x, 0).  Since nu = nu0
everywhere, curl(nu0 curl A_s) = 0 identically, so the source enters ONLY
through the conductor mass term and NO background is needed in the Kelvin
exterior -- the reaction field decays like a dipole, which Kelvin represents
exactly.  (This is the decaying-source case of
test_kelvin_exterior_source_routes.py; the 2/3 loss there applies to
non-decaying backgrounds only.)

Both lanes share one mesh, the verified_recipe elements (nograds=True ->
Periodic, gauge reg on non-Kelvin materials, bonus_intorder=4 on curl-curl),
and a direct solver:

    A*    J = -s sigma (A_s + A_r)                      (A-method, no Phi)
    A-Phi J = -s sigma (A_s + A_r + grad W), W = V/s    (mixed HCurl x H1)
          symmetric form: nu curlcurl + s sigma (A_r + grad W).(A' + grad W')

MEASURED (LAB, 2026-07-25), rel = |m_fem - m_ana| / |m_ana|:

    p-sweep (ne=10000):
    p    A* (A-method)   A-Phi        FES verify (slaved / ratio)
    1    2.889%          2.889%       1659 / 1.000000   (lanes identical)
    2    0.473%          0.053%       3871 / 1.000000
    3    0.442%          0.005%       7189 / 1.000000

    h-sweep at p=2:
    ne       A* (A-method)   A-Phi     factor
    10000    0.473%          0.053%      9x
    28022    0.243%          0.011%     22x
    67850    0.141%          0.004%     35x

    gauge-reg sensitivity at p=2: 1e-5..1e-10 -> A* 0.4734% / A-Phi 0.0529%
    (identical to 4 digits: the defect is NOT the regularization)

Conclusions locked here:

1. A-Phi + Periodic Kelvin at p=2 reproduces the analytic moment to 0.05%:
   the answer to "does A-Phi Kelvin work at p=2" is YES.
2. The plain A-method's error is p-SATURATED (p=3 does not improve on p=2)
   and decays only slowly under h-refinement (~O(h): 0.473 -> 0.243 ->
   0.141%), while A-Phi converges fast in both p and h.  With nograds=True
   -- which the Periodic-Kelvin recipe requires -- the gradient test
   functions that would enforce discrete charge conservation are removed
   from the space; the explicit Phi (W) block restores them.  For this
   axisymmetric problem the CONTINUOUS V is exactly zero, so the entire
   A*/A-Phi gap is discrete charge-conservation error: with A* you can
   only pay for accuracy in h, with A-Phi p-refinement pays.
3. At p=1 the two lanes coincide (the p=1 W coupling is inert here); the
   formulation choice starts to matter at p=2 exactly where the A-method
   stops improving.

The p=3 row is recorded above but not executed here (runtime); the golden
locks p=1 and p=2.

Reference: radia_mcp kelvin_transformation(topic="verified_recipe") for the
five mandatory elements; AV_formulation_plan.md for the A-V weak form.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from netgen.meshing import IdentificationType
from netgen.occ import Glue, OCCGeometry, Pnt, Sphere, Vertex
from ngsolve import (
    BilinearForm, CoefficientFunction as CF, GridFunction, H1, HCurl,
    InnerProduct, Integrate, LinearForm, Mesh, Periodic, TaskManager,
    curl, dx, grad, x, y, z,
)

MU_0 = 4e-7 * math.pi
NU_0 = 1.0 / MU_0

A_C = 0.4            # conductor radius [m]
R_K = 1.0            # Kelvin radius [m]
OFFSET = (3.0, 0.0, 0.0)
SIGMA = 5.8e7        # [S/m]
B_0 = 1.0            # [T]
A_OVER_DELTA = 2.0   # moderate skin effect
MAXH = 0.16
MAXH_COND = 0.07

_CACHE: dict = {}


def _analytic():
    delta = A_C / A_OVER_DELTA
    omega = 2.0 / (MU_0 * SIGMA * delta * delta)
    xk = (1 + 1j) * A_C / delta
    m = -2 * math.pi * A_C**3 * (B_0 / MU_0) * (
        1 - 3 / xk * np.cosh(xk) / np.sinh(xk) + 3 / xk**2)
    return omega, m


def _build_mesh():
    cond = Sphere(Pnt(0, 0, 0), A_C)
    outer = Sphere(Pnt(0, 0, 0), R_K)
    for f in cond.faces:
        f.name = "cond_surf"
    for f in outer.faces:
        f.name = "kelvin_int"
    shell = outer - cond
    shell.mat("inner_air")
    cond.mat("conductor")
    cond.maxh = MAXH_COND

    kext = Sphere(Pnt(*OFFSET), R_K)
    kext.mat("kelvin_air")
    for f in kext.faces:
        f.name = "kelvin_ext"

    gnd = Vertex(Pnt(*OFFSET))
    gnd.name = "GND"

    geo = Glue([cond, shell, kext, gnd])
    int_face = ext_face = None
    for s in geo.solids:
        for f in s.faces:
            if f.name == "kelvin_int":
                int_face = f
            elif f.name == "kelvin_ext":
                ext_face = f
    assert int_face is not None and ext_face is not None
    int_face.Identify(ext_face, "periodic", IdentificationType.PERIODIC)

    with TaskManager():
        mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=MAXH, grading=0.5))
        mesh.Curve(2)
    return mesh


def _material_cfs(mesh):
    ox, oy, oz = OFFSET
    rho2 = (x - ox) ** 2 + (y - oy) ** 2 + (z - oz) ** 2 + 1e-24
    nu_d, sg_d = {}, {}
    for m in mesh.GetMaterials():
        if "kelvin" in m.lower():
            nu_d[m] = NU_0 * rho2 / (R_K * R_K)
            sg_d[m] = 0.0
        elif m == "conductor":
            nu_d[m] = NU_0
            sg_d[m] = SIGMA
        else:
            nu_d[m] = NU_0
            sg_d[m] = 0.0
    return (mesh.MaterialCF(nu_d, default=NU_0),
            mesh.MaterialCF(sg_d, default=0.0))


def _fes_verify(mesh, p):
    """The Verify-First trio: slaved DOF count + kelvin boundary ratio."""
    fb = HCurl(mesh, order=p, complex=True, nograds=True)
    fp = Periodic(fb)
    slaved = sum(fb.FreeDofs()) - sum(fp.FreeDofs())
    h1p = Periodic(H1(mesh, order=max(p, 1)))
    gt = GridFunction(h1p)
    gt.vec[:] = 0.0
    gt.Set(1.0, definedon=mesh.Boundaries("kelvin_int"))
    ratio = (Integrate(gt * gt, mesh, definedon=mesh.Boundaries("kelvin_ext"))
             / Integrate(gt * gt, mesh,
                         definedon=mesh.Boundaries("kelvin_int")))
    return slaved, ratio


def _moment_z(mesh, J):
    return complex(Integrate(0.5 * (x * J[1] - y * J[0]), mesh,
                             definedon=mesh.Materials("conductor"), order=8))


def _solve_astar(mesh, p, omega):
    nu_cf, sg_cf = _material_cfs(mesh)
    s = 1j * omega
    A_s = CF((-y, x, 0)) * (B_0 / 2)
    non_kelvin = "|".join(m for m in mesh.GetMaterials()
                          if "kelvin" not in m.lower())

    fes = Periodic(HCurl(mesh, order=p, complex=True, nograds=True))
    u, v = fes.TnT()
    with TaskManager():
        a = BilinearForm(fes, symmetric=True)
        a += nu_cf * InnerProduct(curl(u), curl(v)) * dx(bonus_intorder=4)
        a += s * sg_cf * InnerProduct(u, v) * dx("conductor")
        a += 1e-6 * NU_0 * InnerProduct(u, v) * dx(non_kelvin)
        f = LinearForm(fes)
        f += -s * sg_cf * InnerProduct(A_s, v) * dx("conductor",
                                                    bonus_intorder=2)
        a.Assemble()
        f.Assemble()
        gf = GridFunction(fes)
        gf.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec

    return _moment_z(mesh, -s * SIGMA * (A_s + gf))


def _solve_aphi(mesh, p, omega):
    nu_cf, sg_cf = _material_cfs(mesh)
    s = 1j * omega
    A_s = CF((-y, x, 0)) * (B_0 / 2)
    non_kelvin = "|".join(m for m in mesh.GetMaterials()
                          if "kelvin" not in m.lower())

    fesA = Periodic(HCurl(mesh, order=p, complex=True, nograds=True))
    fesV = H1(mesh, order=p, complex=True,
              definedon=mesh.Materials("conductor"))
    fes = fesA * fesV
    (u, w), (v, q) = fes.TnT()
    with TaskManager():
        a = BilinearForm(fes, symmetric=True)
        a += nu_cf * InnerProduct(curl(u), curl(v)) * dx(bonus_intorder=4)
        a += s * sg_cf * InnerProduct(u + grad(w),
                                      v + grad(q)) * dx("conductor")
        a += 1e-6 * NU_0 * InnerProduct(u, v) * dx(non_kelvin)
        # pin the additive constant of W = V/s in the conductor
        a += 1e-6 * SIGMA * w * q * dx("conductor")
        f = LinearForm(fes)
        f += -s * sg_cf * InnerProduct(A_s, v + grad(q)) * dx(
            "conductor", bonus_intorder=2)
        a.Assemble()
        f.Assemble()
        gf = GridFunction(fes)
        gf.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec

    gfA, gfW = gf.components
    return _moment_z(mesh, -s * SIGMA * (A_s + gfA + grad(gfW)))


def _results():
    """Solve (p=1,2) x (A*, A-Phi) once and cache."""
    if not _CACHE:
        omega, m_ana = _analytic()
        mesh = _build_mesh()
        _CACHE["m_ana"] = m_ana
        for p in (1, 2):
            slaved, ratio = _fes_verify(mesh, p)
            _CACHE[f"fes_p{p}"] = (slaved, ratio)
            _CACHE[f"astar_p{p}"] = _solve_astar(mesh, p, omega)
            _CACHE[f"aphi_p{p}"] = _solve_aphi(mesh, p, omega)
    return _CACHE


def _rel(m_fem, m_ana):
    return abs(m_fem - m_ana) / abs(m_ana)


@pytest.mark.slow
def test_fes_verify_trio_at_p1_and_p2():
    """Periodic Kelvin identification actually constrains DOFs at p=1 AND p=2."""
    res = _results()
    for p in (1, 2):
        slaved, ratio = res[f"fes_p{p}"]
        print(f"\n  p={p}: slaved={slaved}, kelvin ratio={ratio:.8f}")
        assert slaved > 0, f"Periodic slaved no HCurl DOFs at p={p}"
        assert abs(ratio - 1.0) < 1e-8, \
            f"kelvin_int/ext functional ratio != 1 at p={p}: {ratio}"
    s1, _ = res["fes_p1"]
    s2, _ = res["fes_p2"]
    assert s2 > s1, "p=2 must slave MORE (high-order) DOFs than p=1"


@pytest.mark.slow
def test_aphi_p2_matches_analytic():
    """THE question: A-Phi + Kelvin at p=2 reproduces the analytic moment."""
    res = _results()
    rel = _rel(res["aphi_p2"], res["m_ana"])
    print(f"\n  A-Phi p=2: m = {res['aphi_p2']:.6e}")
    print(f"  analytic : m = {res['m_ana']:.6e}")
    print(f"  rel = {rel*100:.3f}%  (measured 0.053%, band 0.2%)")
    assert rel < 0.002, f"A-Phi p=2 off by {rel*100:.3f}% (band 0.2%)"


@pytest.mark.slow
def test_astar_p2_within_its_floor_band():
    """A-method at p=2 lands near its ~0.45% charge-conservation floor."""
    res = _results()
    rel = _rel(res["astar_p2"], res["m_ana"])
    print(f"\n  A* p=2 rel = {rel*100:.3f}%  (measured 0.473%, band 1.5%)")
    assert rel < 0.015, f"A* p=2 off by {rel*100:.3f}% (band 1.5%)"


@pytest.mark.slow
def test_aphi_beats_astar_at_p2():
    """The Phi block pays: at p=2 A-Phi is several times closer than A*.

    Measured factor ~9 (0.473% / 0.053%); locked conservatively at >= 3.
    If this ever fails with both lanes accurate, the A* floor was fixed
    upstream -- celebrate, re-measure, and retighten the A* band instead.
    """
    res = _results()
    rel_astar = _rel(res["astar_p2"], res["m_ana"])
    rel_aphi = _rel(res["aphi_p2"], res["m_ana"])
    print(f"\n  p=2: A* {rel_astar*100:.3f}% vs A-Phi {rel_aphi*100:.3f}% "
          f"(factor {rel_astar/max(rel_aphi, 1e-30):.1f})")
    assert rel_aphi * 3.0 < rel_astar, (
        f"expected A-Phi to beat A* by >= 3x at p=2, got "
        f"{rel_astar*100:.3f}% vs {rel_aphi*100:.3f}%")


@pytest.mark.slow
def test_p1_both_lanes_agree_and_land_in_band():
    """At p=1 the lanes coincide (W inert) and sit at the p=1 mesh error."""
    res = _results()
    m_ana = res["m_ana"]
    rel_astar = _rel(res["astar_p1"], m_ana)
    rel_aphi = _rel(res["aphi_p1"], m_ana)
    lane_gap = abs(res["astar_p1"] - res["aphi_p1"]) / abs(m_ana)
    print(f"\n  p=1: A* {rel_astar*100:.3f}%, A-Phi {rel_aphi*100:.3f}%, "
          f"lane gap {lane_gap:.2e}")
    assert rel_astar < 0.05 and rel_aphi < 0.05, \
        f"p=1 outside 5% band: {rel_astar*100:.2f}% / {rel_aphi*100:.2f}%"
    assert lane_gap < 1e-4, \
        f"p=1 lanes expected (near-)identical, gap = {lane_gap:.2e}"


if __name__ == "__main__":
    res = _results()
    m_ana = res["m_ana"]
    print(f"m_analytic = {m_ana:.6e}")
    for p in (1, 2):
        slaved, ratio = res[f"fes_p{p}"]
        print(f"p={p}: slaved={slaved} ratio={ratio:.6f} | "
              f"A* {_rel(res[f'astar_p{p}'], m_ana)*100:.3f}% | "
              f"A-Phi {_rel(res[f'aphi_p{p}'], m_ana)*100:.3f}%")
