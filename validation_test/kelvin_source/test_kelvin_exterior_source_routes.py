"""Golden: the external field MUST be carried into the Kelvin exterior, and
the 1-form Convention B rule is the one that does it correctly.

A magnetic sphere (radius a, mu_r) sits in a uniform applied field H0 z_hat.
The problem is solved with the Kelvin open boundary and a REDUCED scalar
potential, H = H_s - grad(Omega), so the analytical interior field is

    H_in = 3 H0 / (mu_r + 2).

All three routes below share one mesh, one bilinear form, one solver.  They
differ ONLY in what the source field H_s is inside the Kelvin exterior:

    Z   H_s' = 0                       (source dropped in kext)
    B1  H_s' = -(rho'/R)^2 H0 z_hat    Convention B, 1-form helper
                                       make_reduced_potential_background_cf
    B0  H_s' = -grad(Omega_s'),        Convention B, 0-form helper
        Omega_s' = -(rho'/R)^2 Omega_s(local)
                                       make_reduced_potential_scalar_cf

MEASURED (LAB, 2026-07-23), showing clean mesh convergence:

    config                 Z          B1         B0
    maxh .14 / order 2   -32.830%   +0.743%   +34.329%
    maxh .11 / order 3   -33.331%   +0.002%   +33.337%

    ratio Z/B1  -> 0.666669  (exactly 2/3)
    ratio B0/B1 -> 1.333338  (exactly 4/3)

Conclusions locked here:

1. B1 converges to the analytical answer (+0.002% at p=3), so the 1-form
   Convention B rule is CORRECT for a field-driven reduced-potential weak
   form.
2. Dropping the source in the Kelvin exterior is NOT a small approximation
   for a globally-defined background: it loses exactly one third of the
   field (a factor 2/3).  This is what a solver does when it excludes the
   Kelvin region from the source projection; for a compact coil the error
   shrinks with H_s at the Kelvin radius, but for a background applied at
   infinity it is maximal.
3. Differentiating the 0-form Convention B potential into a field OVERSHOOTS
   by exactly one third (a factor 4/3).  This is the numerical face of the
   symbolic result locked in tests/test_reduced_potential_background.py:
   curl(H_s' 1-form B) != 0, so the two Convention B rules are separate
   contracts and neither is the derivative of the other.

The ratios 2/3 and 4/3 are essentially mesh-independent and are the
strongest assertions in this file.

Reference: docs/kelvin/KELVIN_TRANSFORMATION.md 7.4 / 7.6.
"""

from __future__ import annotations

import math

import pytest

from netgen.meshing import IdentificationType
from netgen.occ import Glue, OCCGeometry, Pnt, Sphere, Vertex
from ngsolve import (
    BilinearForm, CoefficientFunction as CF, GridFunction, H1, LinearForm,
    Mesh, Periodic, TaskManager, dx, grad, x, y, z,
)

from radia.kelvin_material import (
    make_reduced_potential_background_cf,
    make_reduced_potential_scalar_cf,
)

MU_0 = 4e-7 * math.pi
A_SPH = 0.5          # magnetic sphere radius [m]
R_K = 1.0            # Kelvin radius [m]
OFFSET = (3.0, 0.0, 0.0)
MU_R = 100.0
H_0 = 1.0

MAXH = 0.14          # CI-friendly; the p=3 refinement is recorded above
ORDER = 2

H_ANALYTIC = 3.0 * H_0 / (MU_R + 2.0)

_CACHE: dict[str, float] = {}


def _build_mesh():
    mag = Sphere(Pnt(0, 0, 0), A_SPH)
    outer = Sphere(Pnt(0, 0, 0), R_K)
    for f in mag.faces:
        f.name = "sphere"
    for f in outer.faces:
        f.name = "kelvin_int"
    shell = outer - mag
    shell.mat("inner_air")
    mag.mat("magnetic")

    kext = Sphere(Pnt(*OFFSET), R_K)
    kext.mat("kelvin_air")
    for f in kext.faces:
        f.name = "kelvin_ext"

    gnd = Vertex(Pnt(*OFFSET))
    gnd.name = "GND"

    geo = Glue([mag, shell, kext, gnd])

    int_face = ext_face = None
    for s in geo.solids:
        for f in s.faces:
            if f.name == "kelvin_int":
                int_face = f
            elif f.name == "kelvin_ext":
                ext_face = f
    assert int_face is not None and ext_face is not None, \
        "kelvin_int / kelvin_ext faces not found on the glued geometry"
    int_face.Identify(ext_face, "periodic", IdentificationType.PERIODIC)

    with TaskManager():
        mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=MAXH, grading=0.5))
        mesh.Curve(ORDER)
    return mesh


def _source_cf(mesh, route):
    """Source field H_s for the given Kelvin-exterior route."""
    if route == "Z":
        is_inner = mesh.MaterialCF(
            {m: (0.0 if "kelvin" in m.lower() else 1.0)
             for m in mesh.GetMaterials()}, default=1.0)
        return CF((0, 0, H_0 * is_inner))
    if route == "B1":
        return make_reduced_potential_background_cf(
            mesh, lambda xc, yc, zc: CF((0, 0, H_0)),
            R_K=R_K, offset=OFFSET, kelvin_mats=("kelvin",), dim=3)
    if route == "B0":
        omega_s = make_reduced_potential_scalar_cf(
            mesh, lambda xc, yc, zc: -H_0 * zc,
            R_K=R_K, offset=OFFSET, kelvin_mats=("kelvin",), dim=3)
        return -CF((omega_s.Diff(x), omega_s.Diff(y), omega_s.Diff(z)))
    raise ValueError(f"unknown route {route!r}")


def _solve(mesh, route):
    """Reduced-potential solve; returns H_z at the sphere centre."""
    ox, oy, oz = OFFSET
    rho2 = (x - ox) ** 2 + (y - oy) ** 2 + (z - oz) ** 2 + 1e-24

    mu_dict = {}
    for m in mesh.GetMaterials():
        ml = m.lower()
        if "kelvin" in ml:
            mu_dict[m] = MU_0 * (R_K * R_K) / rho2      # mu' = (R/rho')^2 mu0
        elif "magnetic" in ml:
            mu_dict[m] = MU_0 * MU_R
        else:
            mu_dict[m] = MU_0
    mu = mesh.MaterialCF(mu_dict, default=MU_0)

    h_s = _source_cf(mesh, route)
    fes = Periodic(H1(mesh, order=ORDER, dirichlet_bbnd="GND"))
    u, v = fes.TnT()

    with TaskManager():
        a = BilinearForm(fes, symmetric=True)
        a += mu * grad(u) * grad(v) * dx(bonus_intorder=6)
        f = LinearForm(fes)
        f += mu * h_s * grad(v) * dx(bonus_intorder=6)
        a.Assemble()
        f.Assemble()
        gf = GridFunction(fes)
        gf.vec.data = a.mat.Inverse(
            fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    return H_0 - grad(gf)[2](mesh(0.0, 0.0, 0.0))


def _results():
    """Solve all three routes once and cache."""
    if not _CACHE:
        mesh = _build_mesh()
        fes_plain = H1(mesh, order=ORDER, dirichlet_bbnd="GND")
        slaved = sum(fes_plain.FreeDofs()) - sum(Periodic(fes_plain).FreeDofs())
        assert slaved > 0, \
            f"Kelvin periodic identification constrains no DOF (slaved={slaved})"
        for route in ("Z", "B1", "B0"):
            _CACHE[route] = _solve(mesh, route)
    return _CACHE


@pytest.mark.slow
def test_convention_b_1form_reproduces_analytic_sphere():
    """B1: carrying the background into kext by the 1-form rule is correct."""
    h_z = _results()["B1"]
    err = 100.0 * (h_z / H_ANALYTIC - 1.0)
    print(f"\n  B1 H_z(0) = {h_z:.8f}, analytic = {H_ANALYTIC:.8f}, "
          f"err = {err:+.3f}%")
    assert abs(err) < 1.5, f"1-form Convention B off by {err:+.3f}% (band 1.5%)"


@pytest.mark.slow
def test_dropping_the_exterior_source_loses_one_third():
    """Z: omitting H_s in kext is a 33% error, not a small approximation."""
    h_z = _results()["Z"]
    err = 100.0 * (h_z / H_ANALYTIC - 1.0)
    print(f"\n  Z  H_z(0) = {h_z:.8f}, err = {err:+.3f}%")
    assert err < -20.0, (
        "expected a large deficit when the Kelvin exterior carries no source; "
        f"got {err:+.3f}%")


@pytest.mark.slow
def test_zero_exterior_source_ratio_is_exactly_two_thirds():
    """The deficit is the mesh-independent factor 2/3, not mesh error."""
    res = _results()
    ratio = res["Z"] / res["B1"]
    print(f"\n  Z / B1 = {ratio:.6f}  (2/3 = {2/3:.6f})")
    assert abs(ratio - 2.0 / 3.0) < 1e-3, \
        f"Z/B1 ratio drifted from 2/3: {ratio:.6f}"


@pytest.mark.slow
def test_0form_convention_b_gradient_overshoots_by_one_third():
    """B0: -grad(0-form Convention B) is NOT the 1-form Convention B field."""
    h_z = _results()["B0"]
    err = 100.0 * (h_z / H_ANALYTIC - 1.0)
    print(f"\n  B0 H_z(0) = {h_z:.8f}, err = {err:+.3f}%")
    assert err > 20.0, (
        "expected the differentiated 0-form rule to overshoot; "
        f"got {err:+.3f}%")


@pytest.mark.slow
def test_0form_gradient_ratio_is_exactly_four_thirds():
    """The overshoot is the mesh-independent factor 4/3."""
    res = _results()
    ratio = res["B0"] / res["B1"]
    print(f"\n  B0 / B1 = {ratio:.6f}  (4/3 = {4/3:.6f})")
    assert abs(ratio - 4.0 / 3.0) < 1e-3, \
        f"B0/B1 ratio drifted from 4/3: {ratio:.6f}"


if __name__ == "__main__":
    res = _results()
    print(f"analytic H_in = {H_ANALYTIC:.8f}")
    for route in ("Z", "B1", "B0"):
        err = 100.0 * (res[route] / H_ANALYTIC - 1.0)
        print(f"  {route:3s} H_z(0) = {res[route]:.8f}  err = {err:+.3f}%")
    print(f"  Z / B1  = {res['Z'] / res['B1']:.6f}   (2/3)")
    print(f"  B0 / B1 = {res['B0'] / res['B1']:.6f}   (4/3)")
