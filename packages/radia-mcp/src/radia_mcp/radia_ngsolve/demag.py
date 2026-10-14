# -*- coding: utf-8 -*-
"""Per-element permanent-magnet demagnetization (knee) post-processor.

A demag-distribution capability for radia-ngsolve. The classic 2D approach
checks the demag operating point MANUALLY after a solve ("is every point above
the knee?"), or tracks it with a variable-magnet user-subroutine. radia-ngsolve
does it natively as a post-processor on ANY planar PM solve: it projects the
per-point operating field onto the local magnetization direction and reports the
AREA FRACTION whose H falls below the knee (irreversible-demag risk) -- a full
demag *distribution*, not a single bulk permeance coefficient.

Model (matches ``solve.planar_magnet_source`` convention). The magnet
source is ``s = nu*B_rem = Hc*(cos phi, sin phi)`` and the constitutive law in
the magnet is ``H = nu*B - s``. Projected on the easy axis
``m_hat = (cos phi, sin phi)``:

    B_m = B . m_hat                       flux density along the easy axis
    H_m = B_m/(mu0*mu_r) - Hc             field along the easy axis

with ``mu_r`` the recoil permeability. The recoil line is
``B_m = mu0*mu_r*(H_m + Hc)``, so ``B_m = Br`` at ``H_m = 0`` (Br = mu0*mu_r*Hc)
and ``B_m = 0`` at ``H_m = -Hc``. An element is IRREVERSIBLY demagnetized once
``H_m < H_knee`` (the knee field of the demag curve at the operating
temperature; NdFeB's knee moves toward 0 as temperature rises).

Validated on the transverse-magnetized cylinder (uniform interior operating
point): ``H_m = -mu_r*Hc/(mu_r+1)`` exactly -- see validation_test/radia_mcp/test_demag_knee.py.
"""
from __future__ import annotations

import math

from ngsolve import CoefficientFunction, grad, Integrate, IfPos

MU0 = 4e-7 * math.pi


def _direction_cfs(mesh, magnets):
    """(m_hat_x, m_hat_y, Hc) MaterialCFs over the magnet regions (0 elsewhere)."""
    mx, my, hc = {}, {}, {}
    for mat, (Hc, phi_deg) in magnets.items():
        ph = math.radians(phi_deg)
        mx[mat] = math.cos(ph)
        my[mat] = math.sin(ph)
        hc[mat] = float(Hc)
    return (mesh.MaterialCF(mx, default=0.0),
            mesh.MaterialCF(my, default=0.0),
            mesh.MaterialCF(hc, default=0.0))


def demag_operating_field(gfu, mesh, magnets, mu_r):
    """Return ``(B_m, H_m)`` CoefficientFunctions: flux density and H field
    PROJECTED onto the local magnetization direction, valid inside each magnet
    region.

    Parameters
    ----------
    gfu     : A_z GridFunction from ``solve_planar_magnetostatic[_nonlinear]``.
    mesh    : the 2D mesh.
    magnets : ``{material: (Hc, phi_deg)}`` -- the SAME dict handed to the solve.
    mu_r    : recoil permeability -- a scalar, or ``{material: mu_r}``.
    """
    mx, my, hc = _direction_cfs(mesh, magnets)
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    B_m = B[0] * mx + B[1] * my
    mur_cf = mesh.MaterialCF(mu_r, default=1.0) if isinstance(mu_r, dict) else float(mu_r)
    H_m = B_m / (MU0 * mur_cf) - hc
    return B_m, H_m


def demag_knee_report(gfu, mesh, magnets, mu_r, H_knee, region=None, area_tol=1e-3):
    """Per-element demag-knee safety report over the magnet region(s).

    Parameters
    ----------
    H_knee  : knee field [A/m], NEGATIVE (e.g. ~ -7.5e5 for NdFeB at room
              temperature; rises toward 0 as temperature increases). A point is
              unsafe when its operating ``H_m < H_knee``.
    region  : material name (or ``"a|b"``) to restrict to; default = all magnets.
    area_tol: ``safe`` is True when the below-knee area fraction is below this.

    Returns
    -------
    dict with: ``region``, ``area``, ``mean_Bm``, ``mean_Hm``, ``H_knee``,
    ``area_fraction_below_knee`` (the demag-distribution metric) and
    ``safe`` (bool).
    """
    B_m, H_m = demag_operating_field(gfu, mesh, magnets, mu_r)
    pat = region or "|".join(magnets.keys())
    mats = mesh.Materials(pat)
    area = Integrate(CoefficientFunction(1.0), mesh, definedon=mats)
    mean_Bm = Integrate(B_m, mesh, definedon=mats) / area
    mean_Hm = Integrate(H_m, mesh, definedon=mats) / area
    frac = Integrate(IfPos(H_knee - H_m, 1.0, 0.0), mesh, definedon=mats) / area
    return {
        "region": pat,
        "area": area,
        "mean_Bm": mean_Bm,
        "mean_Hm": mean_Hm,
        "H_knee": H_knee,
        "area_fraction_below_knee": frac,
        "safe": bool(frac < area_tol),
    }
