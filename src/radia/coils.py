"""Filament current-loop GEOMETRY for Radia (circular loop, Helmholtz pair).

WHY THIS EXISTS
---------------
The companion :mod:`radia.round_bodies` builds SOLID magnetized bodies (cylinder,
sphere) from permanent-magnet material.  This module is the CURRENT-SOURCE
complement: thin filament coils built from ``ObjFlmCur`` (a closed polyline
current), for field synthesis, Helmholtz / bias coils and inductance studies.
A circular loop is approximated by a regular ``nseg``-gon filament; the on-axis
faceting error is O(1/nseg^2).

UNITS
-----
Current is in Amperes (Radia ``ObjFlmCur`` convention).  Returned helper fields
(:func:`loop_axial_field`, :func:`helmholtz_center_field`) are in Tesla.

VERIFIED (run via PowerShell -- see tests/test_coils.py, and the ELF/MAGIC
cross-check in C:\\temp\\radia_helmholtz.py):
  * circular_loop on-axis field == single-loop closed form mu0 I R^2 /
    (2 (R^2+z^2)^{3/2}) to <0.05% (nseg=96, convergent in nseg).
  * helmholtz_pair (separation = radius) centre field == (4/5)^{3/2} mu0 I / R
    to <0.05% (nseg=96).
  * the exact 24-gon Helmholtz filament reproduces an independent Biot-Savart
    reference (ELF/MAGIC MOMENT solver) to 2e-12 on the whole axis.
"""

import math
import radia as rad

MU0 = 4.0e-7 * math.pi

# perpendicular-plane axes for each coil axis (the loop lives in this plane)
_PERP = {"z": (0, 1), "x": (1, 2), "y": (2, 0)}
_IA = {"x": 0, "y": 1, "z": 2}


def _loop_points(center, radius, axis, nseg):
    """Closed regular nseg-gon polyline (filament) of circum-radius ``radius``
    centred at ``center`` in the plane perpendicular to ``axis``."""
    iu, iv = _PERP[axis]
    ia = _IA[axis]
    pts = []
    for k in range(nseg):
        th = 2.0 * math.pi * k / nseg
        p = [0.0, 0.0, 0.0]
        p[ia] = center[ia]
        p[iu] = center[iu] + radius * math.cos(th)
        p[iv] = center[iv] + radius * math.sin(th)
        pts.append(p)
    pts.append(list(pts[0]))                 # close the loop
    return pts


def circular_loop(center, radius, current, axis="z", nseg=48):
    """Single circular current loop as a closed ``nseg``-gon filament.

    Args:
        center: [x, y, z] loop centre.
        radius: loop radius (circum-radius of the N-gon).
        current: loop current [A] (sign sets circulation by the right-hand rule
                 about +axis).
        axis: coil axis "x" | "y" | "z".
        nseg: number of polygon sides (>= 3).

    Returns:
        Radia object handle (ObjFlmCur).
    """
    axis = axis.lower()
    if axis not in _PERP:
        raise ValueError("axis must be 'x', 'y' or 'z'")
    if nseg < 3:
        raise ValueError("nseg must be >= 3")
    return rad.ObjFlmCur(_loop_points(center, radius, axis, nseg), current)


def helmholtz_pair(center, radius, current, axis="z", nseg=48):
    """Helmholtz coil pair: two coaxial loops separated by ``radius`` (the
    Helmholtz condition), same circulation, centred on ``center``.

    Args:
        center: [x, y, z] midpoint between the two loops.
        radius: loop radius == loop separation (Helmholtz condition).
        current: current [A] in each loop (same sense -> aiding fields).
        axis: coil axis "x" | "y" | "z".
        nseg: polygon sides per loop.

    Returns:
        Radia container (ObjCnt) of the two loops.
    """
    axis = axis.lower()
    if axis not in _PERP:
        raise ValueError("axis must be 'x', 'y' or 'z'")
    ia = _IA[axis]
    c0 = list(center); c0[ia] = center[ia] - 0.5 * radius
    c1 = list(center); c1[ia] = center[ia] + 0.5 * radius
    return rad.ObjCnt([
        rad.ObjFlmCur(_loop_points(c0, radius, axis, nseg), current),
        rad.ObjFlmCur(_loop_points(c1, radius, axis, nseg), current),
    ])


def solenoid(center, radius, length, current, nturns, axis="z", nseg=48):
    """Finite solenoid as a stack of ``nturns`` coaxial circular-loop filaments of
    ``radius`` evenly spaced over ``length``, each carrying ``current``.

    Args:
        center: [x, y, z] centre of the solenoid.
        radius: winding radius.
        length: axial length.
        current: per-turn current [A].
        nturns: number of turns (loops in the stack).
        axis: coil axis "x" | "y" | "z".
        nseg: polygon sides per loop.

    Returns:
        Radia container (ObjCnt) of the ``nturns`` loops.
    """
    axis = axis.lower()
    if axis not in _PERP:
        raise ValueError("axis must be 'x', 'y' or 'z'")
    nturns = int(nturns)
    if nturns < 1:
        raise ValueError("nturns must be >= 1")
    ia = _IA[axis]
    loops = []
    for k in range(nturns):
        # turn centres span the full length symmetrically about `center`
        frac = 0.0 if nturns == 1 else (k / (nturns - 1) - 0.5)
        c = list(center)
        c[ia] = center[ia] + frac * length
        loops.append(rad.ObjFlmCur(_loop_points(c, radius, axis, nseg), current))
    return rad.ObjCnt(loops)


def loop_axial_field(radius, current, z):
    """On-axis flux density [T] of a single circular loop, distance ``z`` from
    its plane:  B = mu0 I R^2 / (2 (R^2 + z^2)^{3/2})  (closed form)."""
    return MU0 * current * radius * radius / (2.0 * (radius * radius + z * z) ** 1.5)


def solenoid_axial_field(radius, length, nturns, current, z):
    """On-axis flux density [T] of a thin finite solenoid (radius R, length L, N turns,
    per-turn current I), at axial position ``z`` measured from the solenoid centre:

        B(z) = (mu0 N I / 2L) [ (z+L/2)/sqrt(R^2+(z+L/2)^2) - (z-L/2)/sqrt(R^2+(z-L/2)^2) ].

    Centre (z=0): B = mu0 N I (L/2) / (L sqrt(R^2+(L/2)^2)).  Infinite-solenoid limit
    L>>R: B -> mu0 N I / L = mu0 n I."""
    n = nturns / length
    zp, zm = z + 0.5 * length, z - 0.5 * length
    a = zp / math.sqrt(radius * radius + zp * zp)
    b = zm / math.sqrt(radius * radius + zm * zm)
    return 0.5 * MU0 * n * current * (a - b)


def helmholtz_center_field(radius, current):
    """Centre flux density [T] of an ideal Helmholtz pair (separation = radius):
    B = (4/5)^{3/2} mu0 I / R."""
    return (0.8 ** 1.5) * MU0 * current / radius
