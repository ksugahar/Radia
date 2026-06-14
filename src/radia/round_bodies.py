"""Solid round permanent-magnet GEOMETRY for Radia (cylinder, sphere).

WHY THIS EXISTS
---------------
Radia's ``ObjArcPgnMag`` (surface-of-revolution) HEAP-CORRUPTS (Windows
0xC0000374) when the revolved cross-section TOUCHES THE ROTATION AXIS (r=0) --
i.e. for any SOLID body of revolution.  An annular cross-section (r>0) builds,
but it bores a hole on the axis, so the ON-AXIS field then reads the air-bore
demag field (negative) instead of the solid interior field.  Solid cylinders /
spheres are therefore built here from ``ObjThckPgn`` (extruded regular N-gon) --
no revolve, no axis singularity, correct on-axis field.

UNITS
-----
Magnetization is in A/m (Radia 4.32 ObjRecMag/ObjThckPgn convention).  For a
remanence Br [T] pass ``M = Br / mu0`` (mu0 = 4*pi*1e-7), i.e. Br * 7.9577e5.

VERIFIED (R=10, L=20, Br=1.2; run via PowerShell -- see tests/test_round_bodies.py):
  * cyl_mag axial on-axis field == equivalent-solenoid analytic to ~1e-3
    (N-gon faceting, convergent in nseg).
  * sphere_mag center field == (2/3)Br (demag N=1/3) to ~6e-4 (convergent in nz).
"""

import math
import radia as rad

# perpendicular-plane axes for each extrusion axis (ObjThckPgn polygon lives here)
_PERP = {"z": (0, 1), "x": (1, 2), "y": (2, 0)}


def _ngon(cu, cv, R, n):
    """Regular n-gon of circum-radius R centred at (cu, cv) in a 2D plane."""
    return [[cu + R * math.cos(2.0 * math.pi * i / n),
             cv + R * math.sin(2.0 * math.pi * i / n)] for i in range(n)]


def cyl_mag(center, radius, height, magnetization, axis="z", nseg=72):
    """Solid cylinder magnet as an extruded regular ``nseg``-gon.

    Avoids the ObjArcPgnMag axis crash; the polygon faceting error in the field
    is O(1/nseg^2) (nseg=72 -> ~1e-3).

    Args:
        center: [x, y, z] centre of the cylinder.
        radius: cylinder radius (circum-radius of the N-gon).
        height: full length along ``axis``.
        magnetization: [Mx, My, Mz] in A/m (global frame).
        axis: extrusion axis "x" | "y" | "z".
        nseg: number of polygon sides (>= 3).

    Returns:
        Radia object handle.
    """
    axis = axis.lower()
    if axis not in _PERP:
        raise ValueError("axis must be 'x', 'y' or 'z'")
    if nseg < 3:
        raise ValueError("nseg must be >= 3")
    iu, iv = _PERP[axis]
    ia = {"x": 0, "y": 1, "z": 2}[axis]
    poly = _ngon(center[iu], center[iv], radius, nseg)
    return rad.ObjThckPgn(center[ia], height, poly, axis, list(magnetization))


def sphere_mag(center, radius, magnetization, nz=80, nseg=48, axis="z"):
    """Solid sphere magnet as a stack of extruded polygonal disks.

    Each z-slice is an ObjThckPgn N-gon of radius sqrt(R^2 - z_mid^2); the stack
    is wrapped in an ObjCnt.  Staircase error O(1/nz^2) + faceting O(1/nseg^2)
    (nz=80, nseg=48 -> ~6e-4 on the centre field).

    Args:
        center: [x, y, z] centre of the sphere.
        radius: sphere radius.
        magnetization: [Mx, My, Mz] in A/m (global frame).
        nz: number of stacked disks along ``axis``.
        nseg: polygon sides per disk.
        axis: stacking axis "x" | "y" | "z".

    Returns:
        Radia container (ObjCnt) handle.
    """
    axis = axis.lower()
    if axis not in _PERP:
        raise ValueError("axis must be 'x', 'y' or 'z'")
    ia = {"x": 0, "y": 1, "z": 2}[axis]
    iu, iv = _PERP[axis]
    R = float(radius)
    disks = []
    for k in range(nz):
        a0 = -R + 2.0 * R * k / nz
        a1 = -R + 2.0 * R * (k + 1) / nz
        amid = 0.5 * (a0 + a1)
        dz = a1 - a0
        rk = math.sqrt(max(R * R - amid * amid, 0.0))
        if rk < 1.0e-9:
            continue
        poly = _ngon(center[iu], center[iv], rk, nseg)
        disks.append(rad.ObjThckPgn(center[ia] + amid, dz, poly, axis,
                                    list(magnetization)))
    return rad.ObjCnt(disks)
