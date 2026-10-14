"""Verify radia.round_bodies solid round-magnet builders (cylinder, sphere).

These exist because ObjArcPgnMag heap-corrupts on axis-touching (solid) revolves;
round_bodies builds solids from ObjThckPgn (extruded N-gon) instead. Run via
PowerShell (radia imports fine there). Magnetization in A/m (M = Br/mu0)."""
import os
import sys
import math
import pytest

rad = pytest.importorskip("radia")

_SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "src", "radia")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
import round_bodies  # noqa: E402

MU0 = 4.0e-7 * math.pi
BR = 1.2
M = BR / MU0


def _cyl_axial_analytic(R, L, z):
    a = (z + L / 2) / math.sqrt((z + L / 2) ** 2 + R ** 2)
    b = (z - L / 2) / math.sqrt((z - L / 2) ** 2 + R ** 2)
    return 0.5 * BR * (a - b)


def test_cylinder_axial_z():
    rad.UtiDelAll()
    R, L = 10.0, 20.0
    c = round_bodies.cyl_mag([0, 0, 0], R, L, [0, 0, M], axis="z", nseg=72)
    for z in (0.0, 5.0, 10.0, 20.0, 30.0):
        b = rad.Fld(c, "bz", [0, 0, z])
        a = _cyl_axial_analytic(R, L, z)
        assert abs(b - a) / abs(a) < 3e-3, (z, b, a)
    rad.UtiDelAll()


def test_cylinder_axis_x():
    """axis='x' cylinder: on-axis Bx == the same axial closed form (symmetry)."""
    rad.UtiDelAll()
    R, L = 10.0, 20.0
    c = round_bodies.cyl_mag([0, 0, 0], R, L, [M, 0, 0], axis="x", nseg=72)
    for x in (0.0, 10.0, 30.0):
        b = rad.Fld(c, "bx", [x, 0, 0])
        a = _cyl_axial_analytic(R, L, x)
        assert abs(b - a) / abs(a) < 3e-3, (x, b, a)
    rad.UtiDelAll()


def test_sphere_center_demag():
    """Uniformly magnetized sphere center field == (2/3)Br (demag N=1/3)."""
    rad.UtiDelAll()
    s = round_bodies.sphere_mag([0, 0, 0], 10.0, [0, 0, M], nz=80, nseg=48)
    b = rad.Fld(s, "bz", [0, 0, 0])
    a = 2.0 / 3.0 * BR
    assert abs(b - a) / a < 2e-3, (b, a)
    rad.UtiDelAll()
