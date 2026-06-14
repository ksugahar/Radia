"""
radia permanent-magnet field core vs the closed-form on-axis field of a uniformly
magnetized rectangular block (cuboid).  radia's ObjRecMag field IS the exact analytic
field of a uniformly magnetized cuboid, so it must reproduce the closed form to machine
precision.  (An independent FEM reference on the same magnet agreed to ~1-6%, linear-tet
discretization -- stored internally as a regression reference; not part of this open
analytic-gated test.)

KEY radia gotcha (locked by test_magnetization_units_are_A_per_m): ObjRecMag magnetization
is in A/m (SI M), NOT Tesla.  For a remanence Br [T] pass M = Br/mu0 [A/m]; rad.Fld 'b'
returns B in Tesla; lengths are in mm.  A Tesla-valued magnetization gives a field mu0x
too small (the classic mistake).

On-axis B_z of a block (cross-section Wx x Wy, thickness Wz, centered, remanence Br) at
distance zsurf above the pole face (surface-charge / K&J closed form):
    B(zsurf) = (Br/pi)[ atan( Wx Wy / (2 d sqrt(4 d^2 + Wx^2 + Wy^2)) ) ]
               evaluated at d=zsurf minus the same at d=zsurf+Wz.
"""
import math
import pytest

rad = pytest.importorskip("radia")

MU0 = 4e-7 * math.pi
BR = 1.2            # remanence [T]
WX = WY = 20.0      # block cross-section [mm]
WZ = 10.0           # block thickness [mm]


def _analytic(zsurf):
    def term(d):
        return math.atan((WX * WY) / (2 * d * math.sqrt(4 * d * d + WX * WX + WY * WY)))
    return (BR / math.pi) * (term(zsurf) - term(zsurf + WZ))


def test_recmag_onaxis_matches_closed_form():
    """radia ObjRecMag on-axis B_z == closed form to machine precision (it IS analytic)."""
    mag = rad.ObjRecMag([0, 0, 0], [WX, WY, WZ], [0, 0, BR / MU0])   # M in A/m
    worst = 0.0
    for zsurf in (2.0, 5.0, 10.0, 20.0, 40.0, 80.0):
        bz = rad.Fld(mag, "b", [0, 0, WZ / 2.0 + zsurf])[2]
        ba = _analytic(zsurf)
        rel = abs(bz - ba) / abs(ba)
        assert rel < 1e-6, f"zsurf={zsurf}mm: radia {bz:.6e} vs analytic {ba:.6e} (rel {rel:.2e})"
        worst = max(worst, rel)
    assert worst < 1e-6


def test_magnetization_units_are_A_per_m():
    """Guard the A/m (not Tesla) magnetization convention: a Tesla-valued M is mu0x too small."""
    z = WZ / 2.0 + 10.0
    bz_Am = rad.Fld(rad.ObjRecMag([0, 0, 0], [WX, WY, WZ], [0, 0, BR / MU0]), "b", [0, 0, z])[2]
    bz_T = rad.Fld(rad.ObjRecMag([0, 0, 0], [WX, WY, WZ], [0, 0, BR]), "b", [0, 0, z])[2]
    assert abs(bz_T / bz_Am - MU0) / MU0 < 1e-6          # Tesla-input == mu0 x A/m-input
    assert abs(bz_Am - _analytic(10.0)) / _analytic(10.0) < 1e-6
