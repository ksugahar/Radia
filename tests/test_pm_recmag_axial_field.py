"""
radia rectangular permanent-magnet field vs the closed-form on-axis field of a uniformly
magnetized rectangular block (cuboid).

The canonical rectangular PM is now ``rad.magnet_box`` -- an MMMM **surface-charge**
``ObjHexahedron`` (the Python-facing surface-current ``ObjRecMag`` constructor was retired;
``rad.ObjRecMag`` is a thin shim forwarding to ``magnet_box``).  For a uniformly magnetized
block the surface-charge field IS the exact analytic (Kennelly & Joch / surface-charge K&J)
cuboid field, so it reproduces the closed form to machine precision -- in fact tighter than
the old surface-current ``ObjRecMag`` (~1e-16 on-axis), so this test asserts to 1e-12.
(An independent FEM reference on the same magnet agreed to ~1-6%, linear-tet discretization
-- stored internally as a regression reference; not part of this open analytic-gated test.)

KEY radia gotcha (locked by test_magnetization_units_are_A_per_m): magnetization is in A/m
(SI M), NOT Tesla.  For a remanence Br [T] pass M = Br/mu0 [A/m]; rad.Fld 'b' returns B in
Tesla; lengths are in mm here.  A Tesla-valued magnetization gives a field mu0x too small
(the classic mistake).

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

# The MMMM surface-charge magnet_box reproduces the K&J closed form on-axis to ~1e-16
# near the face, growing to ~3e-14 in the far field (atan cancellation).  Assert 1e-12:
# far tighter than the old surface-current ObjRecMag's 1e-6 gate, strict but cancellation-safe.
TOL = 1e-12


def _analytic(zsurf):
    def term(d):
        return math.atan((WX * WY) / (2 * d * math.sqrt(4 * d * d + WX * WX + WY * WY)))
    return (BR / math.pi) * (term(zsurf) - term(zsurf + WZ))


def test_recmag_onaxis_matches_closed_form():
    """radia magnet_box on-axis B_z == closed form to machine precision (MMMM surface-charge
    field IS the exact analytic cuboid field)."""
    mag = rad.magnet_box([0, 0, 0], [WX, WY, WZ], [0, 0, BR / MU0])   # M in A/m
    worst = 0.0
    for zsurf in (2.0, 5.0, 10.0, 20.0, 40.0, 80.0):
        bz = rad.Fld(mag, "b", [0, 0, WZ / 2.0 + zsurf])[2]
        ba = _analytic(zsurf)
        rel = abs(bz - ba) / abs(ba)
        assert rel < TOL, f"zsurf={zsurf}mm: radia {bz:.6e} vs analytic {ba:.6e} (rel {rel:.2e})"
        worst = max(worst, rel)
    assert worst < TOL


def test_magnetization_units_are_A_per_m():
    """Guard the A/m (not Tesla) magnetization convention: a Tesla-valued M is mu0x too small."""
    z = WZ / 2.0 + 10.0
    bz_Am = rad.Fld(rad.magnet_box([0, 0, 0], [WX, WY, WZ], [0, 0, BR / MU0]), "b", [0, 0, z])[2]
    bz_T = rad.Fld(rad.magnet_box([0, 0, 0], [WX, WY, WZ], [0, 0, BR]), "b", [0, 0, z])[2]
    assert abs(bz_T / bz_Am - MU0) / MU0 < 1e-6          # Tesla-input == mu0 x A/m-input
    assert abs(bz_Am - _analytic(10.0)) / _analytic(10.0) < TOL
