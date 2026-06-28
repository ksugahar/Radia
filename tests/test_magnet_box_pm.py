"""Golden: radia.magnet_box -- the MMMM permanent-magnet substitute for the retiring ObjRecMag.

A uniformly magnetized rectangular block has the IDENTICAL external field whether modeled as a
surface-current ObjRecMag or as a surface-charge MMMM ObjHexahedron (sigma = M.n).  magnet_box builds
the ObjHexahedron from (center, dimensions, magnetization) -- same call shape as ObjRecMag -- so
permanent magnets survive the ObjRecMag retirement (CLAUDE.md "Reduce Proprietary API Surface").

The Python-facing surface-current ObjRecMag constructor has now been RETIRED (un-exposed from the
extension); rad.ObjRecMag is a thin shim in radia.magnet that forwards to magnet_box.  So the historical
"surface-charge magnet == surface-current ObjRecMag" cross-check is now a TAUTOLOGY -- it is kept only as
a shim-identity smoke (rad.ObjRecMag is the magnet_box shim and gives bit-identical fields).

Locks (a) magnet_box reproduces stored reference fields (ObjRecMag-INDEPENDENT, so this golden is the real
permanent-magnet correctness lock), (b) the rad.ObjRecMag shim IS magnet_box (identity smoke), (c) a PM
needs no Solve, (d) the on-axis transverse field vanishes by symmetry.
"""
from pathlib import Path
import sys
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import radia as rad

# Canonical box: 20x20x10 mm, M = 954930 A/m (Br = 1.2 T) along +z.
CENTER = [0.0, 0.0, 0.0]
DIMS = [0.02, 0.02, 0.01]
MAGN = [0.0, 0.0, 954930.0]

# Reference B (Tesla), magnet_box, captured + RecMag-cross-validated (rel ~1e-7) 2026-06-28.
REF = {
    (0.0, 0.0, 0.02): [0.0, 0.0, 0.06661383297],
    (0.015, 0.01, 0.012): [0.04986655, 0.030685837, 0.009070943],
    (0.0, 0.0, 0.1): [0.0, 0.0, 7.525758956e-04],
}


def test_magnet_box_reproduces_reference_field():
    """ObjRecMag-INDEPENDENT lock: magnet_box gives the stored permanent-magnet field (no Solve)."""
    rad.UtiDelAll()
    pm = rad.magnet_box(CENTER, DIMS, MAGN)
    for pt, b_ref in REF.items():
        b = np.array(rad.Fld(pm, "b", list(pt)))
        ref = np.array(b_ref)
        assert np.linalg.norm(b - ref) <= 1e-5 * np.linalg.norm(ref) + 1e-6, (
            f"magnet_box B{pt} = {b} != reference {ref}"
        )
    rad.UtiDelAll()


def test_objrecmag_shim_is_magnet_box():
    """Identity smoke: with the C++ ObjRecMag constructor retired, rad.ObjRecMag is the radia.magnet
    shim forwarding to magnet_box -- so it must give a BIT-IDENTICAL field.  (This replaces the old
    surface-charge-vs-surface-current cross-check, which is now a tautology.)"""
    assert hasattr(rad, "ObjRecMag"), "rad.ObjRecMag shim should remain for production/other-solver scripts"
    assert rad.ObjRecMag.__module__ == "radia.magnet", (
        f"rad.ObjRecMag should be the radia.magnet shim, not the C++ ObjRecMag "
        f"(got {rad.ObjRecMag.__module__})"
    )
    # C++ ObjRecMag constructor must be un-exposed from the extension (retirement is real, not shadowed).
    import radia._radia_pybind as _p
    assert not hasattr(_p, "ObjRecMag"), "C++ ObjRecMag is still exposed in _radia_pybind"

    rad.UtiDelAll()
    pm = rad.magnet_box(CENTER, DIMS, MAGN)
    pts = list(REF.keys()) + [(0.03, -0.02, 0.04)]
    b_box = [np.array(rad.Fld(pm, "b", list(p))) for p in pts]
    rad.UtiDelAll()
    rm = rad.ObjRecMag(CENTER, DIMS, MAGN)   # shim -> magnet_box
    b_shim = [np.array(rad.Fld(rm, "b", list(p))) for p in pts]
    rad.UtiDelAll()
    for p, a, b in zip(pts, b_box, b_shim):
        assert np.array_equal(a, b), f"rad.ObjRecMag shim != magnet_box at {p}: {a} vs {b}"


def test_magnet_box_no_solve_needed():
    """A permanent magnet magnetizes without Solve (fixed M; sigma==0 -> analytic field branch)."""
    rad.UtiDelAll()
    pm = rad.magnet_box(CENTER, DIMS, MAGN)
    bz = rad.Fld(pm, "b", [0.0, 0.0, 0.02])[2]
    rad.UtiDelAll()
    assert bz > 0.05, f"on-axis Bz above a +z magnet should be ~0.067 T, got {bz}"


def test_magnet_box_on_axis_symmetry():
    """On the magnetization axis the transverse field vanishes by symmetry (surface-charge model OK)."""
    rad.UtiDelAll()
    pm = rad.magnet_box(CENTER, DIMS, MAGN)
    b = rad.Fld(pm, "b", [0.0, 0.0, 0.05])
    rad.UtiDelAll()
    assert abs(b[0]) < 1e-9 and abs(b[1]) < 1e-9, f"transverse field should vanish on axis: {b}"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
