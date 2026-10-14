"""saddle_coil(): the cylinder-wrapped DC end-turn geometry as an API.

Promotes the tilt-search recipe validated on 2026-08-06/07 (the demo
that also exposed the to_occ_shape Euler-order bug) into
radia.coil_builder.saddle_coil.  The locks:

* the centerline CLOSES exactly and stays on the cylinder except for
  the fillets, whose outward excursion is the geometric
  bend_radius^2/(2 radius) flare of a round corner leaving a cylinder;
* a two-saddle dipole reproduces the independently measured golden
  field B(0) = 0.0745 T (5000 A reference geometry), which itself was
  cross-checked against a filament Biot-Savart integration to 2.7e-3;
* the axis parameter rotates the whole construction consistently;
* geometric nonsense fails loudly instead of building something else.
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from radia.coil_builder import saddle_coil  # noqa: E402

MM = 1e-3
REF = dict(radius=40 * MM, length=120 * MM, span_deg=70.0,
           bend_radius=6 * MM, width=10 * MM, height=8 * MM)


def _path(coil, n_arc=64):
    segs, _ = coil.to_wire_segments(n_arc=n_arc)
    return np.array([s[0] for s in segs] + [segs[-1][1]])


def test_saddle_closes_and_stays_on_the_cylinder():
    coil = saddle_coil(1000.0, **REF)
    assert len(coil.segments) == 8            # 2 wraps + 4 fillets + 2 legs
    p = _path(coil)
    assert float(np.linalg.norm(p[-1] - p[0])) < 1e-9

    r = np.hypot(p[:, 0], p[:, 2])            # axis = y
    assert r.min() == pytest.approx(REF["radius"], abs=1e-9)
    flare = REF["bend_radius"] ** 2 / (2 * REF["radius"])
    assert r.max() - REF["radius"] < 1.05 * flare
    # the axial window is respected
    assert abs(p[:, 1]).max() == pytest.approx(REF["length"] / 2, abs=1e-9)


def test_saddle_pair_reproduces_the_golden_dipole_field():
    rad = pytest.importorskip("radia", reason="radia not installed")

    rad.UtiDelAll()
    try:
        objs = []
        for phi, sign in ((0.0, +1.0), (180.0, -1.0)):
            coil = saddle_coil(sign * 5000.0, phi_center_deg=phi, **REF)
            objs.extend(coil.to_radia(arc_max_segment_length=4 * MM))
        src = rad.ObjCnt(objs)
        B = np.array(rad.Fld(src, 'b', [0.0, 0.0, 0.0]))
    finally:
        rad.UtiDelAll()

    # golden band from the 2026-08-07 demo (cross-checked against an
    # independent filament Biot-Savart integration, rel 2.7e-3)
    assert B[0] == pytest.approx(0.07451, rel=5e-3)
    assert abs(B[1]) < 1e-6 * abs(B[0])
    assert abs(B[2]) < 1e-6 * abs(B[0])


def test_axis_parameter_rotates_the_construction():
    for axis, ax_idx in (("x", 0), ("y", 1), ("z", 2)):
        coil = saddle_coil(1000.0, axis=axis, **REF)
        p = _path(coil)
        perp = [i for i in range(3) if i != ax_idx]
        r = np.hypot(p[:, perp[0]], p[:, perp[1]])
        assert r.min() == pytest.approx(REF["radius"], abs=1e-9), axis
        assert abs(p[:, ax_idx]).max() == pytest.approx(
            REF["length"] / 2, abs=1e-9), axis
        assert float(np.linalg.norm(p[-1] - p[0])) < 1e-9, axis


def test_saddle_rejects_geometric_nonsense():
    with pytest.raises(ValueError, match="axis"):
        saddle_coil(1.0, axis="w", **REF)
    bad = dict(REF)
    bad["span_deg"] = 200.0
    with pytest.raises(ValueError, match="span_deg"):
        saddle_coil(1.0, **bad)
    bad = dict(REF)
    bad["bend_radius"] = 50 * MM              # > radius
    with pytest.raises(ValueError, match="bend_radius"):
        saddle_coil(1.0, **bad)
    bad = dict(REF)
    bad["length"] = 10 * MM                   # < 2 * bend_radius
    with pytest.raises(ValueError, match="length"):
        saddle_coil(1.0, **bad)


def test_current_sign_flips_the_field():
    rad = pytest.importorskip("radia", reason="radia not installed")

    rad.UtiDelAll()
    try:
        Bs = []
        for sign in (+1.0, -1.0):
            objs = saddle_coil(sign * 2000.0, **REF).to_radia(
                arc_max_segment_length=5 * MM)
            src = rad.ObjCnt(objs)
            Bs.append(np.array(rad.Fld(src, 'b', [0.0, 0.0, 0.0])))
            rad.UtiDelAll()
    finally:
        rad.UtiDelAll()
    assert np.allclose(Bs[0], -Bs[1], rtol=1e-12, atol=1e-15)
