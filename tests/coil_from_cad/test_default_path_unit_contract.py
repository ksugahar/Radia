"""Unit-contract + partial-walk locks for the DEFAULT filament path.

Background (2026-07-29): ``_filaments_via_coil_builder`` (the default
``filaments_from_step`` path, ``n_peri=None``) never scaled the walked
centerline by ``1 / cad_units_per_meter`` -- a mm-authored STEP flowed
into PEEC as if metres and returned a silently ~1000x-too-large L
(measured: 90.8 nH vs 90809 nH on an identical 30 mm torus authored in
m vs mm).  Every earlier test of this path authored geometry in metres
with ``cad_units_per_meter=1``, which hid the missing unit boundary.
The production CLI (``calc_inductance`` / ``calc_fem_kelvin``) always
passes ``n_peri`` and therefore routed through the correctly-scaled UV
tiers; only direct API callers of the default path were affected.

The same path also returned a PARTIAL coil silently when the walker
halted mid-solid (observed on a rounded-rectangle sweep: a single
~97 mm leg of a ~600 mm racetrack) because it was the only filament
tier without ``_check_filaments_cover_solid_bbox``.

These tests pin both fixes:
  1. mm/m twin STEP -> identical metre-level filaments and L;
  2. truncated walk -> hard ValueError (tier='coil_builder');
  3. ``filaments_from_shape`` (in-memory entry) shares the same unit
     boundary, orientation-agnostic seed, and coverage check.
"""

import numpy as np
import pytest


def _make_torus_step(path, scale):
    """Write a full-revolution circular torus coil STEP.

    R = 30 mm major, r = 3 mm minor, authored at ``scale`` CAD units
    per metre (1.0 -> metre coordinates, 1000.0 -> mm coordinates).
    """
    from netgen.occ import Axis, Axes, Dir, Pnt, WorkPlane
    R, r = 0.030 * scale, 0.003 * scale
    profile = WorkPlane(Axes(
        p=Pnt(R, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1),
    )).Circle(r).Face()
    coil = profile.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360)
    coil.WriteStep(str(path))


def _flat_points(topo):
    paths = topo.get("filament_paths") or []
    return np.asarray(
        [p for path in paths for seg in path for p in seg], dtype=float)


def _port_L_nH(topo):
    Z = topo["solver"].compute_port_impedance(1.0)
    return float(np.imag(Z) / (2.0 * np.pi * 1.0) * 1e9)


# Analytic single circular loop: L = mu0 * R * (ln(8R/r) - 2)
#   R = 30 mm, r = 3 mm  ->  89.8 nH.
_L_ANALYTIC_NH = 89.8


def test_mm_and_m_twin_torus_identical(tmp_path, monkeypatch):
    """A mm-authored STEP with cad_units_per_meter=1000 must give the
    same metre-level filaments and L as its metre-authored twin."""
    from radia.coil_from_cad import filaments_from_step

    monkeypatch.setenv("RADIA_PEEC_CACHE_DISABLE", "1")
    f_m = tmp_path / "torus_m.step"
    f_mm = tmp_path / "torus_mm.step"
    _make_torus_step(f_m, 1.0)
    _make_torus_step(f_mm, 1000.0)

    topo_m = filaments_from_step(str(f_m), cad_units_per_meter=1.0)
    topo_mm = filaments_from_step(str(f_mm), cad_units_per_meter=1000.0)

    pts_m = _flat_points(topo_m)
    pts_mm = _flat_points(topo_mm)
    # Both in metres: the coil must fit inside a 0.1 m ball.  (Broken
    # unit boundary left the mm case at ~30 "units" = raw mm.)
    assert float(np.max(np.abs(pts_m))) < 0.05
    assert float(np.max(np.abs(pts_mm))) < 0.05

    L_m = _port_L_nH(topo_m)
    L_mm = _port_L_nH(topo_mm)
    # The two walks differ only by OCC numerics at different absolute
    # scales; loose 0.1% tolerance absorbs that, while the pre-fix bug
    # was a factor of 1000.
    assert L_mm == pytest.approx(L_m, rel=1e-3)
    assert L_m == pytest.approx(_L_ANALYTIC_NH, rel=0.10)


def test_partial_walk_fails_loud(tmp_path, monkeypatch):
    """A walker halt mid-solid (partial centerline) must raise, not
    silently feed a partial coil to PEEC."""
    import radia.coil_from_step as cfs
    from radia.coil_from_cad import filaments_from_step

    monkeypatch.setenv("RADIA_PEEC_CACHE_DISABLE", "1")
    f_m = tmp_path / "torus_partial.step"
    _make_torus_step(f_m, 1.0)

    real_extract = cfs.extract_centerline

    def truncated(*args, **kwargs):
        res = real_extract(*args, **kwargs)
        k = max(4, len(res.polyline) // 4)      # keep a quarter arc
        return cfs.CenterlineResult(
            polyline=res.polyline[:k],
            tangents=res.tangents[:k],
            profiles=res.profiles[:k],
            polygons=res.polygons[:k],
            arclen=res.arclen[:k],
            closed=False,
        )

    monkeypatch.setattr(cfs, "extract_centerline", truncated)
    with pytest.raises(ValueError, match="coil_builder"):
        filaments_from_step(str(f_m), cad_units_per_meter=1.0)


def test_filaments_from_shape_oblique_mm(monkeypatch):
    """In-memory build123d entry: a mm-authored oblique torus with
    cad_units_per_meter=1000 must yield metre filaments and a plausible
    L (locks the orientation-agnostic seed + unit boundary + coverage
    check inside ``filaments_from_shape``)."""
    bd = pytest.importorskip("build123d")
    from radia.coil_from_cad import filaments_from_shape

    monkeypatch.setenv("RADIA_PEEC_CACHE_DISABLE", "1")
    solid = bd.Solid.make_torus(30.0, 3.0)     # mm-authored
    solid = solid.rotate(bd.Axis.X, 37).rotate(bd.Axis.Y, 23)

    topo = filaments_from_shape(solid, cad_units_per_meter=1000.0)

    pts = _flat_points(topo)
    assert pts.shape[0] > 20
    assert np.all(np.isfinite(pts))
    assert float(np.max(np.abs(pts))) < 0.05   # metres, not raw mm
    assert _port_L_nH(topo) == pytest.approx(_L_ANALYTIC_NH, rel=0.10)
