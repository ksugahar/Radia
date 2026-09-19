"""Swept fin graph on synthetic stations (no CAD)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from radia.fin_section import analyze_section, resample_outline
from radia.fin_sweep import (
    StationOutline,
    _is_straight_prism,
    align_outline_origin,
    align_station_origins,
    build_fin_graph,
    parallel_transport_frames,
    probe_points_3d,
    station_outlines_from_step,
)


def _beak_outline(scale=1.0, n_arc=48):
    """Beak fixture section with the fin scaled in x by ``scale`` (mm -> m)."""
    tl = (3.8, -(0.25**2 - 0.05**2) ** 0.5)
    tu = (3.8, +(0.25**2 - 0.05**2) ** 0.5)
    centre, radius = (3.75, 0.0), 0.25
    a0 = math.atan2(tl[1], tl[0] - centre[0])
    a1 = math.atan2(tu[1], tu[0] - centre[0])
    pts = [(-4, -2), (1.6, -2), (1.6, -0.7), tl]
    pts += [(centre[0] + radius * math.cos(t), radius * math.sin(t))
            for t in np.linspace(a0, a1, n_arc)[1:-1]]
    pts += [tu, (1.6, 0.7), (1.6, 2), (-4, 2)]
    p = np.array(pts, dtype=float)
    fin = p[:, 0] > 1.6
    p[fin, 0] = 1.6 + (p[fin, 0] - 1.6) * scale
    if scale < 0.05:
        p = np.array([(-4, -2), (1.6, -2), (1.6, 2), (-4, 2)], dtype=float)
    return p * 1e-3


def _stations(centroids, scales, n_outline=1024, roll=None):
    t, u, v = parallel_transport_frames(centroids)
    out = []
    for i, (c, sc) in enumerate(zip(centroids, scales)):
        uv = resample_outline(_beak_outline(sc), n_outline)
        if roll is not None:
            uv = np.roll(uv, roll[i], axis=0)
        per = float(np.sum(np.linalg.norm(np.roll(uv, -1, axis=0) - uv, axis=1)))
        area = 0.5 * abs(np.sum(uv[:, 0] * np.roll(uv[:, 1], -1)
                                - np.roll(uv[:, 0], -1) * uv[:, 1]))
        out.append(StationOutline(i, np.asarray(c, float), u[i], v[i], t[i],
                                  uv, area, per))
    return tuple(out)


def test_parallel_transport_frames_are_orthonormal_and_smooth():
    theta = np.linspace(0, np.pi / 2, 12)
    c = np.c_[0.03 * np.cos(theta), 0.03 * np.sin(theta), 0.002 * theta]
    t, u, v = parallel_transport_frames(c)
    for i in range(len(c)):
        assert abs(t[i] @ u[i]) < 1e-9 and abs(t[i] @ v[i]) < 1e-9
        assert abs(u[i] @ v[i]) < 1e-9
        assert np.allclose([np.linalg.norm(t[i]), np.linalg.norm(u[i]),
                            np.linalg.norm(v[i])], 1.0)
    assert np.all(np.einsum("ij,ij->i", u[1:], u[:-1]) > 0.98)


def test_origin_registration_recovers_a_cyclic_shift():
    ref = resample_outline(_beak_outline(1.0), 1024)
    shifted = np.roll(ref, 137, axis=0)
    aligned, shift = align_outline_origin(ref, shifted)
    assert shift == (1024 - 137) % 1024 or np.allclose(aligned, ref)
    assert np.allclose(aligned, ref)
    # slightly different section (fin shortened) still registers near zero
    other = np.roll(resample_outline(_beak_outline(0.7), 1024), 300, axis=0)
    aligned, _ = align_outline_origin(ref, other)
    assert np.linalg.norm(aligned[0] - ref[0]) < 0.3e-3


def test_straight_tapering_fin_sweep_builds_conservative_graph():
    n_st = 9
    z = np.linspace(0, 0.024, n_st)
    centroids = np.c_[np.zeros(n_st), np.zeros(n_st), z]
    scales = np.linspace(1.0, 0.0, n_st)            # fin vanishes at the end
    rng = np.random.default_rng(0)
    raw = _stations(centroids, scales, roll=rng.integers(0, 1024, n_st))
    stations = align_station_origins(raw)
    result = build_fin_graph(stations, n_lanes=64, lane_grading="auto")
    g = result.graph
    assert g.n_stations == n_st and g.n_lanes == 64
    assert g.mesh_stations == tuple(range(1, n_st - 1))
    # fin detected while it is long enough, absent once it is gone
    fins = [f.fin is not None for f in result.features]
    assert fins[0] and fins[1] and fins[2]
    assert not fins[-1]
    assert result.meta["fin_stations"] == [i for i, f in enumerate(fins) if f]
    # lanes follow material lines: longitudinal branches stay short
    spacing = z[1] - z[0]
    lon = np.linalg.norm(result.rings[1:] - result.rings[:-1], axis=2)
    assert np.all(lon < 1.6 * spacing)
    # graded where a fin exists, uniform elsewhere
    assert result.features[0].graded and not result.features[-1].graded
    # fin lanes shrink with the fin: current-carrying tip lanes vanish
    tip0 = result.features[0].fin.tip_mask(result.features[0].lane_uv).sum()
    assert tip0 >= 8


def test_curved_sweep_lifts_rings_with_the_frame():
    n_st = 10
    theta = np.linspace(0, np.pi / 2, n_st)
    r0 = 0.040
    centroids = np.c_[r0 * np.cos(theta), r0 * np.sin(theta), np.zeros(n_st)]
    stations = align_station_origins(_stations(centroids, np.ones(n_st)))
    result = build_fin_graph(stations, n_lanes=48, lane_grading="auto")
    rings = result.rings
    # every ring lies in its station plane spanned by (u, v)
    for i, st in enumerate(stations):
        off = rings[i] - st.centroid
        normal = np.cross(st.u_hat, st.v_hat)
        assert np.max(np.abs(off @ normal)) < 1e-9
    # planar coil: u = +z, so the section's v (former y) is radial
    assert np.allclose(stations[0].u_hat, [0, 0, 1])
    # the section's +u (fin axis) maps to +z on a planar coil
    f0 = result.features[0]
    ext = stations[0].to_xyz(f0.fin.extremity[None, :])[0]
    assert ext[2] > stations[0].centroid[2] + 3.9e-3
    assert len(result.graph.branches) == (n_st - 1) * 48 + (n_st - 2) * 48
    probe = probe_points_3d(f0, stations[0], offset=1e-3, n=5)
    assert probe.shape == (5, 3)
    normal0 = np.cross(stations[0].u_hat, stations[0].v_hat)
    assert np.max(np.abs((probe - stations[0].centroid) @ normal0)) < 1e-9


def test_build_requires_three_stations_and_known_grading():
    c = np.c_[np.zeros(2), np.zeros(2), [0.0, 0.01]]
    with pytest.raises(ValueError):
        build_fin_graph(_stations(c, np.ones(2)), n_lanes=32)
    c = np.c_[np.zeros(3), np.zeros(3), [0.0, 0.01, 0.02]]
    with pytest.raises(ValueError):
        build_fin_graph(_stations(c, np.ones(3)), n_lanes=32, lane_grading="x")


def test_straight_step_adapter_preserves_wire_order_and_scale():
    fixture = (Path(__file__).parent / "coil_from_cad" / "fixtures"
               / "beak_fin_short.step")
    stations, tag = station_outlines_from_step(
        fixture, n_stations=3, n_outline=256, cad_units_per_meter=1000)
    assert tag == "step_straight_prism_sections"
    assert len(stations) == 3
    assert np.allclose([st.area for st in stations], 24.55223e-6, rtol=2e-5)
    for st in stations:
        edge_lengths = np.linalg.norm(np.roll(st.uv, -1, axis=0) - st.uv,
                                      axis=1)
        assert np.max(edge_lengths) < 2.0 * st.perimeter / len(st.uv)
        assert analyze_section(st.uv).primary is not None


def test_flat_planar_coil_is_not_classified_as_axial_z_sweep():
    from radia._b3d_shim import import_step

    fixture = (Path(__file__).parents[1] / "validation_test" / "panels"
               / "golden" / "rect_torus_lofted_united.step")
    solid = import_step(str(fixture)).solids()[0]
    assert not _is_straight_prism(solid)


# ----------------------------------------------------------------------
# automatic CAD unit + discretisation
# ----------------------------------------------------------------------

FIXTURES = Path(__file__).parent / "coil_from_cad" / "fixtures"
GOLDEN = (Path(__file__).parents[1] / "validation_test" / "panels" / "golden")


def test_cad_unit_scale_is_resolved_for_both_lab_numbering_conventions():
    """Every lab STEP declares millimetre in its header whatever its
    numbers mean, so only the extent can decide -- and it must."""
    from radia._b3d_shim import import_step
    from radia.fin_sweep import resolve_cad_units_per_meter

    millimetre_numbered = import_step(
        str(FIXTURES / "beak_fin_straight.step")).solids()[0]
    metre_numbered = import_step(
        str(GOLDEN / "rect_torus_lofted_united.step")).solids()[0]

    assert resolve_cad_units_per_meter(millimetre_numbered) == 1000.0
    assert resolve_cad_units_per_meter(metre_numbered) == 1.0


def test_cad_unit_scale_refuses_an_ambiguous_extent():
    """A conductor of a few CAD units is plausible at BOTH scales; the
    resolver must say so instead of silently picking one."""
    from radia.fin_sweep import resolve_cad_units_per_meter

    class _Corner:
        def __init__(self, v):
            self.X = self.Y = self.Z = v

    class _Box:
        min = _Corner(0.0)
        max = _Corner(2.0)

    class _Solid:
        def bounding_box(self):
            return _Box()

    with pytest.raises(ValueError, match="ambiguous"):
        resolve_cad_units_per_meter(_Solid())


def test_auto_resolution_measures_the_beak_tip_and_sizes_the_lanes():
    from radia.fin_sweep import (
        LANE_SPACING_PER_TIP_RADIUS,
        OUTLINE_SAMPLES_PER_LANE,
        auto_fin_resolution,
    )

    auto = auto_fin_resolution(FIXTURES / "beak_fin_straight.step")
    measured = auto["measured"]

    # Fixture geometry recovered from CAD alone, with no unit argument.
    assert measured["tip_radius_m_min"] == pytest.approx(0.25e-3, abs=1e-6)
    assert measured["sweep_length_m"] == pytest.approx(60e-3, rel=1e-6)
    assert measured["fin_probe_stations"] == measured["probe_stations"]

    # Lanes follow the measured tip radius, outline follows the lanes.
    spacing = LANE_SPACING_PER_TIP_RADIUS * measured["tip_radius_m_min"]
    assert auto["n_lanes"] == pytest.approx(
        measured["perimeter_m"] / spacing, rel=0.1)
    assert auto["n_outline"] >= OUTLINE_SAMPLES_PER_LANE * auto["n_lanes"]
    assert auto["lane_grading"] == "auto"
    assert auto["rule"]["lane_spacing_driver"] == "fin tip radius"


def test_fin_graph_from_step_needs_only_the_file_and_records_its_choice():
    from radia.fin_sweep import fin_graph_from_step

    sweep = fin_graph_from_step(FIXTURES / "beak_fin_short.step")
    resolution = sweep.meta["resolution"]

    assert sweep.meta["auto_resolution"] is not None
    assert sweep.meta["n_lanes"] == resolution["n_lanes"]
    assert sweep.meta["n_stations"] == resolution["n_stations"]
    # A prismatic beak fin has its fin at every station.
    assert len(sweep.fin_station_indices()) == resolution["n_stations"]
