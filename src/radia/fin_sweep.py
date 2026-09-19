"""STEP -> fin-capable surface PEEC graph for swept conductors.

This is the general-geometry successor of the straight-prism fixture route
in :mod:`radia.peec_fin_topology`.  It takes a swept conductor (any spine,
section allowed to change along the sweep, a fin that may appear, taper or
vanish) and produces the hybrid lane/ring graph that
:func:`radia.peec_fin_topology.assemble_experimental_fin_peec` consumes.

Pipeline
--------
1. ``station_outlines_from_step`` — dense (default 2048-point) outline of
   the cross-section at every station in the station's parallel-transport
   frame, using the existing CAD extractors (``_section_faces_from_solid``
   for general sweeps, direct z-sectioning for straight prisms).  No PEEC
   solver is built here.
2. ``align_station_origins`` — the CAD sampler starts each outline at the
   boundary point nearest +u; when a fin appears that start point jumps.
   Each station is cyclically re-registered to its predecessor so lane k
   follows the same material line along the sweep.
3. ``build_fin_graph`` — every station is analysed with
   :func:`radia.fin_section.analyze_section`; lanes are graded toward the
   detected fin tip (``lane_grading="auto"``) or uniform; the ring points
   are lifted back to 3-D with the station frame; circumferential rings are
   placed at every interior station so longitudinal redistribution is
   available wherever the section changes.

Everything except the STEP adapters is pure numpy and unit-tested on
synthetic stations (tapering fin, curved sweep).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from radia.fin_section import (
    SectionAnalysis,
    _as_closed_ccw,
    analyze_section,
    graded_lane_arclengths,
    outline_arclength,
    outline_points_at,
)

# ----------------------------------------------------------------------
# station container
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class StationOutline:
    """Dense cross-section outline at one spine station (SI units)."""
    index: int
    centroid: np.ndarray          # (3,) m
    u_hat: np.ndarray             # (3,) in-plane axis
    v_hat: np.ndarray             # (3,) in-plane axis, n x u
    tangent: np.ndarray           # (3,) sweep direction
    uv: np.ndarray                # (n, 2) m, CCW in (u, v)
    area: float
    perimeter: float

    def to_xyz(self, uv):
        q = np.asarray(uv, dtype=float)
        return (self.centroid[None, :] + q[:, :1] * self.u_hat[None, :]
                + q[:, 1:2] * self.v_hat[None, :])


@dataclass(frozen=True)
class StationFeatures:
    index: int
    analysis: SectionAnalysis
    lane_fraction: np.ndarray     # (n_lanes,) normalised arc length
    lane_uv: np.ndarray           # (n_lanes, 2)
    graded: bool

    @property
    def fin(self):
        return self.analysis.primary


@dataclass(frozen=True)
class FinSweepGraph:
    graph: object                 # HybridSurfaceTopology
    rings: np.ndarray             # (S, n_lanes, 3)
    stations: tuple[StationOutline, ...]
    features: tuple[StationFeatures, ...]
    meta: dict = field(default_factory=dict)

    def fin_station_indices(self):
        return tuple(f.index for f in self.features if f.fin is not None)


# ----------------------------------------------------------------------
# frames
# ----------------------------------------------------------------------

def parallel_transport_frames(centroids):
    """Tangent / u / v frames along a polyline, transported without twist.

    Planar (constant-z) spines get the axial-radial frame used by the
    existing filament builder so curved planar coils keep u = +z.
    """
    c = np.asarray(centroids, dtype=float)
    n = len(c)
    if n < 2:
        raise ValueError("need at least two stations")
    t = np.zeros_like(c)
    t[0] = c[1] - c[0]
    t[-1] = c[-1] - c[-2]
    if n > 2:
        t[1:-1] = c[2:] - c[:-2]
    t /= np.maximum(np.linalg.norm(t, axis=1), 1e-300)[:, None]

    z_spread = float(np.ptp(c[:, 2]))
    radial_min = float(np.min(np.linalg.norm(c[:, :2], axis=1)))
    if z_spread < 1e-3 * max(radial_min, 1e-30):
        u = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
        v = c.copy()
        v[:, 2] = 0.0
        v /= np.maximum(np.linalg.norm(v, axis=1), 1e-300)[:, None]
        return t, u, v

    u = np.zeros_like(c)
    seed = np.array([0.0, 0.0, 1.0])
    if abs(seed @ t[0]) > 0.9:
        seed = np.array([1.0, 0.0, 0.0])
    u0 = seed - (seed @ t[0]) * t[0]
    u[0] = u0 / np.linalg.norm(u0)
    for i in range(1, n):
        # rotate previous u by the minimal rotation taking t[i-1] to t[i]
        axis = np.cross(t[i - 1], t[i])
        s = np.linalg.norm(axis)
        d = float(np.clip(t[i - 1] @ t[i], -1.0, 1.0))
        if s < 1e-12:
            ui = u[i - 1]
        else:
            axis /= s
            ang = np.arctan2(s, d)
            ui = (u[i - 1] * np.cos(ang) + np.cross(axis, u[i - 1]) * np.sin(ang)
                  + axis * (axis @ u[i - 1]) * (1 - np.cos(ang)))
        ui -= (ui @ t[i]) * t[i]
        u[i] = ui / np.linalg.norm(ui)
    v = np.cross(t, u)
    return t, u, v


# ----------------------------------------------------------------------
# origin registration
# ----------------------------------------------------------------------

def align_outline_origin(reference_uv, uv):
    """Cyclically shift ``uv`` (same sample count) to best match ``reference_uv``.

    Both outlines are dense equal-arc-length CCW samplings.  The shift that
    minimises the mean squared point distance is found with an FFT
    cross-correlation over all integer shifts.  Returns ``(uv_aligned,
    shift)`` where ``uv_aligned = np.roll(uv, -shift, axis=0)``.
    """
    a = np.asarray(reference_uv, dtype=float)
    b = np.asarray(uv, dtype=float)
    if a.shape != b.shape:
        raise ValueError("outlines must have the same sample count")
    n = len(a)
    # sum_i |a_i - b_{i+k}|^2 = const - 2 sum_i a_i . b_{i+k}
    corr = np.zeros(n)
    for dim in range(2):
        fa = np.fft.rfft(a[:, dim])
        fb = np.fft.rfft(b[:, dim])
        corr += np.fft.irfft(np.conj(fa) * fb, n)
    shift = int(np.argmax(corr))
    return np.roll(b, -shift, axis=0), shift


def align_station_origins(stations):
    """Return stations re-registered so arc-length origins are continuous."""
    if not stations:
        return ()
    out = [stations[0]]
    for st in stations[1:]:
        uv_aligned, _shift = align_outline_origin(out[-1].uv, st.uv)
        out.append(StationOutline(st.index, st.centroid, st.u_hat, st.v_hat,
                                  st.tangent, uv_aligned, st.area,
                                  st.perimeter))
    return tuple(out)


# ----------------------------------------------------------------------
# graph construction
# ----------------------------------------------------------------------

def _lane_fraction_for_station(analysis, n_lanes, lane_grading, tip_lanes):
    if lane_grading == "auto" and analysis.fins:
        s = graded_lane_arclengths(analysis, n_lanes, tip_lanes=tip_lanes)
        return s / analysis.perimeter, True
    return (np.arange(n_lanes) + 0.5) / n_lanes, False


def build_fin_graph(stations, *, n_lanes=64, lane_grading="auto",
                    tip_lanes=8, mesh_stations="all", fin_kwargs=None):
    """Analyse every station and build the hybrid lane/ring PEEC graph.

    Args:
        stations: sequence of :class:`StationOutline` (already origin
            aligned; see :func:`align_station_origins`).
        n_lanes: longitudinal lanes (constant along the sweep).
        lane_grading: ``"auto"`` grades toward detected fin tips per station,
            ``"uniform"`` keeps equal arc length everywhere.
        tip_lanes: cells on each tip arc when grading.
        mesh_stations: ``"all"`` (every interior station gets a ring) or an
            explicit iterable of interior station indices.
        fin_kwargs: forwarded to :func:`analyze_section`.
    """
    from radia.peec_fin_topology import build_hybrid_surface_topology

    if lane_grading not in {"uniform", "auto"}:
        raise ValueError("lane_grading must be 'uniform' or 'auto'")
    stations = tuple(stations)
    if len(stations) < 3:
        raise ValueError("need at least three stations")
    fin_kwargs = dict(fin_kwargs or {})
    features = []
    rings = np.empty((len(stations), n_lanes, 3))
    for i, st in enumerate(stations):
        closed = _as_closed_ccw(st.uv)
        analysis = analyze_section(closed, **fin_kwargs)
        frac, graded = _lane_fraction_for_station(
            analysis, n_lanes, lane_grading, tip_lanes)
        if np.any(np.diff(frac) <= 0):
            raise ValueError(f"station {st.index}: lane fractions not increasing")
        _s, per = outline_arclength(closed)
        lane_uv = outline_points_at(closed, frac * per)
        rings[i] = st.to_xyz(lane_uv)
        features.append(StationFeatures(st.index, analysis, frac, lane_uv, graded))
    if mesh_stations == "all":
        mesh = range(1, len(stations) - 1)
    else:
        mesh = tuple(int(m) for m in mesh_stations)
    graph = build_hybrid_surface_topology(rings, mesh)
    fin_idx = [f.index for f in features if f.fin is not None]
    meta = {
        "n_stations": len(stations), "n_lanes": n_lanes,
        "lane_grading": lane_grading, "tip_lanes": tip_lanes,
        "mesh_stations": tuple(graph.mesh_stations),
        "fin_stations": fin_idx,
        "fin_fraction_of_sweep": len(fin_idx) / len(stations),
        "station_fins": [
            None if f.fin is None else f.fin.as_dict() for f in features],
        "n_branches": len(graph.branches),
    }
    return FinSweepGraph(graph, rings, stations, tuple(features), meta)


def probe_points_3d(feature, station, offset, n=41):
    """Workpiece-side probe segment for one fin station, in 3-D."""
    fin = feature.fin
    if fin is None:
        raise ValueError("station has no fin")
    return station.to_xyz(fin.probe_points(offset=offset, n=n))


# ----------------------------------------------------------------------
# CAD unit resolution
# ----------------------------------------------------------------------

#: The two numeric conventions actually present in the tracked STEP
#: corpus.  Every lab STEP declares ``SI_UNIT(.MILLI.,.METRE.)`` in its
#: header regardless of what its numbers mean, so the header carries no
#: information: ``rect_torus_lofted_united`` and ``keiko_outsideline``
#: hold METRE-valued numbers while ``beak_fin_*`` and
#: ``tapered_helix_coil`` hold MILLIMETRE-valued numbers (measured
#: 2026-09-20).  Only the magnitude separates them.
CAD_UNIT_CANDIDATES = (1.0, 1000.0)

#: Conductor extent window used to decide between the candidates.  The
#: upper bound matches ``coil_from_cad``'s plausibility guard; the lower
#: bound rejects a sub-millimetre "coil".
_FIN_MIN_EXTENT_M = 1.0e-3
_FIN_MAX_EXTENT_M = 5.0


def resolve_cad_units_per_meter(solid, source_tag="fin_sweep"):
    """Scale for ``solid`` when the STEP determines it UNIQUELY.

    Tries each entry of :data:`CAD_UNIT_CANDIDATES` and keeps the ones
    that put the solid's largest extent inside the conductor window.
    Exactly one survivor is the answer; zero or two is reported as an
    error naming both, because guessing here silently scales every
    downstream length (and inductance) by 1000.  Pass an explicit
    ``cad_units_per_meter`` to settle an ambiguous case.
    """
    box = solid.bounding_box()
    extent_cad = max(float(box.max.X - box.min.X),
                     float(box.max.Y - box.min.Y),
                     float(box.max.Z - box.min.Z))
    if not np.isfinite(extent_cad) or extent_cad <= 0:
        raise ValueError(
            f"{source_tag}: STEP bounding box is degenerate "
            f"(largest extent {extent_cad!r} CAD units).")
    viable = [c for c in CAD_UNIT_CANDIDATES
              if _FIN_MIN_EXTENT_M <= extent_cad / c <= _FIN_MAX_EXTENT_M]
    if len(viable) == 1:
        return viable[0]
    implied = ", ".join(f"{c:g} -> {extent_cad / c:.4g} m"
                        for c in CAD_UNIT_CANDIDATES)
    if not viable:
        raise ValueError(
            f"{source_tag}: no CAD unit candidate puts the conductor in "
            f"the {_FIN_MIN_EXTENT_M * 1e3:g} mm .. {_FIN_MAX_EXTENT_M:g} m "
            f"window (largest extent {extent_cad:.6g} CAD units; {implied}). "
            f"Pass cad_units_per_meter explicitly.")
    raise ValueError(
        f"{source_tag}: the CAD unit is ambiguous -- {implied} are both "
        f"plausible for a conductor of {extent_cad:.6g} CAD units. "
        f"Pass cad_units_per_meter explicitly (the STEP header says "
        f"millimetre for both lab conventions, so it cannot decide).")


# ----------------------------------------------------------------------
# STEP adapters (need build123d)
# ----------------------------------------------------------------------

def _is_straight_prism(solid, tol=1e-6):
    """True for an axial z extrusion/loft, rather than a flat planar coil.

    Equal cap area is deliberately not required: a straight tapered sweep is
    still best sectioned directly by z planes.  The axial-aspect check keeps a
    flat z-extruded torus from being mistaken for a conductor swept along z.
    """
    from radia._b3d_shim import GeomType

    caps = []
    for face in solid.faces():
        if face.geom_type != GeomType.PLANE:
            continue
        c = face.center()
        nrm = face.normal_at(c)
        if abs(float(nrm.Z)) > 1.0 - tol:
            caps.append(face)
    if len(caps) != 2:
        return False
    bbox = solid.bounding_box()
    z_span = float(bbox.max.Z - bbox.min.Z)
    mean_cap_area = 0.5 * sum(float(face.area) for face in caps)
    return z_span * z_span >= 0.25 * mean_cap_area


def _stations_from_straight_prism(solid, n_stations, n_outline,
                                  cad_units_per_meter):
    from build123d import Plane, section

    from radia.coil_from_cad import _sample_face_perimeter_in_pt_frame

    bbox = solid.bounding_box()
    z0, z1 = float(bbox.min.Z), float(bbox.max.Z)
    if z1 <= z0:
        raise ValueError("zero-length STEP solid")
    e_u = np.array([1.0, 0.0, 0.0])
    e_v = np.array([0.0, 1.0, 0.0])
    e_t = np.array([0.0, 0.0, 1.0])
    stations = []
    for i, z in enumerate(np.linspace(z0, z1, n_stations)):
        probe_z = float(np.clip(z, z0 + 1e-7 * (z1 - z0), z1 - 1e-7 * (z1 - z0)))
        cut = section(solid, section_by=Plane(origin=(0, 0, probe_z),
                                              z_dir=(0, 0, 1)))
        faces = cut.faces()
        if len(faces) != 1:
            raise ValueError(f"station {i}: STEP section must have exactly one face")
        face = faces[0]
        c = face.center()
        origin = np.array([c.X, c.Y, c.Z])
        uv = _sample_face_perimeter_in_pt_frame(face, origin, e_u, e_v, int(n_outline))
        uv_m = uv / cad_units_per_meter
        centroid = np.array([c.X, c.Y, z]) / cad_units_per_meter
        per = float(sum(e.length for e in face.outer_wire().edges())) / cad_units_per_meter
        stations.append(StationOutline(
            i, centroid, e_u, e_v, e_t, uv_m,
            float(face.area) / cad_units_per_meter**2, per))
    return tuple(stations), "step_straight_prism_sections"


def _stations_from_section_planes(solid, n_stations, n_outline,
                                  cad_units_per_meter):
    from radia.coil_from_cad import (
        _sample_face_perimeter_in_pt_frame,
        _section_faces_from_solid,
    )

    sections = _section_faces_from_solid(solid, cad_units_per_meter, n_stations)
    if sections is None:
        raise ValueError(
            "STEP spine sectioning failed (>50 % of stations dropped); "
            "the swept-section route cannot build fin stations for this solid")
    centroids_m, faces = sections
    t, u, v = parallel_transport_frames(centroids_m)
    stations = []
    for i, face in enumerate(faces):
        c = face.center()
        origin = np.array([c.X, c.Y, c.Z])
        uv = _sample_face_perimeter_in_pt_frame(face, origin, u[i], v[i], int(n_outline))
        per = float(sum(e.length for e in face.outer_wire().edges())) / cad_units_per_meter
        stations.append(StationOutline(
            i, np.asarray(centroids_m[i], dtype=float), u[i], v[i], t[i],
            uv / cad_units_per_meter,
            float(face.area) / cad_units_per_meter**2, per))
    return tuple(stations), "step_section_planes"


def station_outlines_from_step(step_path, *, n_stations=20, n_outline=2048,
                               cad_units_per_meter="auto", route="auto"):
    """Dense per-station outlines of one STEP conductor.

    ``route``: ``"auto"`` picks the straight-prism sectioner for a single
    z-extruded prism and the spine sectioner otherwise; ``"straight_prism"``
    or ``"section_planes"`` force one.  ``cad_units_per_meter="auto"``
    defers to :func:`resolve_cad_units_per_meter`, which answers only when
    the STEP settles it.  Returns ``(stations, source_tag)`` with stations
    already origin-aligned.
    """
    from radia._b3d_shim import import_step
    from radia.coil_from_cad import _check_solid_extent_plausible

    if route not in {"auto", "straight_prism", "section_planes"}:
        raise ValueError("route must be auto, straight_prism or section_planes")
    if n_stations < 3 or n_outline < 64:
        raise ValueError("need n_stations >= 3 and n_outline >= 64")
    shape = import_step(str(step_path))
    solids = shape.solids()
    if len(solids) != 1:
        raise ValueError("fin sweep requires exactly one solid (brazed fin "
                         "plates must be united in CAD first)")
    solid = solids[0]
    if isinstance(cad_units_per_meter, str):
        if cad_units_per_meter != "auto":
            raise ValueError("cad_units_per_meter must be a number or 'auto'")
        cad_units_per_meter = resolve_cad_units_per_meter(
            solid, "station_outlines_from_step")
    cad_units_per_meter = _check_solid_extent_plausible(
        solid, cad_units_per_meter, "station_outlines_from_step")
    if route == "auto":
        route = "straight_prism" if _is_straight_prism(solid) else "section_planes"
    if route == "straight_prism":
        # The lightweight shim's section face exposes OCC edges in explorer
        # order, not guaranteed wire order. The established build123d path
        # preserves outer-wire continuity required by the perimeter sampler.
        from build123d import import_step as import_step_build123d

        authored_shape = import_step_build123d(str(step_path))
        authored_solids = authored_shape.solids()
        if len(authored_solids) != 1:
            raise ValueError("straight-prism route requires exactly one solid")
        stations, tag = _stations_from_straight_prism(
            authored_solids[0], n_stations, n_outline, cad_units_per_meter)
    else:
        stations, tag = _stations_from_section_planes(
            solid, n_stations, n_outline, cad_units_per_meter)
    return align_station_origins(stations), tag


#: Body lane spacing as a fraction of the fin tip radius.  The tip
#: radius is the length scale the surface current crowds on, and a
#: straight-beak-fin convergence run at 150 kHz (2026-09-20, mid-station
#: of the 60 mm fixture against the independent 2-D SIBC reference)
#: measured: h_body = r_tip gave R/m +1.3 %, beak loss fraction -1.1 %;
#: h_body = r_tip/2 gave +0.3 % and -0.2 %; h_body = r_tip/3 gave +0.06 %
#: and -0.6 %.  Half the tip radius is where the delivery metrics enter
#: the 1 % band without paying the cubic solve cost of the next step.
LANE_SPACING_PER_TIP_RADIUS = 0.5

#: Dense outline samples per lane.  The same run showed 96 lanes on a
#: 1024-point outline and on a 2048-point outline disagreeing by 2 % on
#: the beak loss fraction purely from where the graded lane centres
#: quantised onto the sampled outline; >= 20 samples per lane removed it.
OUTLINE_SAMPLES_PER_LANE = 24

#: Cells on each tip arc.  Measured at a fixed 96-lane budget: 6-8 cells
#: sit closest to the 2-D reference, and raising the request past ~16
#: silently saturates against the grading's ``max_ratio`` clamp because
#: the lanes it needs come out of the body.
DEFAULT_TIP_LANES = 8


def auto_fin_resolution(step_path, *, cad_units_per_meter="auto",
                        route="auto", probe_stations=5, probe_outline=4096,
                        max_lanes=256, max_stations=41):
    """Measure one STEP conductor and choose its fin discretisation.

    A cheap CAD-only pass sections the solid at ``probe_stations``
    stations, analyses each section, and reports both the measurements
    and the ``n_lanes`` / ``n_outline`` / ``tip_lanes`` / ``n_stations``
    they imply.  Nothing is solved here.  Every returned number carries
    the measurement it came from so a run can be audited without
    re-deriving the rule.
    """
    stations, tag = station_outlines_from_step(
        step_path, n_stations=probe_stations, n_outline=probe_outline,
        cad_units_per_meter=cad_units_per_meter, route=route)
    analyses = [analyze_section(_as_closed_ccw(st.uv)) for st in stations]
    perimeter = max(float(a.perimeter) for a in analyses)
    areas = np.array([float(st.area) for st in stations])
    tip_radii = [float(a.primary.tip_radius)
                 for a in analyses if a.primary is not None]
    centroids = np.array([st.centroid for st in stations])
    sweep_length = float(np.sum(np.linalg.norm(
        centroids[1:] - centroids[:-1], axis=1)))

    if tip_radii:
        driver = "fin tip radius"
        lane_spacing = LANE_SPACING_PER_TIP_RADIUS * min(tip_radii)
    else:
        # No fin anywhere: there is no tip length scale to resolve, so
        # fall back on the section size.  Reported, never silent.
        driver = "section size (no fin detected)"
        lane_spacing = float(np.sqrt(np.mean(areas))) / 8.0
    n_lanes = int(np.clip(8 * math.ceil(perimeter / lane_spacing / 8),
                          32, max_lanes))
    n_outline = int(np.clip(
        2 ** math.ceil(math.log2(OUTLINE_SAMPLES_PER_LANE * n_lanes)),
        1024, 16384))
    section_size = float(np.sqrt(np.mean(areas)))
    n_stations = int(np.clip(
        round(sweep_length / (2.0 * section_size)) + 1, 5, max_stations))

    return {
        "n_lanes": n_lanes, "n_outline": n_outline,
        "tip_lanes": DEFAULT_TIP_LANES, "n_stations": n_stations,
        "lane_grading": "auto" if tip_radii else "uniform",
        "measured": {
            "cad_source": tag,
            "perimeter_m": perimeter,
            "section_area_m2_mean": float(np.mean(areas)),
            "section_area_relative_spread": float(
                (areas.max() - areas.min()) / np.mean(areas)),
            "tip_radius_m_min": min(tip_radii) if tip_radii else None,
            "fin_probe_stations": len(tip_radii),
            "probe_stations": len(stations),
            "sweep_length_m": sweep_length,
        },
        "rule": {
            "lane_spacing_driver": driver,
            "lane_spacing_m": lane_spacing,
            "outline_samples_per_lane": OUTLINE_SAMPLES_PER_LANE,
            "station_spacing_per_section_size": 2.0,
            "clamped_lanes": n_lanes == max_lanes,
            "clamped_stations": n_stations == max_stations,
        },
    }


def fin_graph_from_step(step_path, *, n_lanes=None, n_stations=None,
                        n_outline=None, lane_grading=None, tip_lanes=None,
                        cad_units_per_meter="auto", route="auto"):
    """One call: STEP -> origin-aligned stations -> fin-graded PEEC graph.

    Any discretisation argument left as ``None`` is measured from the
    geometry by :func:`auto_fin_resolution`; the chosen values and their
    basis land in ``meta["auto_resolution"]``.
    """
    requested = {"n_lanes": n_lanes, "n_stations": n_stations,
                 "n_outline": n_outline, "lane_grading": lane_grading,
                 "tip_lanes": tip_lanes}
    auto = None
    if any(v is None for v in requested.values()):
        auto = auto_fin_resolution(
            step_path, cad_units_per_meter=cad_units_per_meter, route=route)
        for key, value in requested.items():
            if value is None:
                requested[key] = auto[key]

    stations, tag = station_outlines_from_step(
        step_path, n_stations=requested["n_stations"],
        n_outline=requested["n_outline"],
        cad_units_per_meter=cad_units_per_meter, route=route)
    result = build_fin_graph(stations, n_lanes=requested["n_lanes"],
                             lane_grading=requested["lane_grading"],
                             tip_lanes=requested["tip_lanes"])
    result.meta["cad_source"] = tag
    result.meta["step_path"] = str(step_path)
    result.meta["auto_resolution"] = auto
    result.meta["resolution"] = dict(requested)
    return result
