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
   frame: direct z-sectioning for straight prisms, otherwise the
   terminal-to-terminal section march of :mod:`radia.fin_cap_march`
   (leads, elbows and narrow terminal gaps included).  No PEEC solver is
   built here.
2. ``align_station_origins`` — the first station starts at a geometric
   origin (:func:`radia.fin_cap_march.outline_origin`); each later outline
   is expressed in its predecessor's frame transported onto its own
   tangent and cyclically re-registered to the predecessor, so lane k
   follows the same material line along the sweep whatever the local
   in-plane frames or wire seams are.
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
from pathlib import Path

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
    """Return stations re-registered so lane origins follow material lines.

    Station ``i`` is compared with station ``i-1`` in the predecessor's
    (u, v) parallel transported onto station ``i``'s tangent, so frames
    rotated independently per station do not shift the registration.  Each
    station keeps its own frame; only its sample order is rolled.
    """
    from radia.fin_cap_march import transport_vector

    if not stations:
        return ()
    out = [stations[0]]
    for st in stations[1:]:
        prev = out[-1]
        t0 = prev.tangent / np.linalg.norm(prev.tangent)
        t1 = st.tangent / np.linalg.norm(st.tangent)
        u_t = transport_vector(prev.u_hat, t0, t1)
        u_t = u_t - float(u_t @ t1) * t1
        u_t /= np.linalg.norm(u_t)
        v_t = np.cross(t1, u_t)
        q = st.uv[:, :1] * st.u_hat[None, :] + st.uv[:, 1:2] * st.v_hat[None, :]
        common = np.column_stack([q @ u_t, q @ v_t])
        _aligned, shift = align_outline_origin(prev.uv, common)
        out.append(StationOutline(st.index, st.centroid, st.u_hat, st.v_hat,
                                  st.tangent, np.roll(st.uv, -shift, axis=0),
                                  st.area, st.perimeter))
    return tuple(out)


# ----------------------------------------------------------------------
# graph construction
# ----------------------------------------------------------------------

def _lane_fraction_for_station(analysis, n_lanes, lane_grading, tip_lanes, lane_phase):
    if lane_grading == "auto" and analysis.fins:
        s = graded_lane_arclengths(analysis, n_lanes, tip_lanes=tip_lanes)
        return s / analysis.perimeter, True
    if lane_phase == "anchored":  # lane 0 on the geometric outline origin
        return np.arange(n_lanes) / n_lanes, False
    return (np.arange(n_lanes) + 0.5) / n_lanes, False


def build_fin_graph(stations, *, n_lanes=64, lane_grading="auto",
                    tip_lanes=8, mesh_stations="all", fin_kwargs=None,
                    lane_phase="centered"):
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
        lane_phase: ``"anchored"`` puts lane 0 on the outline origin (the
            cap-march route starts outlines at a geometric anchor);
            ``"centered"`` offsets uniform lanes by half a spacing.
    """
    from radia.peec_fin_topology import build_hybrid_surface_topology

    if lane_grading not in {"uniform", "auto"}:
        raise ValueError("lane_grading must be 'uniform' or 'auto'")
    if lane_phase not in {"centered", "anchored"}:
        raise ValueError("lane_phase must be 'centered' or 'anchored'")
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
            analysis, n_lanes, lane_grading, tip_lanes, lane_phase)
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


_CENTERLINE_CACHE = {}


def _cap_march_centerline(step_path, solid):
    """Dense terminal-to-terminal centreline, cached per STEP file content."""
    import hashlib

    from radia.fin_cap_march import march_centerline

    digest = hashlib.sha256(Path(step_path).read_bytes()).hexdigest()
    line = _CENTERLINE_CACHE.get(digest)
    if line is None:
        line = march_centerline(solid.wrapped)
        _CENTERLINE_CACHE.clear()
        _CENTERLINE_CACHE[digest] = line
    return line


def cap_march_station_outlines(outline_xyz, centroids, tangents, areas,
                               perimeters, cad_units_per_meter=1.0):
    """Station outlines from sampled section wires (pure numpy, any seam).

    Frames come from :func:`radia.fin_cap_march.section_frames`; every
    outline is put counter-clockwise in its frame and started at
    :func:`radia.fin_cap_march.outline_origin` along the transported ``u``.
    Registration (:func:`align_station_origins`) then carries the first
    station's origin along the sweep.
    """
    from radia.fin_cap_march import outline_in_frame, outline_origin, section_frames

    centroids = np.asarray(centroids, float)
    tangents = np.asarray(tangents, float)
    u, v = section_frames(centroids, tangents)
    stations = []
    for i, xyz in enumerate(outline_xyz):
        uv = outline_in_frame(xyz, centroids[i], u[i], v[i])
        uv = np.roll(uv, -outline_origin(uv, np.array([1.0, 0.0])), axis=0)
        stations.append(StationOutline(
            i, centroids[i] / cad_units_per_meter, u[i], v[i],
            tangents[i] / np.linalg.norm(tangents[i]),
            uv / cad_units_per_meter, float(areas[i]) / cad_units_per_meter**2,
            float(perimeters[i]) / cad_units_per_meter))
    return tuple(stations)


def _stations_from_cap_march(step_path, solid, n_stations, n_outline,
                             cad_units_per_meter):
    from radia.fin_cap_march import sample_face_outline, station_sections

    line = _cap_march_centerline(step_path, solid)
    sections = station_sections(solid.wrapped, line, n_stations)
    sampled = [sample_face_outline(sec[0], int(n_outline)) for sec in sections]
    stations = cap_march_station_outlines(
        [xyz for xyz, _ in sampled], [sec[1] for sec in sections],
        [sec[2] for sec in sections], [sec[3] for sec in sections],
        [per for _, per in sampled], cad_units_per_meter)
    return stations, "step_cap_march"


def station_outlines_from_step(step_path, *, n_stations=20, n_outline=2048,
                               cad_units_per_meter="auto", route="auto"):
    """Dense per-station outlines of one STEP conductor.

    ``route``: ``"auto"`` picks the straight-prism sectioner for a single
    z-extruded prism and the terminal-to-terminal section march
    (:mod:`radia.fin_cap_march`) otherwise; ``"straight_prism"`` or
    ``"cap_march"`` force one.  ``cad_units_per_meter="auto"``
    defers to :func:`resolve_cad_units_per_meter`, which answers only when
    the STEP settles it.  Returns ``(stations, source_tag)`` with stations
    already origin-aligned.
    """
    from radia._b3d_shim import import_step
    from radia.coil_from_cad import _check_solid_extent_plausible

    if route not in {"auto", "straight_prism", "cap_march"}:
        raise ValueError("route must be auto, straight_prism or cap_march")
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
        route = "straight_prism" if _is_straight_prism(solid) else "cap_march"
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
        stations, tag = _stations_from_cap_march(
            step_path, solid, n_stations, n_outline, cad_units_per_meter)
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

#: How much of the section's total change one station may absorb.  A
#: station sweep at 150 kHz (2026-09-20) compared a prismatic beak fin
#: with a tapered one whose fin runs 2.40 -> 0.60 mm over 60 mm.  The
#: prismatic terminal resistance moved 0.17 % between 5 and 33 stations,
#: so a constant section needs nothing beyond the aspect rule.  The
#: tapered one kept climbing -- 330.3, 332.4, 334.0, 336.4, 336.5 uOhm at
#: 9, 13, 17, 25, 33 stations -- and only settled to 0.03 % per step
#: between 25 and 33.  Its fin length spans 75 %, and 0.75 / 0.03 lands
#: at 26 stations, which is where that settling begins.
STATION_RELATIVE_CHANGE = 0.03

#: Dense partial-element matrices cost 16 * N**2 bytes and factor in
#: O(N**3), so the measured discretisation is checked against a branch
#: budget rather than being quietly coarsened.  12000 branches is ~2.3 GB.
DEFAULT_MAX_BRANCHES = 12000


def auto_fin_resolution(step_path, *, cad_units_per_meter="auto",
                        route="auto", probe_stations=5, probe_outline=4096,
                        max_lanes=256, max_branches=DEFAULT_MAX_BRANCHES,
                        n_lanes=None, n_stations=None):
    """Measure one STEP conductor and choose its fin discretisation.

    A cheap CAD-only pass sections the solid at ``probe_stations``
    stations, analyses each section, and reports both the measurements
    and the ``n_lanes`` / ``n_outline`` / ``tip_lanes`` / ``n_stations``
    they imply.  Nothing is solved here.  Every returned number carries
    the measurement it came from so a run can be audited without
    re-deriving the rule.

    ``n_lanes`` / ``n_stations`` given by the caller replace the measured
    values before the outline size and the branch budget are derived, so
    pinning one of them is how a caller resolves a budget refusal.

    Raises when the measured discretisation would exceed ``max_branches``.
    A section that both crowds its current on a small tip radius and
    changes along the sweep wants fine lanes AND many stations, and the
    dense kernel cannot afford both; coarsening one of them quietly would
    hide that, so the caller is told the cost and pins a value itself.
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
    measured_lanes = int(np.clip(8 * math.ceil(perimeter / lane_spacing / 8),
                                 32, max_lanes))
    n_lanes = measured_lanes if n_lanes is None else int(n_lanes)
    n_outline = int(np.clip(
        2 ** math.ceil(math.log2(OUTLINE_SAMPLES_PER_LANE * n_lanes)),
        1024, 16384))
    # Stations follow how much the section CHANGES along the sweep, not
    # how long the sweep is: a prismatic fin is insensitive to the count.
    fin_lengths = [float(a.primary.length)
                   for a in analyses if a.primary is not None]
    if len(fin_lengths) < len(analyses):
        fin_span = 1.0                      # the fin vanishes somewhere
    elif fin_lengths:
        fin_span = (max(fin_lengths) - min(fin_lengths)) / max(fin_lengths)
    else:
        fin_span = 0.0
    area_span = float((areas.max() - areas.min()) / np.mean(areas))
    section_change = max(fin_span, area_span)

    # A curved path also needs stations for its turning, which five probe
    # sections cannot see; the cap march carries the dense centreline.
    total_turn_deg = 0.0
    n_from_turn = 0
    if tag == "step_cap_march":
        from radia._b3d_shim import import_step
        from radia.fin_cap_march import STATION_TURN_DEG, stations_for_turning

        line = _cap_march_centerline(
            step_path, import_step(str(step_path)).solids()[0])
        units = math.sqrt(line.cap_a.area / stations[0].area)
        sweep_length = float(line.arclength[-1]) / units
        total_turn_deg = line.total_turn_deg
        n_from_turn = stations_for_turning(line)

    section_size = float(np.sqrt(np.mean(areas)))
    n_from_aspect = round(sweep_length / (2.0 * section_size)) + 1
    n_from_change = math.ceil(section_change / STATION_RELATIVE_CHANGE) + 1
    measured_stations = max(5, int(n_from_aspect), int(n_from_change), int(n_from_turn))
    n_stations = measured_stations if n_stations is None else int(n_stations)

    branches = (2 * n_stations - 3) * n_lanes
    if branches > max_branches:
        raise ValueError(
            f"the measured discretisation needs {branches} branches "
            f"({n_lanes} lanes x {n_stations} stations), over the "
            f"{max_branches} the dense partial-element kernel is budgeted "
            f"for (~{16 * branches ** 2 / 2**30:.1f} GB).  Lanes come from a "
            f"{min(tip_radii) * 1e3:.4g} mm tip radius on a "
            f"{perimeter * 1e3:.4g} mm perimeter; stations from a "
            f"{section_change:.0%} section change along the sweep.  Pin "
            f"n_lanes or n_stations to choose which one to give up, or "
            f"raise max_branches."
            if tip_radii else
            f"the measured discretisation needs {branches} branches "
            f"({n_lanes} lanes x {n_stations} stations), over the "
            f"{max_branches} budget.  Pin n_lanes or n_stations, or raise "
            f"max_branches.")

    return {
        "n_lanes": n_lanes, "n_outline": n_outline,
        "tip_lanes": DEFAULT_TIP_LANES, "n_stations": n_stations,
        "lane_grading": "auto" if tip_radii else "uniform",
        "measured": {
            "cad_source": tag,
            "perimeter_m": perimeter,
            "section_area_m2_mean": float(np.mean(areas)),
            "section_area_relative_spread": area_span,
            "fin_length_relative_span": fin_span,
            "tip_radius_m_min": min(tip_radii) if tip_radii else None,
            "fin_probe_stations": len(tip_radii),
            "probe_stations": len(stations),
            "sweep_length_m": sweep_length,
            "path_turn_deg": total_turn_deg,
        },
        "rule": {
            "lane_spacing_driver": driver,
            "lane_spacing_m": lane_spacing,
            "outline_samples_per_lane": OUTLINE_SAMPLES_PER_LANE,
            "station_spacing_per_section_size": 2.0,
            "station_relative_change": STATION_RELATIVE_CHANGE,
            "stations_from_aspect": int(n_from_aspect),
            "stations_from_section_change": int(n_from_change),
            "station_turn_deg": STATION_TURN_DEG if n_from_turn else None,
            "stations_from_turning": int(n_from_turn),
            "estimated_branches": int(branches),
            "clamped_lanes": measured_lanes == max_lanes,
            "measured_lanes": measured_lanes,
            "measured_stations": measured_stations,
            "pinned": sorted(k for k, v in (("n_lanes", n_lanes != measured_lanes),
                                            ("n_stations", n_stations != measured_stations)) if v),
        },
    }


def _require_branch_budget(n_lanes, n_stations, max_branches):
    branches = (2 * int(n_stations) - 3) * int(n_lanes)
    if branches > max_branches:
        raise ValueError(
            f"the requested discretisation needs {branches} branches "
            f"({n_lanes} lanes x {n_stations} stations), over the "
            f"{max_branches} the dense partial-element kernel is budgeted "
            f"for (~{16 * branches ** 2 / 2**30:.1f} GB).  Lower n_lanes or "
            f"n_stations, or raise max_branches.")
    return branches


def fin_graph_from_step(step_path, *, n_lanes=None, n_stations=None,
                        n_outline=None, lane_grading=None, tip_lanes=None,
                        cad_units_per_meter="auto", route="auto",
                        max_branches=DEFAULT_MAX_BRANCHES):
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
            step_path, cad_units_per_meter=cad_units_per_meter, route=route,
            n_lanes=n_lanes, n_stations=n_stations, max_branches=max_branches)
        for key, value in requested.items():
            if value is None:
                requested[key] = auto[key]
    # The same budget holds whether the values were measured or all pinned.
    _require_branch_budget(requested["n_lanes"], requested["n_stations"], max_branches)

    stations, tag = station_outlines_from_step(
        step_path, n_stations=requested["n_stations"],
        n_outline=requested["n_outline"],
        cad_units_per_meter=cad_units_per_meter, route=route)
    result = build_fin_graph(stations, n_lanes=requested["n_lanes"],
                             lane_grading=requested["lane_grading"],
                             tip_lanes=requested["tip_lanes"],
                             lane_phase="anchored" if tag == "step_cap_march" else "centered")
    result.meta["cad_source"] = tag
    result.meta["step_path"] = str(step_path)
    result.meta["auto_resolution"] = auto
    result.meta["resolution"] = dict(requested)
    return result
