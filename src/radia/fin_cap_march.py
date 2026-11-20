"""Terminal-to-terminal section chain of one STEP conductor.

The fin-surface route needs ordered cross-sections along the current
path.  The arc-spine sectioner in :mod:`radia.coil_from_cad` assumes the
path is an arc about the z axis, so straight leads, elbows and narrow
terminal gaps break it.  This module recovers the path from the solid
itself:

1. **Terminals.**  A planar face is a terminal cap when exactly one side
   of it is material (probed a tiny distance off the face, so a narrow
   gap to the other terminal is never crossed) and a section march
   started from it travels at least twice the face's own extent into
   the conductor.  The second test rejects flat side faces of bars and
   fins: marching "into" those immediately runs out of material.
   Exactly two caps are required; anything else raises.

2. **March.**  From cap A the march steps along the current tangent.
   At each predicted point the section plane is cut through the solid
   and the *component containing the predicted point* is taken; a
   section whose containing component is missing or ambiguous is never
   replaced by the nearest or the opposite-side component.  The
   tangent is corrected by the minimum-area section through the point,
   which is the perpendicular cut on a tube and stays well defined
   through bends and elbows.  Steps halve on a large turn, area jump or
   centroid jump and the march fails loudly below ``sqrt(A) / 64``.
   It ends when the predicted point crosses the plane of cap B inside
   cap B's footprint.

3. **Stations.**  The dense centreline is resampled at equal arc length.
   Each interior station is sectioned with the dense march's tangent,
   again keeping only the component that contains the centreline point.
   The end stations are the cap faces themselves.

Units: everything here runs in CAD units; the caller converts.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

#: Inward probe distance off a candidate cap face, relative to sqrt(area).
#: Small enough that a 2 mm terminal gap next to a 5 mm cap is never
#: crossed (the 0.3 sqrt(A) probe of the arc-spine detector did cross it).
CAP_PROBE_FRACTION = 1e-3

#: A cap must lead into at least this many face extents of conductor.
CAP_MIN_WALK_EXTENTS = 2.0

#: Predictor step as a fraction of sqrt(section area).
STEP_FRACTION = 0.25

#: Largest accepted tangent turn per march step [deg].
MAX_TURN_DEG = 45.0

#: Accepted section area band relative to the starting cap.
AREA_BAND = (0.4, 2.5)

#: Station-count rule: at most this much tangent turning per station [deg].
STATION_TURN_DEG = 10.0


@dataclass(frozen=True)
class TerminalCap:
    face: object            # TopoDS_Face
    centroid: np.ndarray
    inward: np.ndarray      # unit normal pointing into the conductor
    area: float
    extent: float           # largest bounding-box side of the face


@dataclass(frozen=True)
class Centerline:
    points: np.ndarray      # (n, 3) section centroids, cap A ... cap B
    tangents: np.ndarray    # (n, 3) unit section normals along the path
    areas: np.ndarray       # (n,)
    cap_a: TerminalCap
    cap_b: TerminalCap
    sections: int           # section evaluations spent

    @property
    def arclength(self) -> np.ndarray:
        d = np.linalg.norm(np.diff(self.points, axis=0), axis=1)
        return np.concatenate([[0.0], np.cumsum(d)])

    @property
    def total_turn_deg(self) -> float:
        c = np.clip(np.einsum("ij,ij->i", self.tangents[:-1], self.tangents[1:]), -1, 1)
        return float(np.degrees(np.sum(np.arccos(c))))


# ----------------------------------------------------------------------
# OCC primitives
# ----------------------------------------------------------------------

def _gp_pnt(p):
    from OCP.gp import gp_Pnt
    return gp_Pnt(float(p[0]), float(p[1]), float(p[2]))


def _inside(solid_w, p) -> bool:
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.TopAbs import TopAbs_IN
    cls = BRepClass3d_SolidClassifier(solid_w)
    cls.Perform(_gp_pnt(p), 1e-9)
    return cls.State() == TopAbs_IN


def _face_props(fw):
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    props = GProp_GProps()
    BRepGProp.SurfaceProperties_s(fw, props)
    com = props.CentreOfMass()
    return np.array([com.X(), com.Y(), com.Z()]), float(props.Mass())


def _plane_normal(fw):
    from OCP.BRep import BRep_Tool
    from OCP.BRepTools import BRepTools
    from OCP.GeomLProp import GeomLProp_SLProps
    from OCP.TopAbs import TopAbs_REVERSED
    surf = BRep_Tool.Surface_s(fw)
    umin, umax, vmin, vmax = BRepTools.UVBounds_s(fw)
    props = GeomLProp_SLProps(surf, 0.5 * (umin + umax), 0.5 * (vmin + vmax), 1, 1e-9)
    n = props.Normal()
    nn = np.array([n.X(), n.Y(), n.Z()])
    nn /= np.linalg.norm(nn)
    return -nn if fw.Orientation() == TopAbs_REVERSED else nn


def _point_in_planar_face(fw, p) -> bool:
    from OCP.BRep import BRep_Tool
    from OCP.BRepClass import BRepClass_FaceClassifier
    from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
    from OCP.TopAbs import TopAbs_IN, TopAbs_ON
    from OCP.gp import gp_Pnt2d
    proj = GeomAPI_ProjectPointOnSurf(_gp_pnt(p), BRep_Tool.Surface_s(fw))
    if proj.NbPoints() < 1:
        return False
    u, v = proj.LowerDistanceParameters()
    state = BRepClass_FaceClassifier(fw, gp_Pnt2d(u, v), 1e-7).State()
    return state in (TopAbs_IN, TopAbs_ON)


def _face_extent(fw) -> float:
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    box = Bnd_Box()
    BRepBndLib.Add_s(fw, box)
    lo, hi = box.CornerMin(), box.CornerMax()
    return float(max(hi.X() - lo.X(), hi.Y() - lo.Y(), hi.Z() - lo.Z()))


_DIAGONAL = {}


def _solid_diagonal(solid_w) -> float:
    key = solid_w.__hash__()
    if key not in _DIAGONAL:
        from OCP.Bnd import Bnd_Box
        from OCP.BRepBndLib import BRepBndLib
        box = Bnd_Box()
        BRepBndLib.Add_s(solid_w, box)
        lo, hi = box.CornerMin(), box.CornerMax()
        _DIAGONAL[key] = float(math.dist((lo.X(), lo.Y(), lo.Z()), (hi.X(), hi.Y(), hi.Z())))
    return _DIAGONAL[key]


def section_components(solid_w, point, normal):
    """Planar section of ``solid_w`` split into connected faces.

    A plane through a CAD edge or vertex can give OCC a degenerate (open)
    section; such a cut is retried with a tilt of 1e-5 and 1e-4 rad and
    raises :class:`MarchError` only when every attempt is degenerate.
    """
    n0 = np.asarray(normal, float)
    u, v = _frame(n0 / np.linalg.norm(n0))
    for tilt in (0.0, 1e-5, -1e-5, 1e-4, -1e-4):
        for axis in ((u,) if tilt == 0.0 else (u, v)):
            n = n0 + tilt * axis
            comps = _section_once(solid_w, point, n / np.linalg.norm(n))
            if comps is not None:
                return comps
    raise MarchError(
        f"the planar section through {np.round(point, 6).tolist()} stays degenerate "
        f"(open boundary wire) under small tilts")


def _section_once(solid_w, point, normal):
    """Planar section of ``solid_w`` split into connected faces.

    Returns a list of ``(face, centroid, area)``; inner wires (holes) are
    not supported by the fin route and make the face construction fail,
    which drops that component and is caught by the containment test.
    """
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    from OCP.ShapeAnalysis import ShapeAnalysis_FreeBounds
    from OCP.TopAbs import TopAbs_EDGE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    from OCP.gp import gp_Dir, gp_Pln
    try:
        from OCP.TopTools import TopTools_HSequenceOfShape
    except ImportError:  # OCP builds that expose NCollection sequences here
        from OCP.collections import HSequence_TopoDS_Shape as TopTools_HSequenceOfShape

    plane = gp_Pln(_gp_pnt(point), gp_Dir(*map(float, normal)))
    sec = BRepAlgoAPI_Section(solid_w, plane, False)
    sec.ComputePCurveOn1(True)
    sec.Approximation(False)
    sec.Build()
    if not sec.IsDone():
        return []
    edges = TopTools_HSequenceOfShape()
    exp = TopExp_Explorer(sec.Shape(), TopAbs_EDGE)
    while exp.More():
        edges.Append(exp.Current())
        exp.Next()
    if edges.Length() == 0:
        return []
    from OCP.BRep import BRep_Tool

    # Section edges of B-spline faces leave gaps far above 1e-6 CAD units;
    # join them at a tolerance tied to the solid size and keep closed wires.
    wires = TopTools_HSequenceOfShape()
    ShapeAnalysis_FreeBounds.ConnectEdgesToWires_s(
        edges, 1e-7 * _solid_diagonal(solid_w), False, wires)
    out = []
    for i in range(1, wires.Length() + 1):
        wire = TopoDS.Wire(wires.Value(i))
        if not BRep_Tool.IsClosed_s(wire):
            return None
        maker = BRepBuilderAPI_MakeFace(plane, wire, True)
        if not maker.IsDone():
            continue
        face = maker.Face()
        centroid, area = _face_props(face)
        if area > 0.0:
            out.append((face, centroid, area))
    return out


def _box_corners(fw):
    from OCP.Bnd import Bnd_Box
    from OCP.BRepBndLib import BRepBndLib
    box = Bnd_Box()
    BRepBndLib.Add_s(fw, box)
    lo, hi = box.CornerMin(), box.CornerMax()
    xs, ys, zs = (lo.X(), hi.X()), (lo.Y(), hi.Y()), (lo.Z(), hi.Z())
    return np.array([(x, y, z) for x in xs for y in ys for z in zs])


def _touches_cap(face, cap, tol) -> bool:
    """True when the section ``face`` touches the terminal ``cap`` face.

    Two bounding-box tests rule most calls out before the exact
    face-to-face distance: the section box lying wholly on one side of
    the cap plane, or not overlapping the cap's box.
    """
    corners = _box_corners(face)
    side = (corners - cap.centroid) @ cap.inward
    if float(side.min()) > tol or float(side.max()) < -tol:
        return False                       # wholly on one side of the cap plane
    cap_corners = _box_corners(cap.face)
    if (np.any(corners.min(axis=0) > cap_corners.max(axis=0) + tol)
            or np.any(corners.max(axis=0) < cap_corners.min(axis=0) - tol)):
        return False                       # bounding boxes do not overlap
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    dist = BRepExtrema_DistShapeShape(face, cap.face)
    return dist.IsDone() and dist.Value() <= tol


def _containing(components, point, caps=()):
    """The one component containing ``point``; None if missing or ambiguous.

    A component that touches a terminal cap is a cut through the terminal
    itself, not an interior cross-section, and never qualifies.
    """
    hit = [c for c in components if _point_in_planar_face(c[0], point)]
    if len(hit) != 1:
        return None
    comp = hit[0]
    for cap in caps:
        if (np.linalg.norm(comp[1] - cap.centroid) < 2.0 * cap.extent
                and _touches_cap(comp[0], cap, 1e-6 * cap.extent)):
            return None
    return comp


def _frame(t):
    a = np.array([1.0, 0.0, 0.0]) if abs(t[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(t, a)
    u /= np.linalg.norm(u)
    return u, np.cross(t, u)


def _min_area_section(solid_w, point, t, scale, caps=(), max_tilt_deg=60.0):
    """Normal near ``t`` minimising the area of the section containing ``point``."""
    from scipy.optimize import minimize

    u, v = _frame(t)
    lim = math.tan(math.radians(max_tilt_deg))
    cache = {}

    def evaluate(x):
        key = (round(float(x[0]), 7), round(float(x[1]), 7))
        if key not in cache:
            n = t + x[0] * u + x[1] * v
            n = n / np.linalg.norm(n)
            cache[key] = (_containing(section_components(solid_w, point, n), point, caps), n)
        return cache[key]

    def objective(x):
        if math.hypot(x[0], x[1]) > lim:
            return 1e300
        comp, _n = evaluate(x)
        return comp[2] if comp is not None else 1e300

    res = minimize(objective, np.zeros(2), method="Nelder-Mead",
                   options={"xatol": 1e-3, "fatol": 1e-6 * scale * scale,
                            "initial_simplex": [[0, 0], [0.15, 0], [0, 0.15]],
                            "maxfev": 60})
    comp, n = evaluate(res.x)
    return comp, n, len(cache)


# ----------------------------------------------------------------------
# march
# ----------------------------------------------------------------------

class MarchError(ValueError):
    pass


def _march(solid_w, start, *, cap_end=None, max_length=None, quick=False):
    caps = (start,) if cap_end is None else (start, cap_end)
    """March sections from ``start`` (a TerminalCap) until cap_end or no material."""
    scale = math.sqrt(start.area)
    h0 = STEP_FRACTION * scale
    pts, tans, areas = [start.centroid.copy()], [start.inward.copy()], [start.area]
    c, t = start.centroid.copy(), start.inward.copy()
    length, spent = 0.0, 0
    limit = max_length if max_length is not None else math.inf
    while length < limit:
        h = h0
        while True:
            p = c + h * t
            if cap_end is not None:
                cb = cap_end
                q = p - float(np.dot(p - cb.centroid, cb.inward)) * cb.inward
                if (np.dot(c - cb.centroid, cb.inward) > 0.0
                        and np.dot(p - cb.centroid, cb.inward) <= 0.0
                        and _point_in_planar_face(cb.face, q)):
                    pts.append(cb.centroid.copy())
                    tans.append(-cb.inward)
                    areas.append(cb.area)
                    return np.array(pts), np.array(tans), np.array(areas), "cap", spent
            comps = section_components(solid_w, p, t)
            spent += 1
            if not comps:
                if not _inside(solid_w, p):
                    return np.array(pts), np.array(tans), np.array(areas), "no material", spent
                u, _v = _frame(t)
                comps = section_components(solid_w, p, (t + 1e-4 * u) / np.linalg.norm(t + 1e-4 * u))
                spent += 1
                if not comps:
                    raise MarchError(
                        f"the planar section through {np.round(p, 6).tolist()} is empty "
                        f"although the point lies inside the solid (OCC section failure)")
            comp = _containing(comps, p, caps)
            n = None
            if comp is not None:
                chord = comp[1] - c
                n = chord / max(np.linalg.norm(chord), 1e-300)
                turn = math.degrees(math.acos(float(np.clip(n @ t, -1.0, 1.0))))
                smooth = turn <= 5.0 and abs(comp[2] / areas[-1] - 1.0) <= 0.05
                if not (quick or smooth):
                    comp = None
            if comp is None and not quick:
                comp, n, k = _min_area_section(solid_w, p, t, scale, caps)
                spent += k
            if comp is None:
                h *= 0.5
                if h < scale / 64.0:
                    raise MarchError(
                        f"no unique section contains the predicted point "
                        f"{np.round(p, 6).tolist()} (step reduced below sqrt(A)/64)")
                continue
            if n @ t < 0.0:
                n = -n
            turn = math.degrees(math.acos(float(np.clip(n @ t, -1.0, 1.0))))
            ratio = comp[2] / start.area
            jump = float(np.linalg.norm(comp[1] - c))
            if (turn > MAX_TURN_DEG or not (AREA_BAND[0] <= ratio <= AREA_BAND[1])
                    or jump > 3.0 * h + 0.5 * scale):
                h *= 0.5
                if h < scale / 64.0:
                    raise MarchError(
                        f"section chain breaks near {np.round(c, 6).tolist()}: "
                        f"turn {turn:.1f} deg, area ratio {ratio:.2f}, "
                        f"centroid jump {jump:.4g}")
                continue
            break
        length += jump
        c, t = comp[1], n
        pts.append(c.copy())
        tans.append(t.copy())
        areas.append(comp[2])
    return np.array(pts), np.array(tans), np.array(areas), "length limit", spent


def find_terminal_caps(solid_w):
    """The two planar terminal faces of one conductor solid.

    Raises :class:`MarchError` unless exactly two caps qualify; the
    message lists every planar candidate and how far it led.
    """
    from OCP.BRep import BRep_Tool
    from OCP.GeomAbs import GeomAbs_Plane
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS

    del BRep_Tool
    caps, report = [], []
    exp = TopExp_Explorer(solid_w, TopAbs_FACE)
    while exp.More():
        fw = TopoDS.Face(exp.Current())
        exp.Next()
        if BRepAdaptor_Surface(fw).GetType() != GeomAbs_Plane:
            continue
        centroid, area = _face_props(fw)
        if area <= 0.0:
            continue
        normal = _plane_normal(fw)
        eps = CAP_PROBE_FRACTION * math.sqrt(area)
        plus, minus = _inside(solid_w, centroid + eps * normal), _inside(solid_w, centroid - eps * normal)
        if plus == minus:
            continue
        cand = TerminalCap(fw, centroid, normal if plus else -normal, area, _face_extent(fw))
        try:
            pts, _t, _a, _why, _k = _march(
                solid_w, cand, max_length=(CAP_MIN_WALK_EXTENTS + 0.5) * cand.extent, quick=True)
            walk = float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
        except MarchError:
            walk = 0.0
        ok = walk >= CAP_MIN_WALK_EXTENTS * cand.extent
        report.append(f"area {area:.6g} at {np.round(centroid, 6).tolist()}: "
                      f"walk {walk:.4g} vs extent {cand.extent:.4g} -> {'cap' if ok else 'side'}")
        if ok:
            caps.append(cand)
    if len(caps) != 2:
        raise MarchError(
            f"expected exactly two terminal caps, found {len(caps)}; planar faces:\n  "
            + "\n  ".join(report or ["none"])
            + "\nThe fin route needs one conductor with two planar terminal faces.")
    return tuple(caps)


def march_centerline(solid_w) -> Centerline:
    """Cap A -> cap B section centroids of one conductor solid."""
    cap_a, cap_b = find_terminal_caps(solid_w)
    pts, tans, areas, why, spent = _march(solid_w, cap_a, cap_end=cap_b)
    if why != "cap":
        raise MarchError(
            f"section march from the terminal at {np.round(cap_a.centroid, 6).tolist()} "
            f"ended by '{why}' {np.linalg.norm(pts[-1] - cap_b.centroid):.4g} CAD units "
            f"short of the other terminal")
    return Centerline(pts, tans, areas, cap_a, cap_b, spent)


# ----------------------------------------------------------------------
# stations
# ----------------------------------------------------------------------

def sample_face_outline(fw, origin, u_hat, v_hat, n):
    """(n, 2) CCW equal-arc-length samples of a face's outer wire in (u, v)."""
    from OCP.BRepAdaptor import BRepAdaptor_CompCurve
    from OCP.BRepTools import BRepTools
    from OCP.GCPnts import GCPnts_AbscissaPoint, GCPnts_UniformAbscissa

    curve = BRepAdaptor_CompCurve(BRepTools.OuterWire_s(fw))
    total = GCPnts_AbscissaPoint.Length_s(curve)
    if total <= 0.0:
        raise MarchError("section face has a zero-length outer wire")
    ua = GCPnts_UniformAbscissa(curve, int(n) + 1, 1e-9)
    if not ua.IsDone() or ua.NbPoints() != n + 1:
        raise MarchError("arc-length sampling of a section outline failed")
    xyz = np.empty((n, 3))
    for k in range(n):
        q = curve.Value(ua.Parameter(k + 1))
        xyz[k] = (q.X(), q.Y(), q.Z())
    d = xyz - origin
    uv = np.column_stack([d @ u_hat, d @ v_hat])
    signed = float(np.sum(uv[:, 0] * np.roll(uv[:, 1], -1) - uv[:, 1] * np.roll(uv[:, 0], -1)))
    if signed < 0.0:
        uv = uv[::-1].copy()
    start = int(np.argmax(uv[:, 0] - 1e-9 * np.abs(uv[:, 1])))
    return np.roll(uv, -start, axis=0), float(total)


def section_frames(points, tangents):
    """(u, v) per station, orthogonal to each station tangent.

    A coplanar path (all leads and rings in one plane) uses u = the path
    plane normal, so the frame never twists; otherwise u is parallel
    transported.  v = t x u in both cases.
    """
    pts = np.asarray(points, float)
    t = np.asarray(tangents, float)
    centered = pts - pts.mean(axis=0)
    _s, sv, vt = np.linalg.svd(centered, full_matrices=False)
    if sv[-1] <= 1e-6 * max(sv[0], 1e-300):
        seed = vt[-1]
        if seed[np.argmax(np.abs(seed))] < 0:
            seed = -seed
        u = seed[None, :] - (t @ seed)[:, None] * t
    else:
        u = np.empty_like(t)
        seed = np.array([0.0, 0.0, 1.0]) if abs(t[0, 2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        u[0] = seed - (seed @ t[0]) * t[0]
        for i in range(1, len(t)):
            u[i] = u[i - 1] - (u[i - 1] @ t[i]) * t[i]
    u /= np.linalg.norm(u, axis=1)[:, None]
    return u, np.cross(t, u)


STATION_AREA_TOL = 0.25


def _station_section(solid_w, p, t, band, caps, index):
    """The station section through ``p``, consistent with the dense march.

    Its area must lie in ``band`` (the neighbouring march sections widened
    by STATION_AREA_TOL) and its centroid within sqrt(area) of ``p``.  A
    plane through CAD edges can make OCC return one merged face (2026-11-19:
    a 139 mm^2 face across a ring of 24.6 mm^2 sections); such a cut is
    retried with tilts up to 1e-3 rad, then raises.
    """
    u, v = _frame(t)
    seen = []
    for tilt in (0.0, 1e-4, -1e-4, 1e-3, -1e-3):
        for axis in ((u,) if tilt == 0.0 else (u, v)):
            n = t + tilt * axis
            n = n / np.linalg.norm(n)
            comp = _containing(section_components(solid_w, p, n), p, caps)
            if comp is None:
                seen.append(None)
                continue
            seen.append(round(float(comp[2]), 9))
            if (band[0] <= comp[2] <= band[1]
                    and np.linalg.norm(comp[1] - p) <= math.sqrt(band[1])):
                return comp, n
    raise MarchError(
        f"station {index}: no section through {np.round(p, 6).tolist()} matches the "
        f"march (area band {band[0]:.6g}..{band[1]:.6g}; got {seen})")


FEATURE_AREA_STEP = 0.03   # relative section-area change that marks a feature end (fin start/end)
FEATURE_PAIR_FRACTION = 1.0 / 16.0   # feature station offset, fraction of the grid spacing


def station_targets(line: Centerline, n_stations):
    """Arc-length station positions aligned to the section-area features.

    An equal-arc grid of ``n_stations`` lets a feature that starts or ends
    between stations (a fin along the sweep) gain or lose up to one spacing
    of length depending on the count, so refinement is not monotone.  Each
    feature end (where the dense-march area crosses ``FEATURE_AREA_STEP``
    relative to the median area) gets a station pair at +/- h/16 around it
    (h = grid spacing); grid points within h/3 of a pair are dropped.  The
    returned count therefore differs from ``n_stations``.
    """
    if n_stations < 3:
        raise ValueError("need n_stations >= 3")
    s = line.arclength
    a = np.asarray(line.areas, dtype=float)
    grid = np.linspace(0.0, s[-1], n_stations)
    h = s[-1] / (n_stations - 1)
    ref = float(np.median(a))
    if not np.isfinite(ref) or ref <= 0.0:
        raise MarchError("station targets: non-positive median section area")
    inside = np.abs(a / ref - 1.0) > FEATURE_AREA_STEP
    ends = [0.5 * (s[k] + s[k + 1]) for k in np.nonzero(inside[1:] != inside[:-1])[0]]
    d = FEATURE_PAIR_FRACTION * h
    # A feature end within h/3 of a cap is resolved by the cap station itself;
    # a pair there would be one-sided or nearly coincide with the cap.
    ends = [e for e in ends if h / 3 < e < s[-1] - h / 3]
    if not ends:
        return grid
    if len(ends) > 1 and np.min(np.diff(ends)) <= 2 * h / 3:
        raise MarchError(f"station targets: section-area features {np.round(ends, 6).tolist()} "
                         f"are closer than 2h/3 = {2 * h / 3:.6g}; raise n_stations")
    pairs = np.asarray([x for e in ends for x in (e - d, e + d)])
    keep = [g for g in grid if g in (0.0, s[-1]) or np.min(np.abs(pairs - g)) > h / 3]
    t = np.unique(np.concatenate([keep, pairs]))
    if np.min(np.diff(t)) < 2 * d * (1 - 1e-9):
        raise MarchError(f"station targets: spacing {np.min(np.diff(t)):.6g} below 2h/16")
    return t


def station_sections(solid_w, line: Centerline, n_stations):
    """Feature-aligned stations (see station_targets): list of (face, centroid, tangent, area)."""
    s = line.arclength
    targets = station_targets(line, n_stations)
    out = [(line.cap_a.face, line.cap_a.centroid, line.cap_a.inward, line.cap_a.area)]
    for i, st in enumerate(targets[1:-1], start=1):
        k = int(np.clip(np.searchsorted(s, st) - 1, 0, len(s) - 2))
        w = (st - s[k]) / max(s[k + 1] - s[k], 1e-300)
        p = (1 - w) * line.points[k] + w * line.points[k + 1]
        t = (1 - w) * line.tangents[k] + w * line.tangents[k + 1]
        t /= np.linalg.norm(t)
        band = (min(line.areas[k], line.areas[k + 1]) / (1 + STATION_AREA_TOL),
                max(line.areas[k], line.areas[k + 1]) * (1 + STATION_AREA_TOL))
        comp, n = _station_section(solid_w, p, t, band, (line.cap_a, line.cap_b), i)
        out.append((comp[0], comp[1], n, comp[2]))
    out.append((line.cap_b.face, line.cap_b.centroid, -line.cap_b.inward, line.cap_b.area))
    return out


def stations_for_turning(line: Centerline) -> int:
    """Station count that keeps the per-station tangent turn <= STATION_TURN_DEG."""
    return int(math.ceil(line.total_turn_deg / STATION_TURN_DEG)) + 1
