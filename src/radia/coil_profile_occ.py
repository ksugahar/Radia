"""coil_profile_occ.py

build123d / OCC bridge for the filament-approximation PEEC profile
abstractions in coil_profile.py.

Exposes:
    profile_from_build123d_face(face, n_boundary=64, u_axis=None) ->
        PolygonProfile projected onto the face's plane. Sign convention:
        the face's normal is the +y_hat of the coil local frame (flow
        direction), and (u_axis, v_axis) span the cross-section plane.

    polygon_profile_from_points_3d(points, normal, origin, u_axis=None) ->
        Lower-level helper that projects ordered 3D points onto a plane
        and returns (PolygonProfile, u_axis, v_axis).

Usage (CAD-first coil design):

    from build123d import Rectangle, Circle, loft, Plane
    from radia.coil_profile_occ import profile_from_build123d_face
    from radia.coil_builder import CoilBuilder

    # Design the shape in build123d
    p0 = Rectangle(5e-3, 5e-3).face()
    p1 = Circle(2.5e-3).face()

    # Convert CAD faces to PEEC profiles
    prof0 = profile_from_build123d_face(p0)
    prof1 = profile_from_build123d_face(p1)

    # Assemble coil for PEEC filament analysis
    coil = (CoilBuilder(current=100)
            .set_start([0, 0, 0])
            .set_profile(prof0)
            .add_loft_straight(profile_end=prof1, length=40e-3))
    paths, _ = coil.to_filaments(nw=5, nh=8, frequency=7e3)
    # ... feed into peec_bundle ...
"""

import numpy as np

from radia.coil_profile import PolygonProfile, CircleProfile, RectProfile


def _pick_orthogonal(n):
    """Pick a unit vector orthogonal to n."""
    n = np.asarray(n, dtype=float)
    n = n / (np.linalg.norm(n) + 1e-30)
    # Use the smallest absolute component to seed the cross-product
    i = int(np.argmin(np.abs(n)))
    axis = np.zeros(3)
    axis[i] = 1.0
    u = np.cross(n, axis)
    u /= np.linalg.norm(u) + 1e-30
    return u


def polygon_profile_from_points_3d(points, normal, origin, u_axis=None):
    """Project ordered 3D boundary points onto a plane and return
    (PolygonProfile, u_axis, v_axis).

    Args:
        points: (n, 3) array-like of ordered (CCW looking down the +normal
            direction) 3D boundary points on a plane.
        normal: (3,) face normal. Defines which side of the plane is +y
            in the coil local frame.
        origin: (3,) origin point on the plane (projects to (u, v)=(0,0)).
        u_axis: optional (3,) axis in the plane; defines the +u direction.
            Defaults to a canonical axis orthogonal to normal.

    Returns:
        (PolygonProfile, u_axis, v_axis) where u_axis, v_axis are the
        chosen 3D axes spanning the plane (unit vectors, mutually
        orthogonal, v = normal x u).
    """
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3 or pts.shape[0] < 3:
        raise ValueError(
            f"points must be (n>=3, 3); got {pts.shape}")
    n = np.asarray(normal, dtype=float)
    n = n / (np.linalg.norm(n) + 1e-30)
    o = np.asarray(origin, dtype=float)

    if u_axis is None:
        u = _pick_orthogonal(n)
    else:
        u = np.asarray(u_axis, dtype=float)
        u = u - np.dot(u, n) * n
        u /= np.linalg.norm(u) + 1e-30
    v = np.cross(n, u)

    rel = pts - o[None, :]
    uv = np.column_stack([rel @ u, rel @ v])
    return PolygonProfile(uv), u, v


_CURVED_GEOM_TYPES = {
    'CIRCLE', 'ELLIPSE', 'HYPERBOLA', 'PARABOLA', 'BSPLINE', 'BEZIER',
    'OFFSET', 'OTHER',
}


def _n_samples_for_edge(edge, n_curved, n_straight):
    """Pick boundary sample count for an edge based on its geometry.

    Straight edges (LINE) need only the endpoints worth of parameter
    points; curved edges (CIRCLE, SPLINE, ...) need many to approximate
    the arc faithfully.
    """
    try:
        gt = edge.geom_type
        # build123d >= 0.9 returns GeomType enum; convert to name.
        name = getattr(gt, 'name', str(gt)).upper()
    except Exception:
        name = ''
    if name == 'LINE':
        return n_straight
    if name in _CURVED_GEOM_TYPES:
        return n_curved
    # Unknown -> err on the side of more samples
    return n_curved


def profile_from_build123d_face(face, n_per_edge=None, u_axis=None,
                                  deduplicate=True, n_curved=32,
                                  n_straight=2):
    """Convert a build123d Face (assumed planar) into a PolygonProfile.

    The outer wire is walked in order and each edge sampled at
    n_samples(edge) parameter values. Straight LINE edges are sampled
    with only ``n_straight`` points by default (the endpoints are
    enough) while curved edges (CIRCLE, ELLIPSE, BSPLINE, ...) get
    ``n_curved`` points to keep the polygonized boundary close to the
    true curve. Passing an explicit ``n_per_edge`` overrides both.

    Boundary points are projected onto the face plane using
    face.center() as origin and face.normal_at() as +y_hat normal.

    Args:
        face: build123d Face (or any object exposing .edges(),
            .normal_at(), .center()).
        n_per_edge: if given, used for every edge irrespective of type.
            Default None -> auto-pick per edge via geom_type.
        u_axis: optional 3D (plane-projected) vector defining +u in the
            profile frame. Defaults to canonical.
        deduplicate: drop consecutive coincident boundary points.
        n_curved: samples per curved edge (default 32).
        n_straight: samples per LINE edge (default 2; endpoints only).

    Returns:
        PolygonProfile in the face's local (u, v) plane -- OR a
        CircleProfile when the face is a single closed circle (detected
        via one CIRCLE edge whose length matches 2*pi*r). For CAD
        workflows that author round wire in build123d, this yields the
        exact area-preserving polar parametrization that PEEC filament
        sampling expects (same as authoring CircleProfile(r) directly).
    """
    try:
        outer = face.outer_wire()
    except AttributeError:
        outer = face.wires()[0]

    try:
        edges = outer.edges()
    except AttributeError:
        edges = face.edges()

    # --- Fast path: detect a pure-circle face and return CircleProfile.
    # Triggered when the outer wire has exactly one CIRCLE edge whose
    # length equals 2*pi*r (full circle within tol).
    edge_list = list(edges)
    if len(edge_list) == 1:
        e0 = edge_list[0]
        if _edge_geom_name(e0) == 'CIRCLE':
            try:
                r0 = float(e0.radius)
                if r0 > 0 and abs(float(e0.length) - 2 * np.pi * r0) < 1e-9 * r0:
                    return CircleProfile(r=r0)
            except AttributeError:
                pass

    # --- Fast path: detect a rectangle (4 LINE edges, right angles,
    # opposite sides equal) aligned with build123d's (X, Y) axes and
    # return RectProfile(w=X-extent, h=Y-extent).
    # Convention: for Rectangle(w, h) in Plane.XY, we want RectProfile(w, h)
    # with width along the FACE's local +X (= world +X) and height along
    # +Y (= world +Y). This matches the visual expectation of a CAD user
    # drawing the cross-section in standard XY.
    if len(edge_list) == 4 and all(_edge_geom_name(e) == 'LINE'
                                    for e in edge_list):
        verts = [_vec_np(e.position_at(0.0)) for e in edge_list]
        sides = np.array([verts[(i + 1) % 4] - verts[i] for i in range(4)])
        lens = np.linalg.norm(sides, axis=1)
        if lens.min() > 1e-12:
            tol = 1e-9 * lens.max()
            sides_equal = (abs(lens[0] - lens[2]) < tol
                           and abs(lens[1] - lens[3]) < tol)
            d01 = abs(np.dot(sides[0], sides[1]))
            d12 = abs(np.dot(sides[1], sides[2]))
            perp = (d01 < tol * lens[0] * lens[1]
                    and d12 < tol * lens[1] * lens[2])
            if sides_equal and perp:
                # Use X-extent and Y-extent of the 4 vertices directly.
                vs = np.array(verts)
                w_xy = float(vs[:, 0].max() - vs[:, 0].min())
                h_xy = float(vs[:, 1].max() - vs[:, 1].min())
                # If the face lies in Plane.XY (normal = +-Z), this is
                # the intended (w, h). For non-XY planes, fall through
                # to the polygon path.
                try:
                    n = _vec_np(face.normal_at())
                    n /= np.linalg.norm(n) + 1e-30
                    if abs(abs(n[2]) - 1.0) < 1e-6:
                        return RectProfile(w=w_xy, h=h_xy)
                except Exception:
                    pass

    pts3d = []
    for e in edges:
        if n_per_edge is None:
            n = _n_samples_for_edge(e, n_curved, n_straight)
        else:
            n = int(n_per_edge)
        n = max(2, n)
        for t in np.linspace(0.0, 1.0, n, endpoint=False):
            p = e.position_at(float(t))
            pts3d.append((p.X, p.Y, p.Z))

    pts3d = np.asarray(pts3d, dtype=float)

    if deduplicate and pts3d.shape[0] >= 2:
        keep = np.ones(pts3d.shape[0], dtype=bool)
        diffs = np.linalg.norm(np.diff(pts3d, axis=0, append=pts3d[:1]), axis=1)
        keep[diffs < 1e-12] = False
        pts3d = pts3d[keep]

    center = face.center()
    normal = face.normal_at()
    origin = np.array([center.X, center.Y, center.Z])
    n = np.array([normal.X, normal.Y, normal.Z])

    profile, _u, _v = polygon_profile_from_points_3d(
        pts3d, normal=n, origin=origin, u_axis=u_axis)
    return profile


def _edge_geom_name(edge):
    gt = getattr(edge, 'geom_type', '')
    return getattr(gt, 'name', str(gt)).upper()


def _vec_np(v):
    return np.array([v.X, v.Y, v.Z], dtype=float)


def coil_from_build123d_sweep(profile, path, current, u_axis=None,
                               n_curved=32, n_straight=2, plane='XY',
                               atol=1e-6):
    """Build a CoilBuilder from a (build123d Face profile, Wire path)
    pair -- the CAD-first "sweep spec" to PEEC one-liner.

    The profile is projected into a PolygonProfile. The path is walked
    edge-by-edge and mapped to CoilBuilder.add_straight / .add_arc calls.

    Supported edges in this MVP:
        LINE   -> add_straight(length)
        CIRCLE -> add_arc(radius, signed arc_angle)    (planar only)

    The path is assumed to lie in the plane given by ``plane`` (default
    'XY'). For non-planar paths or BSPLINE edges, raise
    NotImplementedError -- use coil_filaments_from_sweep (future) for
    those.

    Args:
        profile: build123d Face for the cross-section.
        path: build123d Wire, or a list of Edges in the correct order.
        current: total coil current [A].
        u_axis: optional u_axis override for profile_from_build123d_face.
        n_curved, n_straight: profile boundary sampling density.
        plane: 'XY' (only supported value currently). Selects the
            right-handed plane normal used to disambiguate arc signs.
        atol: tolerance [m] for checking path planarity.

    Returns:
        CoilBuilder (already set_start + set_profile + all segments
        appended), ready for to_filaments() / IHWorkpieceContext.evaluate().
    """
    from radia.coil_builder import CoilBuilder

    if plane != 'XY':
        raise NotImplementedError(f"plane='{plane}' not supported; use 'XY'.")

    prof = profile_from_build123d_face(
        profile, n_curved=n_curved, n_straight=n_straight, u_axis=u_axis)

    if hasattr(path, 'edges'):
        edges = list(path.edges())
    else:
        edges = list(path)
    if not edges:
        raise ValueError("path has no edges")

    # Planarity check: every endpoint must have |z| < atol.
    for e in edges:
        for t in (0.0, 1.0):
            p = _vec_np(e.position_at(t))
            if abs(p[2]) > atol:
                raise NotImplementedError(
                    f"Path not in the XY plane (z={p[2]:.3e} > atol={atol:.1e})."
                    " MVP supports planar XY paths only.")

    # Initial frame from first edge
    e0 = edges[0]
    p0 = _vec_np(e0.position_at(0.0))
    t0 = _vec_np(e0.tangent_at(0.0))
    norm_t0 = np.linalg.norm(t0)
    if norm_t0 < 1e-15:
        raise ValueError("Zero tangent at start of path.")
    y_axis = t0 / norm_t0
    z_axis = np.array([0.0, 0.0, 1.0])
    x_axis = np.cross(y_axis, z_axis)  # right-hand: X = Y x Z
    nx = np.linalg.norm(x_axis)
    if nx < 1e-12:
        raise ValueError("Start tangent parallel to plane normal.")
    x_axis /= nx
    # Rebuild z_axis to enforce orthonormality (robust to rounding)
    z_axis = np.cross(x_axis, y_axis)
    orientation = np.array([x_axis, y_axis, z_axis])

    coil = (CoilBuilder(current=current)
            .set_start(p0, orientation=orientation)
            .set_profile(prof))

    for e in edges:
        gt = _edge_geom_name(e)
        if gt == 'LINE':
            coil.add_straight(float(e.length))
            continue
        if gt == 'CIRCLE':
            radius = float(e.radius)
            if radius <= 0:
                raise ValueError(f"Non-positive circle radius: {radius}")
            # Signed arc_angle about +Z:
            # magnitude from length/radius; sign from tangent vs 90-deg CCW
            # rotation of (start - center). tangent_at(0) aligned with
            # CCW direction implies positive arc.
            c = _vec_np(e.arc_center)
            ps = _vec_np(e.position_at(0.0))
            ts = _vec_np(e.tangent_at(0.0))
            v0 = ps - c
            # Rotate v0 by +90 deg around +Z -> expected tangent if CCW
            rot_v0 = np.array([-v0[1], v0[0], 0.0])
            ccw = np.dot(rot_v0, ts) > 0
            mag_deg = float(np.degrees(e.length / radius))
            arc_angle = mag_deg if ccw else -mag_deg
            coil.add_arc(radius=radius, arc_angle=arc_angle)
            continue
        raise NotImplementedError(
            f"Edge geom_type '{gt}' not supported in "
            "coil_from_build123d_sweep. Supported: LINE, CIRCLE.")

    return coil
