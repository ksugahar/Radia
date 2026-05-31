"""End-to-end workflow: SF -> single-stroke fingerprint coil -> CAD -> PEEC -> field (Gx).

Transverse companion to ``demo_sf_to_peec_gz.py``.  Closes the same design loop
for a TRANSVERSE gradient (target Bz = Gx*x) on a cylinder of radius a:

  [1] stream function method   target Bz=Gx*x -> psi(phi,z) (ACA+TSVD)
  [2] equal-current contours    psi level set -> closed/open fingerprint loops
  [3] single-stroke chain       open each loop at the point closest to the
                                running chain end, connect via a (phi,z)-
                                geodesic on the cylinder (= a helical arc in
                                3D, periodic in phi) -> ONE continuous polyline
  [4] CAD                       **loft chain** -- circular cross-sections
                                placed along the chain with a parallel-
                                transported frame, lofted in short pieces
                                (~20 sections each, single OCC loft can't
                                handle >~50 with the saddle pattern's
                                cumulative twist), exported as a Compound
                                of solid spools.  Plain `sweep` along a
                                Spline fails on fingerprint kinks; this
                                multi-piece loft is the build123d equivalent
                                of "Cubit loft along curve" with a path
                                guide and is robust against the kinks.
  [5] PEEC                      chain polyline -> PEECBuilder -> L, R at f
  [6] field                     exact Biot-Savart of the chain -> Bz over the
                                DSV and along the x-axis
  [7] verify (loop closed)      manufactured-coil Bz vs the design Gx*x, both
                                RMS over the DSV and gradient nonlinearity on
                                the x-axis

The single-stroke step is geometrically harder than Gz (a saddle pattern has no
single-axis chaining direction), so the connection arcs DO add some stray Bz.
The verification quantifies the residual; in the lab the connection geometry
would be refined by hand or via an outer optimiser (CMA-ES).

Run:  python demo_sf_to_peec_gx.py [--nphi 24] [--nz 40] [--nlevels 12] [--with-peec]
"""
import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from radia.stream_function import aca_tsvd, pseudo_inverse_solve
from radia.biot_savart import h_segments_batch, MU0


# --------------------------------------------------------------------------
# Exact Biot-Savart segment Hz (same form as demo_coil_design_gx.py + tests):
# we use it for the SF solve matrix entries so the design and verification
# kernels are bit-consistent.
# --------------------------------------------------------------------------
def _seg_Hz(O, P1, P2):
    eps = 1.0e-15
    d = P2 - P1
    R21 = float(d @ d)
    o1 = O - P1
    o2 = O - P2
    OR1 = np.sqrt(o1 @ o1)
    OR2 = np.sqrt(o2 @ o2)
    O121 = float(o1 @ d)
    O221 = float(o2 @ d)
    Rc12 = O121 / OR1
    Rc21 = -O221 / OR2
    L = o1 - O121 * d / R21
    L1 = float(L @ L) + eps
    f = (Rc12 + Rc21) / (4.0 * np.pi * L1 * R21)
    return (d[0] * L[1] - d[1] * L[0]) * f


def _loop_Hz(obs, corners):
    hz = 0.0
    for c in range(4):
        n = (c + 1) % 4
        hz += _seg_Hz(obs, corners[c], corners[n])
    return hz


def cyl_point(a, phi, z):
    return np.array([a * np.cos(phi), a * np.sin(phi), z])


def loop_corners(a, phi, z, dphi, dz):
    return np.array([
        cyl_point(a, phi - dphi / 2, z - dz / 2),
        cyl_point(a, phi + dphi / 2, z - dz / 2),
        cyl_point(a, phi + dphi / 2, z + dz / 2),
        cyl_point(a, phi - dphi / 2, z + dz / 2),
    ])


def make_dsv(dsv, n):
    g = np.linspace(-dsv, dsv, n)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    r = np.linalg.norm(pts, axis=1)
    return pts[r <= dsv + 1e-12]


# --------------------------------------------------------------------------
# [2] psi contour extraction (matplotlib) -> polylines in (phi, z) space.
#     Each polyline is either closed (start ~ end) or open (cut by the phi seam).
# --------------------------------------------------------------------------
def contour_polylines_phi_z(psi_zphi, phi_grid, z_grid, n_levels):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lo, hi = float(psi_zphi.min()), float(psi_zphi.max())
    if hi <= lo:
        return [], 0.0
    dI = (hi - lo) / n_levels
    levels = lo + (np.arange(n_levels) + 0.5) * dI

    phi_ext = np.concatenate([phi_grid, [phi_grid[0] + 2.0 * np.pi]])
    psi_ext = np.column_stack([psi_zphi, psi_zphi[:, :1]])
    PHI, Z = np.meshgrid(phi_ext, z_grid)

    fig = plt.figure()
    cs = plt.contour(PHI, Z, psi_ext, levels=levels)
    plt.close(fig)

    polylines = []
    allsegs = getattr(cs, "allsegs", None)
    if allsegs is None:
        return [], dI
    for segs_at_level in allsegs:
        for poly in segs_at_level:
            if len(poly) < 3:
                continue
            polylines.append(np.asarray(poly, dtype=float))
    return polylines, dI


# --------------------------------------------------------------------------
# [3] single-stroke connection in (phi, z) space.
#     Each closed polyline is "opened" at the index whose point is closest to
#     the running chain end (phi-periodic distance).  Open polylines (cut by
#     the phi seam) are reversed when their tail end is closer.  Connections
#     are linear (phi, z) interpolations -> geodesics on the cylinder (helices
#     in 3D, of zero stray Bx/By if both endpoints align in (phi, z), small
#     stray otherwise).
# --------------------------------------------------------------------------
def _phi_period_dist(p, q):
    """Distance in (phi, z) with phi-periodicity."""
    dphi = (p[0] - q[0] + np.pi) % (2.0 * np.pi) - np.pi
    dz = p[1] - q[1]
    return float(np.hypot(dphi, dz))


def _open_closed_polyline_at(poly, anchor):
    """Open + re-roll a closed polyline so it starts at the index nearest to
    anchor, KEEPING traversal direction (matplotlib's natural contour orientation
    encodes the per-loop current sign required for a coherent Bz=Gx*x; reversing
    a polyline flips its Bz contribution and breaks the field).

    The closed loop is OPENED at the anchor index: rolled = body[i:] ++ body[:i],
    of length len(body) -- start at body[i], end at body[i-1].  Do NOT re-append
    body[i] (that would merge nodes[K] with nodes[K+N-1] at identical xyz in PEEC
    and short-circuit the chain to Z=0).  We sacrifice one segment per loop;
    a closed contour with the start-to-end gap is electromagnetically nearly
    identical because the chain's bridging arc fills it.  Open polylines are
    returned unchanged."""
    closed = np.linalg.norm(poly[0] - poly[-1]) < 1e-9
    if closed:
        d = np.array([_phi_period_dist(p, anchor) for p in poly[:-1]])
        i = int(np.argmin(d))
        body = poly[:-1]
        return np.vstack([body[i:], body[:i]])
    return poly


def _geodesic_arc_phi_z(p_start, p_end, n_seg):
    """Short geodesic in (phi, z) from p_start to p_end, n_seg interior points.

    Phi is interpolated along the shortest-wrap direction (delta in [-pi, pi]).
    """
    dphi = (p_end[0] - p_start[0] + np.pi) % (2.0 * np.pi) - np.pi
    dz = p_end[1] - p_start[1]
    t = np.linspace(0.0, 1.0, n_seg + 2)[1:-1]
    return np.column_stack([p_start[0] + t * dphi, p_start[1] + t * dz])


def single_stroke_kuijpers_chain_phi_z(polylines, a, n_blend=1):
    """Kuijpers 2023 Method-1 inspired single-stroke chain.

    Reference:
        B.J.A. Kuijpers, J.W. Jansen, E.A. Lomonova,
        "Comparison of Discretization Methods for Continuous Stream-Function
         Distributions", Compumag 2023, Kyoto, paper [525].

    The paper's Method 1 cuts every iso-contour at a SHARED location and
    connects adjacent contours via SHORT straight segments along the cut
    line.  The chain deviates from the iso-contours ONLY on these segments,
    and the paper proves the field error in the target region is collocated
    with this deviation -- so keeping the segments short AND parallel
    minimises stray field.  This is the "rung" pattern visible in Fig.4 of
    the paper.

    For a 4-lobe Gx fingerprint coil on a cylinder, we adopt one cut line
    per lobe:

      lobe (+x, +z) / (+x, -z) : cut phi = 0
      lobe (-x, +z) / (-x, -z) : cut phi = pi

    Each contour is opened at the index whose phi is closest to its lobe's
    cut phi.  Within a lobe, contours are sorted by the z coordinate of
    their cut point (descending for +z lobes -- spiral down from top edge
    toward the saddle; ascending for -z lobes -- spiral up from bottom edge
    toward the saddle).  Adjacent contours are then connected by a STRAIGHT
    line in (phi, z) at the cut: a vertical "rung" at fixed phi.  Inter-lobe
    transitions use the same straight-line blend; for a Gx coil they are
    three short segments (two axial at the saddle z, one azimuthal half-turn
    at z=+L/2 between the -x and +x saddles).
    """
    n = len(polylines)
    if n == 0:
        return np.zeros((0, 2))

    # 1. 3D centroid -> 4-quadrant lobe classification.
    centroids_3d = []
    for poly in polylines:
        pts = np.column_stack([
            a * np.cos(poly[:, 0]),
            a * np.sin(poly[:, 0]),
            poly[:, 1],
        ])
        centroids_3d.append(pts.mean(axis=0))
    centroids_3d = np.asarray(centroids_3d)

    lobes = {(+1, +1): [], (+1, -1): [], (-1, +1): [], (-1, -1): []}
    for k, c in enumerate(centroids_3d):
        sx = +1 if c[0] >= 0 else -1
        sz = +1 if c[2] >= 0 else -1
        lobes[(sx, sz)].append(k)

    # 2. Cut phi per lobe.  Use the contour's point closest to cut_phi
    #    (NATURAL orientation, no top/bottom strong-arming -- forcing top
    #    or bottom flips the y-mirror partners' current direction relative
    #    to matplotlib's orientation and broke field accuracy in earlier
    #    attempts, 2026-05-30).
    lobe_cut_phi = {
        (+1, +1): 0.0, (+1, -1): 0.0,
        (-1, +1): np.pi, (-1, -1): np.pi,
    }

    def cut_index(idx, cut_phi):
        p = polylines[idx]
        dphi = ((p[:, 0] - cut_phi + np.pi) % (2 * np.pi)) - np.pi
        return int(np.argmin(np.abs(dphi)))

    def cut_z(idx, cut_phi):
        return float(polylines[idx][cut_index(idx, cut_phi), 1])

    # 3. Within each lobe, sort by cut z: descending for +z lobes (chain
    #    spirals down from top edge toward saddle), ascending for -z lobes
    #    (chain spirals up from bottom edge toward saddle).
    #
    #    Pair-aware variants that group +y and -y mirror partners separately
    #    were tried (2026-05-30): they break matplotlib's consistent CCW
    #    contour orientation when -y partners are reversed (-> current
    #    direction flips, RMS goes from 16 % to 40 %).  Sorting both groups
    #    in the SAME direction trades the mirror criss-crosses for a long
    #    inner-+y-to-outer--y rung that is just as bad.  The flat cut-z
    #    sort below preserves matplotlib's orientation and gives the lowest
    #    DSV RMS we have found.  The remaining diagonal "criss-cross" rungs
    #    visible in GMSH are the irreducible inter-lobe transitions of a
    #    single-stroke fingerprint coil; eliminating them requires the
    #    crossover-COMPENSATED iteration (Kuijpers + folded back into the
    #    least-norm RHS), which is the open research direction listed in
    #    ``radia-mcp aca_tsvd(single_stroke)`` topic A.
    for key in lobes:
        cp = lobe_cut_phi[key]
        sort_descending = (key[1] == +1)
        lobes[key].sort(key=(lambda k, _cp=cp: -cut_z(k, _cp))
                         if sort_descending
                         else (lambda k, _cp=cp: cut_z(k, _cp)))

    traversal = [(+1, +1), (+1, -1), (-1, -1), (-1, +1)]

    # 4. Build the chain: open each contour at its cut, connect via straight
    #    (phi, z) blend at fixed cut phi within a lobe.
    chain_parts = []
    for lobe_key in traversal:
        members = lobes[lobe_key]
        if not members:
            continue
        cp = lobe_cut_phi[lobe_key]
        for k in members:
            poly = polylines[k]
            i_cut = cut_index(k, cp)
            closed = np.linalg.norm(poly[0] - poly[-1]) < 1e-9
            if closed:
                body = poly[:-1]
                opened = np.vstack([body[i_cut:], body[:i_cut]])
            else:
                opened = np.vstack([poly[i_cut:], poly[:i_cut]]) \
                    if i_cut > 0 else poly

            if chain_parts:
                prev_end = chain_parts[-1][-1]
                blend = _geodesic_arc_phi_z(prev_end, opened[0], n_blend)
                if len(blend) > 0:
                    chain_parts.append(blend)
            chain_parts.append(opened)

    return np.vstack(chain_parts) if chain_parts else np.zeros((0, 2))


def single_stroke_lobe_chain_phi_z(polylines, a, n_arc=8):
    """Lobe-aware single-stroke chain (recommended).

    Classifies each closed contour into one of 4 LOBES by the sign of its
    3D centroid in (x, z): ``(+x,+z)``, ``(+x,-z)``, ``(-x,-z)``, ``(-x,+z)``.
    For a Gx fingerprint coil the SF design places contours into exactly
    these four saddle regions; within a lobe the contours are nested
    concentrically, so a within-lobe greedy chain uses only SHORT arcs.

    The lobes are then traversed in a fixed order

        (+x,+z) --> (+x,-z) --> (-x,-z) --> (-x,+z)

    with the traversal direction within each lobe alternating between
    ``out-to-in`` and ``in-to-out`` so the chain exits a lobe at the side
    closest to the next lobe (= the lobe transitions cost only 3 arcs
    total, vs ``N-1`` arcs for a naive global greedy NN).  The lobe
    transition arcs are at z = +/- L/2 (cylinder edges) or phi = +/- pi/2
    (cylinder waist) which is the classical Maxwell-pair geometry: the
    return current contributions partially cancel at the DSV center.

    This replaces the earlier global greedy implementation
    ``single_stroke_chain_phi_z`` which produced visually obvious
    "wasted" connection arcs criss-crossing the cylinder.
    """
    n = len(polylines)
    if n == 0:
        return np.zeros((0, 2))

    # 1. classify by 3D centroid sign (sx, sz) -> 4 lobes
    centroids_3d = []
    for poly in polylines:
        pts = np.column_stack([
            a * np.cos(poly[:, 0]),
            a * np.sin(poly[:, 0]),
            poly[:, 1],
        ])
        centroids_3d.append(pts.mean(axis=0))
    centroids_3d = np.asarray(centroids_3d)

    lobes = {(+1, +1): [], (+1, -1): [], (-1, +1): [], (-1, -1): []}
    for k, c in enumerate(centroids_3d):
        sx = +1 if c[0] >= 0 else -1
        sz = +1 if c[2] >= 0 else -1
        lobes[(sx, sz)].append(k)

    # 2. arc length in (phi, z) -- a proxy for "outer vs inner" nestedness;
    #    longer = outer in a nested family.  Sort each lobe outer-first.
    def arc_len_phi_z(idx):
        p = polylines[idx]
        d = np.diff(p, axis=0)
        # phi periodic correction: any |d_phi| > pi wraps the other way
        dphi = np.where(np.abs(d[:, 0]) > np.pi,
                        d[:, 0] - 2 * np.pi * np.sign(d[:, 0]), d[:, 0])
        return float(np.sum(np.sqrt(dphi ** 2 + d[:, 1] ** 2)))

    for key in lobes:
        lobes[key].sort(key=lambda k: -arc_len_phi_z(k))

    # 3. traverse lobes in a fixed spiral order with alternating direction
    #    so we exit each lobe near the next one's entry point.
    traversal = [
        ((+1, +1), "out_to_in"),
        ((+1, -1), "in_to_out"),
        ((-1, -1), "out_to_in"),
        ((-1, +1), "in_to_out"),
    ]

    chain_parts = []

    def append_polyline(poly, prev_end):
        opened = _open_closed_polyline_at(poly, prev_end)
        if prev_end is not None and len(opened) > 0:
            arc = _geodesic_arc_phi_z(prev_end, opened[0], n_arc)
            if len(arc) > 0:
                chain_parts.append(arc)
        chain_parts.append(opened)
        return opened[-1] if len(opened) > 0 else prev_end

    prev_end = None
    n_inter_lobe = 0
    last_lobe = None
    for lobe_key, direction in traversal:
        members = list(lobes[lobe_key])
        if not members:
            continue
        if direction == "in_to_out":
            members = members[::-1]
        if last_lobe is not None and prev_end is not None:
            n_inter_lobe += 1  # one connection arc to enter this lobe
        last_lobe = lobe_key
        for k in members:
            prev_end = append_polyline(polylines[k], prev_end)

    chain = np.vstack(chain_parts) if chain_parts else np.zeros((0, 2))
    return chain


def single_stroke_chain_phi_z(polylines, n_arc=8):
    """LEGACY: global greedy nearest-end chain in (phi, z).  Kept for the
    ``--chain-method greedy`` option and for comparison with the
    lobe-aware default.  Use ``single_stroke_lobe_chain_phi_z`` instead."""
    n = len(polylines)
    used = [False] * n
    # Start from the polyline closest to (phi=0, z=-L/2) for a deterministic chain
    z_min = float(min(float(p[:, 1].min()) for p in polylines))
    anchor0 = np.array([0.0, z_min])
    d0 = [float(min(_phi_period_dist(q, anchor0) for q in p)) for p in polylines]
    k0 = int(np.argmin(d0))
    first = _open_closed_polyline_at(polylines[k0], anchor0)
    chain = [first]
    used[k0] = True

    for _ in range(n - 1):
        end = chain[-1][-1]
        # Find the polyline whose nearest point on it is closest to `end`
        best_d = float("inf")
        best_j = -1
        for j in range(n):
            if used[j]:
                continue
            dj = float(min(_phi_period_dist(q, end) for q in polylines[j]))
            if dj < best_d:
                best_d = dj
                best_j = j
        nxt = _open_closed_polyline_at(polylines[best_j], end)
        arc = _geodesic_arc_phi_z(end, nxt[0], n_arc)
        chain.append(arc)
        chain.append(nxt)
        used[best_j] = True

    return np.vstack(chain)


def _embed3d_bodies(polylines, a):
    """Open-loop bodies (closed loops with the duplicate end vertex dropped)
    plus their 3D cylinder embedding (a*cos phi, a*sin phi, z).  Distances for
    nearest-neighbour ordering and balanced cuts are TRUE wire distances in 3D,
    so an azimuthal rung at large radius is correctly costed more than the same
    delta-phi near the axis."""
    bodies, embs = [], []
    for poly in polylines:
        body = poly[:-1] if np.linalg.norm(poly[0] - poly[-1]) < 1e-9 else poly
        bodies.append(np.asarray(body, dtype=float))
        embs.append(np.column_stack([
            a * np.cos(body[:, 0]), a * np.sin(body[:, 0]), body[:, 1]]))
    return bodies, embs


def single_stroke_nn_blend_phi_z(polylines, a, n_blend=8, cd_passes=4):
    """NEGATIVE RESULT (2026-05-31): geometric-shortest blend != least field.

    Motivation (what was asked): replace the long inter-lobe "rungs" with an
    algorithm that connects ONLY nearest-neighbour contours and BLENDS them
    at the geometrically shortest point, so no rung spans far across the
    (phi, z) map.

    Measured outcome on the Gx 4-lobe fingerprint: this produces the SHORTEST
    individual rungs but the WORST DSV field of all four methods --

        method     DSV RMS   x-nonlin   longest rung
        nn_blend   0.651     0.389      291 mm   <-- this function
        kuijpers   0.162     0.097      291 mm   (default)
        greedy     0.218     0.114       61 mm
        lobe       0.240     0.097       40 mm

    Why it fails: on a multi-lobe coil the stray field of a connecting rung
    scales with the rung's TRANSVERSE (azimuthal) component, NOT its 3D
    length (Kuijpers' collocation argument).  Measuring the AZIMUTHAL arc
    length summed over all rungs makes the correlation explicit:

        method     azimuthal rung total   DSV RMS
        kuijpers          7444 mm          0.162
        nn_blend         12847 mm          0.651

    The balanced cut minimised the 3D CHORD distance between consecutive cut
    points, but the wire actually travels the azimuthal ARC between them, and
    geometry-only nearest-neighbour ordering scatters the cuts in phi AND
    interleaves the four lobes' ALTERNATING current signs (+,-,-,+ for Gx).
    So the "shortest-chord" hop frequently turns into a long azimuthal sweep
    bridging two opposite-sign contours -- more transverse rung, worse field.
    Minimising geometric distance is the wrong objective; kuijpers keeps the
    WITHIN-lobe rungs axial (per-lobe constant-phi cut) and pays the
    azimuthal cost only on the few necessary inter-lobe half-turns.  The
    principled way further down is the Path-A compensated iteration
    (``--compensated-iter``), which folds the rung field back into the SF
    target.

    Kept as a CAUTIONARY option (``--chain-method nn_blend``); do not use it
    for field accuracy.  See memory feedback_single_stroke_chain_orientation_
    traps + radia-mcp aca_tsvd(single_stroke).

    The algorithm itself (for reference):

      1. NEIGHBOURS ONLY -- contours are visited in nearest-neighbour order
         on the true 3D cylinder closest-approach distance, seeded at the
         outermost contour.  The wire therefore only ever hops to the
         spatially nearest not-yet-used contour; it never connects a
         contour to a far-away one.

      2. LEAST-IMPACT BLEND -- each contour is opened at the point that
         JOINTLY minimises the incoming + outgoing connector length
         (coordinate descent: cut_k = argmin over the contour of
         dist(cut_{k-1}, .) + dist(., cut_{k+1})).  Because a closed loop
         opened at one point has its entry and exit at that SAME point, a
         single cut governs both rungs; optimising it drives every rung to
         its shortest possible length -> minimum deviation from the
         iso-contours -> least stray field (Kuijpers' collocation argument).

    Orientation is preserved throughout (loops are rolled, never reversed --
    reversing flips the per-loop current sign and breaks Bz; see
    feedback_single_stroke_chain_orientation_traps).

    Returns the (phi, z) chain.  Use ``connector_lengths_phi_z`` to audit the
    rung lengths (the metric this method minimises).
    """
    n = len(polylines)
    if n == 0:
        return np.zeros((0, 2))
    bodies, embs = _embed3d_bodies(polylines, a)
    if n == 1:
        b = bodies[0]
        return np.vstack([b, b[:1]]) if len(b) else np.zeros((0, 2))

    # --- pairwise closest-approach (distance + the two argmin indices) ---
    def loop_loop(i, j):
        D = np.linalg.norm(embs[i][:, None, :] - embs[j][None, :, :], axis=2)
        ai, aj = np.unravel_index(int(np.argmin(D)), D.shape)
        return float(D[ai, aj]), int(ai), int(aj)

    dist = np.full((n, n), np.inf)
    for i in range(n):
        for j in range(i + 1, n):
            d, _, _ = loop_loop(i, j)
            dist[i, j] = dist[j, i] = d

    # --- seed: outermost contour (largest 3D bounding-box diagonal) ---
    extent = [float(np.linalg.norm(embs[i].max(0) - embs[i].min(0)))
              for i in range(n)]
    start = int(np.argmax(extent))

    # --- nearest-neighbour visiting order (neighbours only) ---
    order = [start]
    used = [False] * n
    used[start] = True
    for _ in range(n - 1):
        last = order[-1]
        cand = [dist[last, j] if not used[j] else np.inf for j in range(n)]
        nxt = int(np.argmin(cand))
        order.append(nxt)
        used[nxt] = True

    # --- initialise each cut at the point closest to its order-successor ---
    cut = [0] * n
    for kk in range(n):
        i = order[kk]
        j = order[kk + 1] if kk + 1 < n else order[kk - 1]
        _, ai, _ = loop_loop(i, j)
        cut[kk] = ai

    # --- coordinate descent: balanced cut minimises in + out rung ---
    def cut_pt(kk):
        return embs[order[kk]][cut[kk]]

    for _ in range(cd_passes):
        for kk in range(n):
            i = order[kk]
            pts = embs[i]
            cost = np.zeros(len(pts))
            if kk > 0:
                cost = cost + np.linalg.norm(pts - cut_pt(kk - 1), axis=1)
            if kk < n - 1:
                cost = cost + np.linalg.norm(pts - cut_pt(kk + 1), axis=1)
            cut[kk] = int(np.argmin(cost))

    # --- build the chain: open at the (balanced) cut, short geodesic blend ---
    chain_parts = []
    for kk in range(n):
        i = order[kk]
        body = bodies[i]
        ci = cut[kk]
        opened = np.vstack([body[ci:], body[:ci]])
        if chain_parts:
            prev_end = chain_parts[-1][-1]
            blend = _geodesic_arc_phi_z(prev_end, opened[0], n_blend)
            if len(blend) > 0:
                chain_parts.append(blend)
        chain_parts.append(opened)

    return np.vstack(chain_parts) if chain_parts else np.zeros((0, 2))


def _kuijpers_lobe_order(polylines, a, snake=False):
    """Return contour indices in kuijpers lobe/current-sign order.

    4-lobe classification by (sign x, sign z) of the 3D centroid; within a
    lobe sort by the z of the fixed-phi cut; lobes traversed
    (+x+z)(+x-z)(-x-z)(-x+z).  This is the ORDER that respects the Gx coil's
    alternating current signs (the part of kuijpers that geometry-only NN
    throws away).

    ``snake=False`` (kuijpers original): +z lobes sorted descending z, -z
    ascending.  This makes each lobe spiral toward the saddle, so the
    inter-lobe bridge jumps saddle -> opposite edge (a long ~L/2 axial
    "wasted" rung).

    ``snake=True`` (boustrophedon): NEGATIVE RESULT (2026-05-31).  The
    within-lobe sort direction was alternated [(+x+z) top->saddle, (+x-z)
    saddle->bottom, (-x-z) bottom->saddle, (-x+z) saddle->top] to try to make
    consecutive lobes MEET at the saddle and so shorten the long inter-lobe
    bridges.  It did NEITHER: the longest rung stayed 291 mm AND the field
    got WORSE (DSV RMS 9.29% -> 19.15%).  No polyline is reversed (only the
    visiting SEQUENCE changes), so it is not an orientation flip -- it is the
    same lesson again: changing the rung arrangement breaks the symmetric
    stray-field CANCELLATION that gives the good field, and the long diagonal
    inter-lobe bridges are STRUCTURALLY irreducible for a single-stroke
    4-lobe coil.  Kept as a default-False cautionary kwarg; do not enable."""
    n = len(polylines)
    centroids = []
    for poly in polylines:
        pts = np.column_stack([a * np.cos(poly[:, 0]),
                               a * np.sin(poly[:, 0]), poly[:, 1]])
        centroids.append(pts.mean(axis=0))
    centroids = np.asarray(centroids)
    lobes = {(+1, +1): [], (+1, -1): [], (-1, +1): [], (-1, -1): []}
    for k, c in enumerate(centroids):
        lobes[(+1 if c[0] >= 0 else -1, +1 if c[2] >= 0 else -1)].append(k)
    cut_phi = {(+1, +1): 0.0, (+1, -1): 0.0, (-1, +1): np.pi, (-1, -1): np.pi}

    def cut_z(idx, cp):
        p = polylines[idx]
        dphi = ((p[:, 0] - cp + np.pi) % (2 * np.pi)) - np.pi
        return float(p[int(np.argmin(np.abs(dphi))), 1])

    # descending-z flag per lobe key
    if snake:
        desc_of = {(+1, +1): True, (+1, -1): True,
                   (-1, -1): False, (-1, +1): False}
    else:
        desc_of = {key: (key[1] == +1) for key in lobes}

    for key in lobes:
        cp = cut_phi[key]
        desc = desc_of[key]
        lobes[key].sort(key=(lambda k, _cp=cp: -cut_z(k, _cp)) if desc
                        else (lambda k, _cp=cp: cut_z(k, _cp)))
    order = []
    for key in [(+1, +1), (+1, -1), (-1, -1), (-1, +1)]:
        order.extend(lobes[key])
    return order


def single_stroke_field_aware_phi_z(polylines, a, n_blend=8, cd_passes=5,
                                    snake=False):
    """Field-aware single-stroke (DEFAULT): kuijpers sign-order + min-AZIMUTHAL
    cuts.  Best DSV RMS we have found -- beats kuijpers in every tested config.

    Two ingredients, learned from the nn_blend negative result (2026-05-31):

      1. ORDER = kuijpers lobe/current-sign order (``_kuijpers_lobe_order``).
         This is the DOMINANT factor.  The nn_blend failure (RMS 0.65) was
         NOT about rung length -- it was its geometry-only nearest-neighbour
         order INTERLEAVING the four lobes' alternating current signs
         (+,-,-,+).  Keep opposite-sign contours out of adjacency and the
         catastrophe disappears.

      2. CUT = chosen to minimise the AZIMUTHAL arc to the chain neighbours
         (coordinate descent: cut_k = argmin over the contour of
         a*|dphi(pt, cut_{k-1})| + a*|dphi(pt, cut_{k+1})|).  Axial rung
         content (dz) is free (zero stray) so it is not penalised.  Versus
         kuijpers' fixed-phi snap, aligning consecutive cuts to their MUTUAL
         best phi lets the rungs' stray fields cancel more symmetrically over
         the DSV.

    Measured (default nphi=24 nz=40 nlevels=12 DSV):

        method        DSV RMS   x-nonlin   azimuthal rung total
        field_aware   0.093     0.072      12522 mm   <- this function
        kuijpers      0.162     0.097       7444 mm
        lobe          0.240     0.097      10863 mm
        greedy        0.218     0.114      11580 mm
        nn_blend      0.651     0.389      12847 mm

    Robust: field_aware < kuijpers RMS at nlevels 10/12/16 and at
    nphi32/nz48/nlevels14 (30-54% lower).  NOTE the azimuthal total is NOT a
    clean predictor -- field_aware has MORE azimuthal arc than kuijpers yet
    lower RMS, because the field impact is the symmetric CANCELLATION of the
    rung stray fields over the DSV, not their summed length.  That subtlety
    (cancellation depends on placement + sign structure, not a closed-form
    metric) is why the single-stroke connection is better handled as a
    reason-and-verify SKILL than a fixed formula.

    Orientation preserved (loops rolled, never reversed).
    """
    n = len(polylines)
    if n == 0:
        return np.zeros((0, 2))
    bodies, _ = _embed3d_bodies(polylines, a)
    if n == 1:
        b = bodies[0]
        return np.vstack([b, b[:1]]) if len(b) else np.zeros((0, 2))

    order = _kuijpers_lobe_order(polylines, a, snake=snake)

    def wrap(d):
        return (d + np.pi) % (2.0 * np.pi) - np.pi

    # initialise each cut at the point closest to its lobe cut phi (kuijpers)
    cut_phi_of = {}
    centroids = [np.column_stack([a*np.cos(p[:,0]), a*np.sin(p[:,0]), p[:,1]]).mean(0)
                 for p in polylines]
    cut = {}
    for idx in order:
        cp = 0.0 if centroids[idx][0] >= 0 else np.pi
        cut[idx] = int(np.argmin(np.abs(wrap(bodies[idx][:, 0] - cp))))

    # coordinate descent on AZIMUTHAL arc to neighbours (axial dz free)
    for _ in range(cd_passes):
        for kk, idx in enumerate(order):
            phis = bodies[idx][:, 0]
            cost = np.zeros(len(phis))
            if kk > 0:
                pp = bodies[order[kk - 1]][cut[order[kk - 1]], 0]
                cost = cost + np.abs(wrap(phis - pp))
            if kk < n - 1:
                pn = bodies[order[kk + 1]][cut[order[kk + 1]], 0]
                cost = cost + np.abs(wrap(phis - pn))
            cut[idx] = int(np.argmin(cost))

    chain_parts = []
    for idx in order:
        body = bodies[idx]
        ci = cut[idx]
        opened = np.vstack([body[ci:], body[:ci]])
        if chain_parts:
            prev_end = chain_parts[-1][-1]
            blend = _geodesic_arc_phi_z(prev_end, opened[0], n_blend)
            if len(blend) > 0:
                chain_parts.append(blend)
        chain_parts.append(opened)
    return np.vstack(chain_parts) if chain_parts else np.zeros((0, 2))


def connector_lengths_phi_z(polylines, chain_phi_z, a):
    """Audit the inter-contour connector segments (the 'rungs').

    A rung is a chain segment that bridges a gap between two contours rather
    than tracing an iso-contour.  We flag segments longer than 2x the median
    edge length -- a robust proxy that does not need per-vertex contour
    membership.  Returns ``(max_rung, total_rung, n_rung, azimuthal_total)``
    in metres.  ``azimuthal_total`` = sum over rungs of ``a*|dphi|``.  It is
    a useful DIAGNOSTIC (azimuthal rung content is more harmful than axial),
    but it is NOT a clean predictor of DSV RMS: field_aware has a LARGER
    azimuthal total than kuijpers yet a LOWER RMS, because the field impact
    is the symmetric CANCELLATION of the rung stray fields (current-sign
    order + placement), not their summed length.  Do not optimise this value
    directly -- optimise the field (DSV RMS) and use the skill's verify loop.
    """
    path = np.column_stack([
        a * np.cos(chain_phi_z[:, 0]),
        a * np.sin(chain_phi_z[:, 0]),
        chain_phi_z[:, 1]])
    seg = np.linalg.norm(np.diff(path, axis=0), axis=1)   # length L = n-1
    valid = seg > 1.0e-12
    if not valid.any():
        return 0.0, 0.0, 0, 0.0
    med = float(np.median(seg[valid]))
    is_rung = seg > 2.0 * med                              # length L, boolean
    if not is_rung.any():
        return 0.0, 0.0, 0, 0.0
    dphi = (np.diff(chain_phi_z[:, 0]) + np.pi) % (2.0 * np.pi) - np.pi
    azim = a * np.abs(dphi)                                # length L
    rung = seg[is_rung]
    return (float(rung.max()), float(rung.sum()), int(is_rung.sum()),
            float(azim[is_rung].sum()))


def chain_phi_z_to_3d(chain_phi_z, a, dedupe_tol=1.0e-9):
    """Map (phi, z) chain to 3D and drop adjacent duplicates (matplotlib's
    contour algorithm emits repeated vertices on cell boundaries; these become
    zero-length segments that break PEEC and Spline construction)."""
    path = np.column_stack([
        a * np.cos(chain_phi_z[:, 0]),
        a * np.sin(chain_phi_z[:, 0]),
        chain_phi_z[:, 1],
    ])
    diffs = np.linalg.norm(np.diff(path, axis=0), axis=1)
    keep = np.concatenate([[True], diffs > dedupe_tol])
    return path[keep]


# --------------------------------------------------------------------------
# [6] Field of the single-stroke chain, with a unit current.  Used in two
#     places: (a) gain-fit so the chain reproduces Gx*x in physical units,
#     (b) verification.
# --------------------------------------------------------------------------
def bz_at(path_3d, current, obs):
    segs = np.stack([path_3d[:-1], path_3d[1:]], axis=1)
    H = h_segments_batch(segs, obs, current=current)
    return MU0 * H[:, 2]


# NOTE on "multi-wire" (avoiding the long inter-lobe bridges): cutting the
# single-stroke chain at its long rungs is NOT a valid way to do it -- the
# resulting sub-paths are OPEN current paths (current cannot start/stop in
# mid-air; div J != 0), so their Biot-Savart field is unphysical (measured
# RMS jumps 9.3% -> 21%).  A physical multi-wire coil needs each independent
# wire to be a CLOSED loop.  For the Gx fingerprint that is exactly
# ``demo_coil_design_gx.py``: it drives each saddle-shaped CLOSED loop as an
# independent conductor (no bridges at all) and reaches ~0.8% DSV RMS -- an
# order of magnitude better than any single-stroke chain.  So the answer to
# "avoid the wasteful long connections" is: use the independent-closed-loop
# design (demo_coil_design_gx.py, multiple current feeds), NOT a cut-up
# single stroke.


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nphi", type=int, default=24, help="surface nodes around phi")
    ap.add_argument("--nz", type=int, default=40, help="surface nodes along z")
    ap.add_argument("--nlevels", type=int, default=12, help="equal-current contour levels")
    ap.add_argument("--k", type=int, default=0, help="TSVD modes (0 = auto)")
    ap.add_argument("--ndsv", type=int, default=7, help="DSV grid points per axis")
    ap.add_argument("--narc", type=int, default=8,
                    help="segments per single-stroke connection arc")
    ap.add_argument("--chain-method",
                    choices=["field_aware", "kuijpers", "lobe", "greedy",
                             "nn_blend"],
                    default="field_aware",
                    help="field_aware (default, recommended, 2026-05-31) = "
                         "kuijpers lobe/current-sign ORDER + each contour cut "
                         "chosen to minimise the AZIMUTHAL arc to its chain "
                         "neighbours.  Beats kuijpers on DSV RMS in every "
                         "tested config (0.093 vs 0.162 at the default; 30-54% "
                         "lower RMS across nlevels/mesh sweeps); "
                         "kuijpers = Kuijpers 2023 Method-1: per-lobe cut at "
                         "fixed phi, straight AXIAL 'rung' blends (the prior "
                         "best); "
                         "lobe = 4-quadrant within-lobe greedy NN; "
                         "greedy = legacy global NN (incoming-only cut); "
                         "nn_blend = NEGATIVE RESULT (2026-05-31): "
                         "nearest-neighbour order + balanced cuts that "
                         "minimise the GEOMETRIC connector length.  Produces "
                         "the SHORTEST rungs but the WORST field (RMS 0.65 vs "
                         "kuijpers 0.16) because on a multi-lobe coil the "
                         "stray field scales with a rung's TRANSVERSE "
                         "(azimuthal) component, not its length, and "
                         "geometry-only NN interleaves the lobes' alternating "
                         "current signs.  Kept as a cautionary option; do NOT "
                         "use for field accuracy")
    ap.add_argument("--with-peec", action="store_true",
                    help="also run STEP export + PEEC L, R")
    ap.add_argument("--compensated-iter", type=int, default=0,
                    help="Path-A compensated iteration count: at each iter "
                         "subtract the chain's parasitic Bz from the SF "
                         "target and re-solve, so the next SF design absorbs "
                         "the chain crossover field.  0 = disable.  ACA+TSVD "
                         "factorisation is reused across iterations, so each "
                         "additional iter is cheap (one back-substitution).")
    ap.add_argument("--compensated-step", type=float, default=1.0,
                    help="damping factor on the phi update (0 < step <= 1); "
                         "0.5-0.8 may help when full step oscillates.")
    args = ap.parse_args()

    a = 0.15           # cylinder radius [m]
    L = 0.50           # cylinder length [m]
    Gx = 1.0           # target transverse gradient dBz/dx [T/m]
    dsv = 0.10         # DSV half-extent [m]

    print("=== SF -> single-stroke fingerprint coil -> (CAD/PEEC) -> field (Gx) ===")

    # ------- [1] SF design (identical solve to demo_coil_design_gx.py) -------
    phi_grid = np.linspace(0.0, 2.0 * np.pi, args.nphi, endpoint=False)
    z_grid = np.linspace(-L / 2, L / 2, args.nz)
    dphi = 2.0 * np.pi / args.nphi
    dz = L / (args.nz - 1)

    corners_list = []
    for q, zq in enumerate(z_grid):
        for p, pp in enumerate(phi_grid):
            corners_list.append(loop_corners(a, pp, zq, dphi, dz))
    N = len(corners_list)

    obs_dsv = make_dsv(dsv, args.ndsv)
    M = obs_dsv.shape[0]
    B_target = Gx * obs_dsv[:, 0]            # Hz units (scale-free solve)

    def entry(i, j):
        return _loop_Hz(obs_dsv[i], corners_list[j])

    modes = min(M, N) if args.k <= 0 else min(args.k, M, N)
    res = aca_tsvd(M, N, entry, modes=modes, kmax=min(M, N),
                   aca_eps=1.0e-8, method=3)
    k_use = res.modes if args.k <= 0 else min(args.k, res.modes)
    psi = pseudo_inverse_solve(res, B_target, k_mode=k_use)
    psi_zphi = psi.reshape(args.nz, args.nphi)

    print(f"[1] SF design: surface {args.nphi}x{args.nz} -> N={N} basis loops,"
          f" M={M} DSV obs (radius={dsv*1e3:.0f}mm)")
    print(f"    ACA+ k_aca = {res.k_aca} (of min(M,N)={min(M, N)}),"
          f" TSVD modes used = {k_use}")

    # ------- [2] equal-current contours of psi -------
    polylines, dI = contour_polylines_phi_z(
        psi_zphi, phi_grid, z_grid, args.nlevels)
    n_wire = len(polylines)
    if n_wire == 0:
        print("[2] no usable contours (psi is constant). aborting.")
        return
    print(f"[2] contours: {n_wire} fingerprint polylines,"
          f" dI = {dI:.4g} per wire (Hz units)")

    # ------- [3] single-stroke chain via cylinder-surface geodesics -------
    def build_chain_from_polylines(_polys):
        if args.chain_method == "field_aware":
            return single_stroke_field_aware_phi_z(_polys, a, n_blend=args.narc)
        elif args.chain_method == "nn_blend":
            return single_stroke_nn_blend_phi_z(_polys, a, n_blend=args.narc)
        elif args.chain_method == "kuijpers":
            return single_stroke_kuijpers_chain_phi_z(_polys, a)
        elif args.chain_method == "lobe":
            return single_stroke_lobe_chain_phi_z(_polys, a, n_arc=args.narc)
        return single_stroke_chain_phi_z(_polys, n_arc=args.narc)

    # ---- Path A: compensated iteration --------------------------------------
    # At each iteration, fit the chain's best single current I_w, compute
    # the residual r = B_target - I_w * Bz_chain_unit, and update
    # phi += step * pseudo_inverse(r).  The factorisation `res` is reused,
    # so each iteration costs only the back-substitution + the chain rebuild.
    # The chain's parasitic Bz at the DSV is absorbed into the new SF design.
    if args.compensated_iter > 0:
        print(f"[3a] Path-A compensated iteration "
              f"({args.compensated_iter} iters, step={args.compensated_step})")
        # Fixed-point iteration phi <- phi + step * pseudo_inverse(B_target -
        # I_w * Bz_chain_unit) is NOT a contraction on the discrete Gx chain
        # (the chain construction is too nonlinear in phi -- level-set
        # topology jumps at every iteration).  We keep the iteration here for
        # research purposes and track the BEST phi seen so we at least return
        # the best chain even if we don't converge.
        best_psi = psi.copy()
        best_polylines = polylines
        best_res_norm = float("inf")
        first_res_norm = None
        for it in range(args.compensated_iter):
            chain_it = build_chain_from_polylines(polylines)
            path_it = chain_phi_z_to_3d(chain_it, a)
            Bz_unit_it = bz_at(path_it, 1.0, obs_dsv)
            denom_it = float(np.dot(Bz_unit_it, Bz_unit_it))
            I_w_it = float(np.dot(Bz_unit_it, B_target) / denom_it) \
                if denom_it > 0 else 0.0
            chain_field_it = I_w_it * Bz_unit_it
            residual = B_target - chain_field_it
            res_norm = float(np.linalg.norm(residual)
                             / (np.linalg.norm(B_target) + 1e-30))
            if first_res_norm is None:
                first_res_norm = res_norm
            tag = ""
            if res_norm < best_res_norm:
                best_res_norm = res_norm
                best_psi = psi.copy()
                best_polylines = polylines
                tag = " <-- best"
            print(f"     iter {it+1}: I_w = {I_w_it:.3e},"
                  f" |B_target - chain|/|B_target| = {res_norm:.3e}{tag}")
            delta_phi = pseudo_inverse_solve(res, residual, k_mode=k_use)
            psi = psi + args.compensated_step * delta_phi
            psi_zphi = psi.reshape(args.nz, args.nphi)
            polylines, dI = contour_polylines_phi_z(
                psi_zphi, phi_grid, z_grid, args.nlevels)
            if not polylines:
                print("     iter aborted: psi became flat")
                break
        # Restore the best psi seen for the downstream chain build + eval.
        psi = best_psi
        polylines = best_polylines
        psi_zphi = psi.reshape(args.nz, args.nphi)
        print(f"     final: best residual = {best_res_norm:.3e}"
              f" (iter-1 baseline was {first_res_norm:.3e},"
              f" method={args.chain_method})")

    chain_phi_z = build_chain_from_polylines(polylines)
    path = chain_phi_z_to_3d(chain_phi_z, a)
    seg_lens = np.linalg.norm(np.diff(path, axis=0), axis=1)
    length = float(seg_lens.sum())
    max_rung, total_rung, n_rung, azim_rung = connector_lengths_phi_z(
        polylines, chain_phi_z, a)
    print(f"[3] single-stroke chain ({args.chain_method}"
          f"{' +compensated' if args.compensated_iter > 0 else ''}):"
          f" {len(path)} points,"
          f" {len(path) - 1} segments, wire length = {length:.3f} m"
          f" (one continuous conductor)")
    print(f"    connectors (inter-contour rungs): {n_rung} rungs,"
          f" longest = {max_rung*1e3:.1f} mm, total = {total_rung*1e3:.1f} mm;"
          f" AZIMUTHAL total = {azim_rung*1e3:.0f} mm"
          f" (<- field-impact correlator)")

    # ------- [6][7] field of single-stroke coil over the DSV and on x-axis ----
    # Fit one global current that minimises ||I*Bz_unit - Gx*x||_DSV; this is
    # the analogue of the "gain g" in demo_coil_design_gx.py and absorbs the
    # SF solve being scale-free (Hz vs T).
    Bz_unit_dsv = bz_at(path, 1.0, obs_dsv)
    denom = float(np.dot(Bz_unit_dsv, Bz_unit_dsv))
    I_w = float(np.dot(Bz_unit_dsv, B_target) / denom) if denom > 0 else 0.0
    Bz_dsv = I_w * Bz_unit_dsv
    rms = float(np.linalg.norm(Bz_dsv - B_target)
                / (np.linalg.norm(B_target) + 1e-30))
    print(f"[6] best-fit single current I_w = {I_w:.4g} A,"
          f" full-DSV RMS = {rms:.3e}  (Bz vs Gx*x)")

    xv = np.linspace(-dsv, dsv, 41)
    obs_xaxis = np.column_stack([xv, np.zeros_like(xv), np.zeros_like(xv)])
    Bz_xaxis = bz_at(path, I_w, obs_xaxis)
    G_fit = float(np.dot(xv, Bz_xaxis) / np.dot(xv, xv))
    nonlin = float(np.max(np.abs(Bz_xaxis - G_fit * xv))
                   / (np.max(np.abs(Bz_xaxis)) + 1e-30))
    print(f"[7] on x-axis: fitted dBz/dx = {G_fit:.4g},"
          f" nonlinearity over DSV = {nonlin:.3e}")
    print(f"    (target Gx = {Gx:.4g}; design Bz = Gx*x in the DSV)")

    if args.with_peec:
        run_peec_chain(path, I_w)

    # ------- plot -------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(12, 4))
        ax = fig.add_subplot(1, 3, 1, projection="3d")
        ax.plot(path[:, 0] * 1e3, path[:, 1] * 1e3, path[:, 2] * 1e3, lw=0.5)
        ax.set_title("Single-stroke Gx coil")
        ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]"); ax.set_zlabel("z [mm]")

        ax2 = fig.add_subplot(1, 3, 2)
        ax2.plot(xv * 1e3, Bz_xaxis, "o-", ms=3, label="single-stroke coil Bz")
        ax2.plot(xv * 1e3, G_fit * xv, "k--", label="linear fit")
        ax2.set_xlabel("x [mm]"); ax2.set_ylabel("Bz [T]")
        ax2.set_title(f"On x-axis (nonlin {nonlin:.1e})")
        ax2.legend(); ax2.grid(alpha=0.3)

        ax3 = fig.add_subplot(1, 3, 3)
        ax3.scatter(obs_dsv[:, 0] * 1e3, Bz_dsv, s=12, alpha=0.6,
                    label="coil Bz over DSV")
        xs = np.linspace(obs_dsv[:, 0].min(), obs_dsv[:, 0].max(), 50)
        ax3.plot(xs * 1e3, Gx * xs, "k--", label="target Gx*x")
        ax3.set_xlabel("x [mm]"); ax3.set_ylabel("Bz [T]")
        ax3.set_title(f"Full DSV (RMS {rms:.1e})")
        ax3.legend(); ax3.grid(alpha=0.3)

        out = HERE / "demo_sf_to_peec_gx.png"
        fig.tight_layout()
        fig.savefig(out, dpi=110)
        print(f"Saved plot: {out.name}")
    except ImportError:
        print("(matplotlib not installed; skipped plot)")


def run_peec_chain(path, I_w):
    """[4] CAD STEP (Frenet-swept round wire along the chain) + [5] PEEC L, R."""
    # Aggressively subsample (PEEC NxN mutual inductance scales poorly,
    # and the fingerprint chain has many nearly-parallel contour runs).
    target_n = 250
    stride = max(1, len(path) // target_n)
    cpath = path[::stride]
    if not np.allclose(cpath[-1], path[-1]):
        cpath = np.vstack([cpath, path[-1]])
    # Drop residual short segments: PEEC self-impedance diverges and OCC
    # Spline produces NaN tangents at near-zero-length steps.
    diffs = np.linalg.norm(np.diff(cpath, axis=0), axis=1)
    keep = np.concatenate([[True], diffs > 1.0e-3])   # 1 mm tol
    cpath = cpath[keep]

    step_path = str(HERE / "sf_coil_gx.step")
    try:
        import build123d as b3d
        from radia.coil_from_cad import export_step_with_labels
        mm = 1000.0
        pmm = cpath * mm
        seg_lens = np.linalg.norm(np.diff(pmm, axis=0), axis=1)
        r_wire = float(max(0.2, min(1.0, 0.25 * np.median(seg_lens))))

        # Multi-piece loft chain.  A SINGLE loft of >~50 cross-sections defeats
        # OCC (StdFail_NotDone) because the cumulative twist of the saddle
        # pattern overwhelms BRepOffsetAPI_ThruSections.  We split the chain
        # into short pieces (~20 sections each) and loft each piece, sharing
        # the boundary section so pieces visually connect; the result is a
        # Compound of solid "spools" that together cover the chain.
        # x_dir is propagated by PARALLEL TRANSPORT along the chain (each new
        # section's x_dir is the previous section's x_dir rotated minimally
        # to be perpendicular to the new tangent).  This is the twist-free
        # frame: it eliminates the negative volumes seen with both the Frenet
        # frame (auto from sweep) and the cylinder-radial frame (which winds
        # 360 deg per polyline traversal and gives the loft too much twist).
        n_sections_total = min(400, max(120, len(pmm) // 4))
        idx = np.unique(np.linspace(0, len(pmm) - 1, n_sections_total).astype(int))

        # Compute tangents and parallel-transported x_dir frame at each idx.
        tans = np.zeros((len(idx), 3))
        for k, i in enumerate(idx):
            if i == 0:
                d = pmm[1] - pmm[0]
            elif i == len(pmm) - 1:
                d = pmm[-1] - pmm[-2]
            else:
                d = pmm[i + 1] - pmm[i - 1]
            nrm = np.linalg.norm(d)
            tans[k] = d / nrm if nrm > 1e-12 else np.array([1.0, 0.0, 0.0])

        # Seed x_dir orthogonal to tan[0] using world-z projection
        world_z = np.array([0.0, 0.0, 1.0])
        seed = world_z - np.dot(world_z, tans[0]) * tans[0]
        if np.linalg.norm(seed) < 1e-3:
            seed = np.array([1.0, 0.0, 0.0]) - np.dot([1.0, 0.0, 0.0], tans[0]) * tans[0]
        xds = np.zeros((len(idx), 3))
        xds[0] = seed / np.linalg.norm(seed)
        for k in range(1, len(idx)):
            t0v, t1v = tans[k - 1], tans[k]
            ax = np.cross(t0v, t1v)
            s = np.linalg.norm(ax)
            if s < 1.0e-9:
                xds[k] = xds[k - 1]
            else:
                ax = ax / s
                c = float(np.clip(np.dot(t0v, t1v), -1.0, 1.0))
                ang = np.arccos(c)
                # Rodrigues rotation of xds[k-1] about ax by ang
                v = xds[k - 1]
                xds[k] = (v * c + np.cross(ax, v) * np.sin(ang)
                          + ax * np.dot(ax, v) * (1.0 - c))
            # Re-orthogonalise against current tangent (parallel transport drift)
            xds[k] = xds[k] - np.dot(xds[k], t1v) * t1v
            nn = np.linalg.norm(xds[k])
            xds[k] = xds[k] / nn if nn > 1e-9 else seed

        def make_section(k):
            i = int(idx[k])
            pl = b3d.Plane(origin=b3d.Vector(*pmm[i]),
                           x_dir=b3d.Vector(*xds[k]),
                           z_dir=b3d.Vector(*tans[k]))
            return pl.location * b3d.Circle(radius=r_wire).face()

        sections_per_piece = 20
        sec_cache = {}
        t0 = time.time()
        pieces = []
        n_fail = 0
        for j in range(0, len(idx) - 1, sections_per_piece - 1):
            piece_ks = list(range(j, min(j + sections_per_piece, len(idx))))
            if len(piece_ks) < 2:
                break
            piece_secs = []
            for k in piece_ks:
                if k not in sec_cache:
                    sec_cache[k] = make_section(k)
                piece_secs.append(sec_cache[k])
            try:
                piece_solid = b3d.loft(piece_secs, ruled=False)
                pieces.append(piece_solid)
            except Exception:
                n_fail += 1

        try:
            if not pieces:
                raise RuntimeError("all loft pieces failed")
            compound = b3d.Compound(pieces)
            try:
                compound.label = "coil"
            except Exception:
                pass
            export_step_with_labels([compound], step_path)
            total_vol = float(sum(p.volume for p in pieces))
            print(f"[4] CAD: round wire r={r_wire:.2f}mm loft-chained through "
                  f"{len(idx)} sections in {len(pieces)} pieces "
                  f"({n_fail} skipped) -> {os.path.basename(step_path)} "
                  f"({os.path.getsize(step_path)//1024} KB, vol={total_vol:.0f} mm^3,"
                  f" {time.time()-t0:.1f}s)")
        except Exception as se:
            poly = b3d.Polyline(*[b3d.Vector(float(x), float(y), float(z)) for x, y, z in pmm])
            b3d.export_step(poly, step_path)
            print(f"[4] CAD: loft chain failed ({type(se).__name__}); exported "
                  f"centerline Polyline WIRE -> {os.path.basename(step_path)} "
                  f"({os.path.getsize(step_path)//1024} KB)")
    except Exception as e:
        print(f"[4] CAD STEP export FAILED: {type(e).__name__}: {str(e)[:120]}")

    try:
        from radia.peec_matrices import PEECBuilder
        from radia.peec_topology import PEECCircuitSolver
        b = PEECBuilder()
        nodes = [b.add_node_at(float(p[0]), float(p[1]), float(p[2])) for p in cpath]
        # Wire cross-section: small enough that adjacent contour wires (at
        # neighbouring psi levels on the cylinder surface) don't overlap.
        # Estimate the local inter-wire spacing from segment lengths.
        seg_lens = np.linalg.norm(np.diff(cpath, axis=0), axis=1)
        w = hh = float(max(2.0e-4, min(2.0e-3, 0.2 * np.median(seg_lens))))
        for i in range(len(nodes) - 1):
            b.add_connected_segment(nodes[i], nodes[i + 1], w, hh, sigma=5.8e7)
        b.add_port(nodes[0], nodes[-1])
        topo = b.build_topology()
        solver = PEECCircuitSolver(topo)
        freq = 1.0e3
        t0 = time.time()
        Z = solver.compute_port_impedance(freq)
        Zc = complex(np.asarray(Z).ravel()[0]) if np.ndim(Z) else complex(Z)
        L = Zc.imag / (2 * np.pi * freq)
        R = Zc.real
        length = float(np.sum(np.linalg.norm(np.diff(cpath, axis=0), axis=1)))
        print(f"[5] PEEC ({len(nodes) - 1} filament segs, {length:.2f} m @ {freq:.0f} Hz, "
              f"{time.time() - t0:.1f}s): L = {L * 1e6:.3f} uH, R = {R * 1e3:.3f} mOhm")
        print(f"    (port = chain start to chain end; the single I_w = {I_w:.3g} A "
              f"realises the Gx*x field above)")
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"[5] PEEC FAILED: {type(e).__name__}: {str(e)[:120]}")


if __name__ == "__main__":
    main()
