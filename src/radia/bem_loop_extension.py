"""Genus-1 loop extension of the scalar P1 BIE with passive SIBC.

A single-valued surface scalar potential cannot carry the net circulating
current through the conductor handle. A cut-open unit-carrier potential g
adds that mode: phi_total = phi_u + alpha*g. Its jump and the weak normal
trace q_g = projection(-H_carrier.n) use the same oriented closed surface.

The ordinary scalar BIE rows and mean-zero phi_u gauge are unchanged.
The loop equation tests electric work against the distributed j_g=n x
(-grad_s g). Magnetic work uses the same distributed current PLUS the
normal-flux term of the exterior Calderon energy:

    Q_g,u = j_g^T S0 C - q_g^T (.5 M-DL)
    Q_g,g = j_g^T S0 j_g + q_g^T SL q_g.

S0 is the singular Galerkin scalar single layer for constant Cartesian
face currents; C maps closed P1 coefficients to those currents. The
kernel normalization is 1/(4*pi*distance), as in SL. The loop row is
K_Z(g,phi_total)+i*omega*mu0*Q(g,phi_scat)=-i*omega*<A_inc,j_g>.
Magnetic work acts on the scattered field; electric work and heat use
the total field. Hence the total-unknown source includes Q_g,u*phi_inc.
A current-only magnetic row is insufficient for finite impedance.

This mixed BIE/work system is not a symmetric matrix. Electric and
magnetic WORK pairings are reciprocal; a changed lift has a discrete
Calderon-identity defect that must converge with surface refinement.
Carrier, seam, singular quadrature and thin-ring approximation errors
remain separate validation duties. No finite-element accuracy claim.

Supports one flux-linked handle, dense undeformed flat triangular P1
surfaces, passive scalar impedance, and positive frequency. Automatic
carrier search proves enclosure/clearance; an explicit carrier must lie
inside the material wall. The alpha=0 frozen diagnostic reproduces the
ordinary scalar solver with its original gauge.

References
==========
* K. Sugahara, "Investigation of a Boundary Integral Equation
  n x H = J_s on Torus-Shaped Perfect Conductors," IEEE Trans.
  Antennas Propag., vol. 56, no. 3, pp. 722-727, Mar. 2008.  The dual
  defect of the same H^1(S) != 0 topology: on a torus the n x H = J_s
  BIE admits a spurious solution violating B . n = 0 (a harmonic
  null-space mode), closed there by ONE virtual-magnetic-current DOF +
  a one-point B . n = 0 constraint -- here the single-valued scalar
  potential LACKS the harmonic mode (representation deficit) and is
  closed by ONE loop DOF + the Faraday constraint.  Per handle, both
  formulations need exactly one extra DOF and one extra condition.
* M. Schoebinger and K. Hollaus, "An Effective Interface Approach for
  Multiply Connected Electromagnetic Shields," IEEE CEFC 2026
  conference digest, Thessaloniki, June 2026.  The FEM thin-shell
  sibling of the same topology treatment: the shielding sheet is
  reduced to a 2-D interface (1-D through-thickness analytic solution
  + nonlinear mu_eff lookup) and each hole gets ONE extra cohomology
  jump unknown (constant ``T = s_i e_z`` in the i-th hole), replacing
  the earlier non-physical auxiliary conductivity inside the holes.
* P. Dlotko, B. Kapidani, S. Pitassi, and R. Specogna, "Fake
  Conductivity or Cohomology: Which to Use When Solving Eddy Current
  Problems With h-Formulations?", IEEE Trans. Magn., vol. 55, no. 6,
  pp. 1-4, 2019.  The case AGAINST auxiliary-conductivity workarounds
  and FOR the cohomology treatment this module implements (via
  ``radia.cohomology``).

"""
from __future__ import annotations

import math
from collections import defaultdict, deque

import numpy as np

MU_0 = 4e-7 * math.pi


def A_from_filaments(points, filament_paths, currents):
    """Vector potential of straight filament segments (exact per segment).

    ``A_seg(x) = (mu0 I / 4 pi) t_hat * ln((L - xi + R2) / (-xi + R1))``
    with ``xi = (x - p1) . t_hat``, ``R1 = |x - p1|``, ``R2 = |x - p2|``.
    Supports complex per-filament currents.  Used as the ``A_inc_fn`` of
    ``solve_loop_extended`` for the PEEC (filament) coil source.

    Args:
        points: (n, 3) observation points [m].
        filament_paths: list of K filaments, each a list of (p1, p2)
            segment tuples (the ``coil_data["paths"]`` convention).
        currents: length-K complex per-filament currents [A].

    Returns:
        (n, 3) complex vector potential [T m].
    """
    pts = np.asarray(points, dtype=float)
    n = len(pts)
    A = np.zeros((n, 3), dtype=complex)
    for fil, Ik in zip(filament_paths, currents):
        Ik = complex(Ik)
        if Ik == 0:
            continue
        for (p1, p2) in fil:
            p1 = np.asarray(p1, dtype=float)
            p2 = np.asarray(p2, dtype=float)
            dl = p2 - p1
            L = float(np.linalg.norm(dl))
            if L < 1e-15:
                continue
            t_hat = dl / L
            xi = (pts - p1[None, :]) @ t_hat
            R1 = np.linalg.norm(pts - p1[None, :], axis=1)
            R2 = np.linalg.norm(pts - p2[None, :], axis=1)
            num = np.maximum(L - xi + R2, 1e-30)
            den = np.maximum(-xi + R1, 1e-30)
            val = (MU_0 * Ik / (4.0 * math.pi)) * np.log(num / den)
            A += val[:, None] * t_hat[None, :]
    return A


# ----------------------------------------------------------------------
# topology
# ----------------------------------------------------------------------
def cut_open_surface(pts, tris, cut_loop):
    """Duplicate the cut-loop vertices and reattach the R-side vertex
    fans (per-vertex side resolution -- robust to zig-zag cuts).
    Requires a CONSISTENTLY oriented ``tris``.  Returns
    ``(pts_open, tris_open, dup)`` with ``dup[orig] = duplicate id``."""
    pts = np.asarray(pts, dtype=float)
    tris = np.asarray(tris, dtype=np.int64)
    nv = len(pts)
    n_c = len(cut_loop)
    ced = [(cut_loop[i], cut_loop[(i + 1) % n_c]) for i in range(n_c)]
    ces = {(min(a, b), max(a, b)) for a, b in ced}
    edge_tris = defaultdict(list)
    tov = defaultdict(list)
    for ti, t in enumerate(tris):
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            edge_tris[(min(a, b), max(a, b))].append(ti)
        for v in t:
            tov[v].append(ti)

    def cdir(t, a, b):
        return any(t[k] == a and t[(k + 1) % 3] == b for k in range(3))

    spv = {}
    for (a, b) in ced:
        for ti in edge_tris[(min(a, b), max(a, b))]:
            s = 'L' if cdir(tris[ti], a, b) else 'R'
            spv[(ti, a)] = s
            spv[(ti, b)] = s
    for v in cut_loop:
        ring = tov[v]
        changed = True
        while changed:
            changed = False
            for ti in ring:
                if (ti, v) in spv:
                    continue
                for tj in ring:
                    if (tj, v) not in spv:
                        continue
                    sh = set(map(int, tris[ti])) & set(map(int, tris[tj]))
                    if len(sh) == 2 and v in sh and \
                            tuple(sorted(sh)) not in ces:
                        spv[(ti, v)] = spv[(tj, v)]
                        changed = True
                        break
        missing = [ti for ti in ring if (ti, v) not in spv]
        if missing:
            raise ValueError(
                f"cut_open_surface: fan side propagation incomplete at "
                f"vertex {v} (tris {missing}) -- is the cut a simple "
                f"closed loop on a manifold surface?")

    dup = {v: nv + i for i, v in enumerate(cut_loop)}
    pts_o = np.vstack([pts, pts[cut_loop]])
    tris_o = tris.copy()
    for (ti, v), s in spv.items():
        if s == 'R':
            for k in range(3):
                if tris_o[ti, k] == v:
                    tris_o[ti, k] = dup[v]
    return pts_o, tris_o, dup


# ----------------------------------------------------------------------
# Theta carrier
# ----------------------------------------------------------------------
def _midwall_ring(pts, pts_o, tris_o, cut_loop, vnorm, n_smooth=5):
    """Ring source inside the wall: ray-cast from each cut vertex along
    the inward normal to the opposite wall, take the mid point, then
    lightly smooth.  Same homology class as the cut by construction."""
    v0 = pts_o[tris_o[:, 0]]
    e1 = pts_o[tris_o[:, 1]] - v0
    e2 = pts_o[tris_o[:, 2]] - v0

    def ray_second_hit(p, d, t_min=2e-4):
        h = np.cross(d[None, :], e2)
        a = np.einsum('ij,ij->i', e1, h)
        ok = np.abs(a) > 1e-14
        f = np.where(ok, 1.0 / np.where(ok, a, 1.0), 0.0)
        s = p[None, :] - v0
        u_b = f * np.einsum('ij,ij->i', s, h)
        q = np.cross(s, e1)
        v_b = f * np.einsum('j,ij->i', d, q)
        t = f * np.einsum('ij,ij->i', e2, q)
        hit = ok & (u_b >= -1e-9) & (v_b >= -1e-9) \
            & (u_b + v_b <= 1 + 1e-9) & (t > t_min)
        return float(np.min(t[hit])) if hit.any() else None

    ring = []
    for v in cut_loop:
        t_hit = ray_second_hit(pts[v], -vnorm[v])
        if t_hit is None:
            raise ValueError(
                f"loop extension: no opposite wall found from cut vertex "
                f"{v} along the inward normal -- cannot place the Theta "
                f"ring source inside the material.")
        ring.append(pts[v] - 0.5 * t_hit * vnorm[v])
    ring = np.array(ring)
    for _ in range(n_smooth):
        ring = 0.5 * ring + 0.25 * (np.roll(ring, 1, axis=0)
                                    + np.roll(ring, -1, axis=0))
    return ring


def _ring_H_factory(ring_pts):
    rp1 = np.asarray(ring_pts, dtype=float)
    rp2 = np.roll(rp1, -1, axis=0)

    def H_ring(robs):
        robs = np.asarray(robs, dtype=float)
        out = np.zeros((len(robs), 3))
        for s0 in range(0, len(robs), 400):
            s1 = min(len(robs), s0 + 400)
            a = rp1[None, :, :] - robs[s0:s1, None, :]
            b = rp2[None, :, :] - robs[s0:s1, None, :]
            na = np.linalg.norm(a, axis=2)
            nb = np.linalg.norm(b, axis=2)
            cr = np.cross(a, b)
            d2 = np.einsum('smd,smd->sm', cr, cr)
            ab = np.einsum('smd,smd->sm', a, b)
            coef = np.where(d2 > 1e-30,
                            (na + nb) * (1 - ab / (na * nb))
                            / np.maximum(d2, 1e-30), 0.0)
            out[s0:s1] = np.einsum('sm,smd->sd', coef, cr) / (4 * math.pi)
        return out

    return H_ring


def _project_ring_neumann(pts, tris, areas, normals, H_ring, M_inv, *, quadrature_n=6):
    """L2-project -H_ring.n_triangle into the closed surface P1 space.

    Tensor Gauss on the Duffy triangle has normalized weights summing to
    one. Batch panels so ring-field evaluation does not allocate all
    panel/quadrature/carrier pairs at once. ``quadrature_n`` is exposed
    here for convergence diagnostics; the solve uses six points per axis.
    """
    x, w = np.polynomial.legendre.leggauss(quadrature_n)
    x, w = (x + 1) / 2, w / 2
    u, v = np.meshgrid(x, x, indexing="ij")
    wu, wv = np.meshgrid(w, w, indexing="ij")
    bary = np.stack([1-u, u*(1-v), u*v], axis=-1).reshape(-1, 3)
    weights = (2*u*wu*wv).ravel()
    rhs = np.zeros(len(pts))
    for first in range(0, len(tris), 100):
        ts = tris[first:first+100]
        gp = np.einsum("qk,tkd->tqd", bary, pts[ts])
        field = H_ring(gp.reshape(-1, 3)).reshape(gp.shape)
        values = -np.einsum("tqd,td->tq", field, normals[first:first+100])
        local = areas[first:first+100, None] * (values*weights) @ bary
        np.add.at(rhs, ts.ravel(), local.ravel())
    return M_inv @ rhs


def _theta_by_path_integration(pts_o, tris_o, dup, cut_loop, H_ring,
                               jump_tol=5e-3):
    """Theta on the open mesh: spanning-tree path integration of
    -H_ring . dl (4-pt Gauss per edge).  Verifies the uniform +-1 jump
    across the cut and raises otherwise."""
    nv_o = len(pts_o)
    adj = defaultdict(list)
    for t in tris_o:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            adj[a].append(b)
            adj[b].append(a)
    Theta = np.full(nv_o, np.nan)
    root = int(tris_o[0, 0])
    Theta[root] = 0.0
    par = {root: None}
    dq = deque([root])
    order = [root]
    while dq:
        u = dq.popleft()
        for w in adj[u]:
            if w not in par:
                par[w] = u
                order.append(w)
                dq.append(w)
    gl_x, gl_w = np.polynomial.legendre.leggauss(4)
    gl_t = 0.5 * (gl_x + 1)
    gl_w = 0.5 * gl_w
    for w in order[1:]:
        u = par[w]
        p1, p2 = pts_o[u], pts_o[w]
        gp = p1[None, :] + gl_t[:, None] * (p2 - p1)[None, :]
        Hs = H_ring(gp)
        Theta[w] = Theta[u] - float(np.sum(gl_w * (Hs @ (p2 - p1))))
    if np.isnan(Theta).any():
        raise ValueError("loop extension: open mesh is disconnected -- "
                         "Theta path integration could not reach all "
                         "vertices.")
    jumps = np.array([Theta[dup[v]] - Theta[v] for v in cut_loop])
    if jumps.std() > jump_tol or abs(abs(jumps.mean()) - 1.0) > jump_tol:
        raise ValueError(
            f"loop extension: Theta jump across the cut is "
            f"{jumps.mean():+.4f} +- {jumps.std():.1e} (expected uniform "
            f"+-1).  The mid-wall ring source construction failed for "
            f"this geometry (it must stay inside the material wall).")
    return Theta, float(jumps.mean())


def _segment_surface_distance(a, b, triangles):
    """Minimum Euclidean distance from a closed segment to flat triangles.

    Includes face piercing, coplanar overlap, and edge/vertex contact. No
    finite sampling of the carrier is used as a containment certificate.
    """
    p, q, r = triangles[:, 0], triangles[:, 1], triangles[:, 2]
    e, f, d = q-p, r-p, b-a
    n = np.cross(e, f)
    nn = np.einsum('ij,ij->i', n, n)
    ee = np.einsum('ij,ij->i', e, e)
    ef = np.einsum('ij,ij->i', e, f)
    ff = np.einsum('ij,ij->i', f, f)

    def face_distance(x):
        w = x-p
        en = np.einsum('ij,ij->i', w, e)
        fn = np.einsum('ij,ij->i', w, f)
        u, v = (ff*en-ef*fn)/nn, (ee*fn-ef*en)/nn
        inside = (u >= 0) & (v >= 0) & (u+v <= 1)
        h = np.einsum('ij,ij->i', w, n)
        return np.where(inside, h*h/nn, np.inf)

    best = np.minimum(face_distance(a), face_distance(b))
    den = n @ d
    t = np.divide(np.einsum('ij,ij->i', p-a, n), den,
                  out=np.full(len(p), np.inf), where=den != 0)
    valid = (t >= 0) & (t <= 1)
    x = a + np.where(valid, t, 0)[:, None]*d
    # Barycentric inclusion at the plane crossing.
    w = x-p
    en, fn = np.einsum('ij,ij->i', w, e), np.einsum('ij,ij->i', w, f)
    u, v = (ff*en-ef*fn)/nn, (ee*fn-ef*en)/nn
    if np.any(valid & (u >= 0) & (v >= 0) & (u+v <= 1)):
        return 0.0
    dd = float(d @ d)
    for c, end in ((p,q), (q,r), (r,p)):
        g, w = end-c, a-c
        gg = np.einsum('ij,ij->i', g, g)
        dg, dw = g @ d, w @ d
        gw = np.einsum('ij,ij->i', g, w)
        det = dd*gg-dg*dg
        ss = np.divide(dg*gw-gg*dw, det, out=np.zeros_like(det),
                       where=det > dd*gg*1e-14)
        tt = (gw+ss*dg)/gg
        interior = (det > dd*gg*1e-14) & (ss >= 0) & (ss <= 1) & (tt >= 0) & (tt <= 1)
        delta = w+ss[:,None]*d-tt[:,None]*g
        best = np.minimum(best, np.where(interior, np.einsum('ij,ij->i', delta, delta), np.inf))
        for endpoint in (a,b):
            t_edge = np.clip(np.einsum('ij,ij->i', endpoint-c, g)/gg, 0, 1)
            delta = endpoint-c-t_edge[:,None]*g
            best = np.minimum(best, np.einsum('ij,ij->i', delta, delta))
        for endpoint in (c,end):
            t_seg = np.clip((endpoint-a) @ d/dd, 0, 1)
            delta = endpoint-a-t_seg[:,None]*d
            best = np.minimum(best, np.einsum('ij,ij->i', delta, delta))
    return float(np.sqrt(max(0., best.min())))


def _surface_winding(point, triangles):
    """Oriented solid-angle winding; caller supplies a closed oriented mesh."""
    a, b, c = (triangles[:,i]-point for i in range(3))
    la, lb, lc = (np.linalg.norm(x, axis=1) for x in (a,b,c))
    numerator = np.einsum('ij,ij->i', a, np.cross(b,c))
    denominator = (la*lb*lc + np.einsum('ij,ij->i', a,b)*lc
                   + np.einsum('ij,ij->i', b,c)*la + np.einsum('ij,ij->i', c,a)*lb)
    return float(np.sum(2*np.arctan2(numerator, denominator))/(4*math.pi))


def _auto_material_ring(pts, tris, *, return_diagnostics=False):
    """Find a circular carrier inside a closed, oriented, z-axis bore wall.

    Edge/vertex hits and coplanar or tangent rays are ambiguous, never
    interpreted as an odd crossing proving exterior. Retry at other
    azimuths and shifted heights. Winding at one point plus positive
    segment-to-surface clearance certifies the entire connected polygon
    inside the flat surface, including non-axisymmetric walls. This is
    deliberately limited to circular carriers surrounding the z axis.
    """
    import json
    pts, tris = np.asarray(pts, float), np.asarray(tris, int)
    if (pts.ndim != 2 or pts.shape[1] != 3 or len(pts) < 4
            or tris.ndim != 2 or tris.shape[1] != 3 or len(tris) < 4
            or not np.all(np.isfinite(pts))
            or np.any(tris < 0) or np.any(tris >= len(pts))):
        raise ValueError("Automatic carrier requires finite points and valid triangle indices")
    directed = np.concatenate([tris[:, [0,1]], tris[:, [1,2]], tris[:, [2,0]]])
    _, inverse, counts = np.unique(np.sort(directed, axis=1), axis=0,
                                   return_inverse=True, return_counts=True)
    orientation = np.bincount(inverse, weights=np.where(directed[:,0] < directed[:,1], 1, -1))
    if np.any(counts != 2) or np.any(orientation != 0):
        raise ValueError("Automatic carrier requires a closed consistently oriented surface")
    triangles = pts[tris]
    origin = triangles[:,0]
    edge1, edge2 = triangles[:,1]-origin, triangles[:,2]-origin
    scale = float(np.ptp(pts, axis=0).max())
    tol = max(scale*1e-9, 1e-14)
    cross = np.cross(edge1, edge2)
    if not np.all(np.isfinite(pts)) or np.any(np.linalg.norm(cross,axis=1) == 0):
        raise ValueError("Automatic carrier requires finite nondegenerate triangles")
    diagnostics = dict(candidate_count=0, ray_attempts=[], candidates=[],
                       minimum_wall_distance=None, clearance_tolerance=tol)

    def intervals(z, angle):
        direction = np.array([math.cos(angle), math.sin(angle), 0.])
        h = np.cross(direction, edge2)
        determinant = np.einsum('ij,ij->i', edge1, h)
        valid = np.abs(determinant) > np.linalg.norm(cross,axis=1)*1e-12
        inverse = np.divide(1., determinant, out=np.zeros_like(determinant), where=valid)
        offset = np.array([0., 0., z])-origin
        u = inverse*np.einsum('ij,ij->i', offset, h)
        q = np.cross(offset, edge1)
        v = inverse*(q @ direction)
        distance = inverse*np.einsum('ij,ij->i', edge2, q)
        hit = valid & (u >= -1e-10) & (v >= -1e-10) & (u+v <= 1+1e-10) & (distance > tol)
        boundary = hit & ((u <= 1e-10) | (v <= 1e-10) | (u+v >= 1-1e-10))
        coplanar = (~valid) & (np.abs(np.einsum('ij,ij->i', offset,cross)) <= tol*np.linalg.norm(cross,axis=1))
        hits = np.sort(distance[hit])
        ambiguous = bool(np.any(boundary) or np.any(coplanar) or np.any(np.diff(hits) <= tol))
        reason = ('ambiguous_tangent_edge_or_vertex' if ambiguous else
                  'odd_crossings_axis_in_material' if len(hits)%2 else
                  'no_crossings' if not len(hits) else 'ok')
        diagnostics['ray_attempts'].append(dict(z=float(z), angle=float(angle), crossings=len(hits), reason=reason))
        return hits.reshape(-1,2) if reason == 'ok' else np.empty((0,2))

    lo, hi = float(pts[:,2].min()), float(pts[:,2].max())
    dz = (hi-lo)/32
    candidates = []
    # Include slab midplanes so a thin axial shelf cannot fall between
    # all uniform samples. Horizontal triangle faces identify actual shelves.
    horizontal = np.linalg.norm(cross[:,:2],axis=1) <= np.linalg.norm(cross,axis=1)*1e-12
    levels = np.unique(np.r_[lo, hi, np.round(origin[horizontal,2]/tol)*tol])
    heights = np.unique(np.r_[.5*(levels[:-1]+levels[1:]),
        *[np.linspace(lo,hi,33)[1:-1]+shift*dz for shift in (0., .271, -.193)]])
    for z in heights[(heights > lo+tol) & (heights < hi-tol)]:
        spans = []
        common = None
        for angle in (.137, .731, 1.913, 2.719):
            intervals_at_angle = intervals(z,angle)
            spans.extend(intervals_at_angle.tolist())
            if common is None and len(intervals_at_angle):
                common = intervals_at_angle
        # Intersection across azimuths can find a thin common annulus
        # even when each individual interval midpoint is outside it.
        if common is None:
            common = np.empty((0,2))
        for angle in np.arange(8)*2*math.pi/8+.319:
            other = intervals(z,angle)
            if len(other) and len(common):
                common = np.array([(max(a,c),min(b,d)) for a,b in common
                                   for c,d in other if max(a,c)<min(b,d)]).reshape(-1,2)
        spans.extend(common.tolist())
        for inner, outer in spans:
            radius = .5*(inner+outer)
            score = min(.5*(outer-inner),z-lo,hi-z)
            candidates.append((score,radius,float(z)))
    candidates = sorted(set(candidates),reverse=True)
    diagnostics['candidate_count'] = len(candidates)
    theta = np.arange(256)*2*math.pi/256
    for score, radius, z in candidates:
        record = dict(radius=radius,z=z,minimum_wall_distance=None)
        diagnostics['candidates'].append(record)
        ring = np.c_[radius*np.cos(theta),radius*np.sin(theta),np.full(256,z)]
        if abs(_surface_winding(np.array([0.,0.,z]),triangles)) > 1e-6:
            record['reason'] = 'axis_not_in_bore'; continue
        if abs(abs(_surface_winding(ring[0],triangles))-1) > 1e-6:
            record['reason'] = 'carrier_start_outside'; continue
        clearance = math.inf
        for a,b in zip(ring,np.roll(ring,-1,axis=0)):
            clearance = min(clearance, _segment_surface_distance(a,b,triangles))
            if clearance == 0.:
                break
        record['minimum_wall_distance'] = clearance
        diagnostics['minimum_wall_distance'] = clearance
        if clearance <= tol:
            record['reason'] = 'surface_contact_or_crossing'; continue
        record['reason'] = 'accepted'
        result = ((float(radius),float(z)),ring)
        return (*result,diagnostics) if return_diagnostics else result
    # Unvisited candidates do not exist on failure; every rejection is retained.
    raise ValueError("Automatic workpiece cohomology could not find an interior carrier around the z axis. "
        "Supported geometry: one through-hole around z (tube, ring, stepped shaft). "
        "Supply a validated section_anchor and carrier_ring for other shapes; "
        "no uncorrected heating result is returned. Diagnostics: " + json.dumps(diagnostics))


# ----------------------------------------------------------------------
# main entry
# ----------------------------------------------------------------------
def solve_loop_extended(bem_solver, phi_inc_nodal, Z_s, omega, A_inc_fn, *, section_anchor=None, carrier_ring=None):
    """Solve the loop-extended scalar BIE + SIBC on a genus-1 workpiece.

    Args:
        bem_solver: a ``ScalarBIESIBCSolver`` built on the CLOSED
            workpiece surface mesh with ``assemble_dense=True`` and the
            intree P1 path (attributes ``M/K/SL/DL/M_inv/mesh`` used).
        phi_inc_nodal: (ndof,) complex incident scalar potential at the
            H1 nodes (e.g. the surface-Poisson reconstruction, or an
            exact expression for a uniform field).
        carrier_ring: optional ordered closed-loop vertices inside the
            material wall, in metres (closing vertex is implicit). Use when
            normal-ray construction fails on a stepped wall. The unit-jump
            check still applies. If both carrier_ring and section_anchor
            are omitted, an interior circular carrier is found automatically
            from the closed surface for a single through-hole around z.
        section_anchor: optional (rho, z) point strictly inside the material
            meridional section, in metres. Supply this for stepped profiles
            whose surface-node mean lies outside the material.
        Z_s: complex Leontovich surface impedance (global scalar).
        omega: angular frequency [rad/s].
        A_inc_fn: callable ``A_inc_fn(points (n,3)) -> (n,3) complex`` --
            incident vector potential, used for the linked-flux term of
            the Faraday closure on the cut loop.

    Returns dict with ``alpha`` (net circulating current, complex [A]),
    ``P_total`` [W], ``H_t_rms`` [A/m], ``phi_u``, ``phi_open``,
    ``theta_jump``, ``cut_n_vertices``, plus ``P_frozen`` / ``Ht_frozen``
    (the alpha=0 sub-solve, == the plain solver, for diagnostics).
    ``carrier_diagnostics`` records automatic search rejections and the
    accepted polygon's minimum wall clearance (None for a supplied carrier).
    """
    from ngsolve import BND

    if (not np.isfinite(omega) or omega <= 0 or np.ndim(Z_s) != 0
            or not np.isfinite(Z_s) or complex(Z_s).real < 0):
        raise ValueError("Loop SIBC requires positive frequency and finite passive scalar Z_s")
    mesh = bem_solver.mesh
    if mesh.GetCurveOrder() > 1 or mesh.deformation is not None:
        raise ValueError("Loop SIBC requires undeformed flat surface triangles")
    M, K = bem_solver.M, bem_solver.K
    SL, DL, M_inv = bem_solver.SL, bem_solver.DL, bem_solver.M_inv
    if SL is None or DL is None:
        raise ValueError("loop extension needs assemble_dense=True "
                         "(dense SL/DL) on the ScalarBIESIBCSolver.")
    gamma = Z_s / (1j * omega * MU_0) if omega > 0 else 0.0

    pts = np.array([[mesh.vertices[i].point[j] for j in range(3)]
                    for i in range(mesh.nv)])
    tris = np.array([[v.nr for v in el.vertices]
                     for el in mesh.Elements(BND)], dtype=np.int64)
    nv, nt = len(pts), len(tris)
    if bem_solver.ndof != nv:
        raise ValueError(
            f"loop extension supports the P1 nodal path only "
            f"(ndof={bem_solver.ndof} != nv={nv}).")

    # verify consistent orientation (the extractors now guarantee it; a
    # mesh from another route must be oriented BEFORE the solver is
    # built, because the BEM operators bake the winding in).
    dcount = defaultdict(int)
    for t in tris:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            dcount[(int(a), int(b))] += 1
    conflicts = sum(1 for c in dcount.values() if c > 1)
    if conflicts:
        raise ValueError(
            f"loop extension: surface winding is inconsistent "
            f"({conflicts} directed-edge conflicts).  Orient the mesh "
            f"(surface_mesh_extract.orient_surface_triangles) BEFORE "
            f"building the BEM solver.")

    # genus / cut
    E = {(min(a, b), max(a, b)) for t in tris for a, b in
         ((t[0], t[1]), (t[1], t[2]), (t[2], t[0]))}
    chi = nv - len(E) + nt
    genus = (2 - chi) // 2
    if genus != 1:
        raise ValueError(
            f"loop extension supports genus-1 surfaces (got genus="
            f"{genus}, chi={chi}).  genus 0 needs no extension; "
            f"genus >= 2 needs one DOF per flux-linked handle (not "
            f"implemented).")
    # Cut selection through the repo's single cohomology engine
    # (radia.cohomology, the gmsh-free port).  The flux-linked cut must
    # be the PURE toroidal class (w_tor, w_pol) = (+-1, 0): a mixed
    # representative like (1, -1) is a valid basis element but forces
    # alpha to carry an equal poloidal component (wrong physics) and its
    # zig-zag geometry breaks the mid-wall Theta ray-cast.  A fixed
    # 2-loop generator basis is NOT guaranteed class-pure, so classify
    # ALL cotree fundamental
    # cycles instead: the class map is a linear image of the harmonic
    # period vector, calibrated on one independent (QR-pivoted) pair, so
    # every cycle costs O(1) -- then take the geometrically shortest
    # pure-toroidal simple loop.  Winding numbers of a CLOSED loop are
    # exact integers: toroidal about the z axis, poloidal about the
    # wall-section centroid (rho0, z0) -- valid for any z-axis genus-1
    # workpiece (torus, tube, ring), whose section centroid lies inside
    # the material wall.
    import scipy.linalg as sla
    from radia.cohomology import surface_fundamental_cycles

    _b1, Pi, expand, _cotree = surface_fundamental_cycles(tris, nv=nv)
    automatic_carrier = section_anchor is None and carrier_ring is None
    carrier_diagnostics = None
    if automatic_carrier:
        section_anchor, carrier_ring, carrier_diagnostics = _auto_material_ring(
            pts, tris, return_diagnostics=True)
    rho_all = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)
    rho0, z0 = float(rho_all.mean()), float(pts[:, 2].mean())

    if section_anchor is not None:
        anchor = np.asarray(section_anchor, dtype=float)
        if anchor.shape != (2,) or not np.all(np.isfinite(anchor)) or anchor[0] <= 0:
            raise ValueError("section_anchor must be a finite (rho, z) pair with rho > 0, inside the material section")
        rho0, z0 = map(float, anchor)

    def _windings(p):
        q = pts[list(p) + [p[0]]]
        th = np.unwrap(np.arctan2(q[:, 1], q[:, 0]))
        rho = np.sqrt(q[:, 0] ** 2 + q[:, 1] ** 2)
        ph = np.unwrap(np.arctan2(q[:, 2] - z0, rho - rho0))
        return np.array([(th[-1] - th[0]) / (2 * math.pi),
                         (ph[-1] - ph[0]) / (2 * math.pi)])

    _, _, piv = sla.qr(Pi.T, pivoting=True, mode="economic")
    sel = piv[:2]
    W_sel = np.array([_windings(expand(k)) for k in sel])      # (2, 2)
    if abs(np.linalg.det(W_sel)) < 0.5:
        raise ValueError(
            f"loop extension: winding map is singular on the generator "
            f"pair (windings {np.rint(W_sel).astype(int).tolist()}) -- "
            f"the workpiece is not a z-axis genus-1 body in the assumed "
            f"sense (toroidal about z, poloidal about the wall section).")
    # windings(cycle e) = Wmap @ Pi[e]  (linear class map, calibrated once)
    Wmap = W_sel.T @ np.linalg.inv(Pi[sel].T)
    W_all = Pi @ Wmap.T                                        # (nc, 2)
    W_int = np.rint(W_all)
    cands = np.where((np.abs(W_all - W_int) < 0.2).all(axis=1)
                     & (np.abs(W_int[:, 0]) == 1)
                     & (W_int[:, 1] == 0))[0]
    if len(cands) == 0:
        raise ValueError(
            f"loop extension: no PURE toroidal fundamental cycle found "
            f"among {Pi.shape[0]} cotree candidates (generator windings "
            f"{np.rint(W_sel).astype(int).tolist()}).  Re-meshing "
            f"usually resolves this.")

    def _geo_len(p):
        q = pts[list(p) + [p[0]]]
        return float(np.sum(np.linalg.norm(np.diff(q, axis=0), axis=1)))

    cut_k = min((int(k) for k in cands),
                key=lambda k: _geo_len(expand(k)))
    cut = expand(cut_k)
    # A loop current coefficient changes sign when the cut orientation
    # changes.  Canonicalize the selected representative to positive
    # toroidal winding so alpha phase, and therefore the reported regime,
    # is independent of the spanning-tree orientation.
    if int(W_int[cut_k, 0]) < 0:
        cut = list(reversed(cut))
    n_cut = len(cut)

    pts_o, tris_o, dup = cut_open_surface(pts, tris, cut)
    nv_o = len(pts_o)
    Tmap = np.arange(nv_o)
    for v, vd in dup.items():
        Tmap[vd] = v

    # geometry tables
    areas = np.zeros(nt)
    normals = np.zeros((nt, 3))
    gvecs = np.zeros((nt, 3, 3))
    for ti, t in enumerate(tris_o):
        P = pts_o[t]
        nvec = np.cross(P[1] - P[0], P[2] - P[0])
        A2 = np.linalg.norm(nvec)
        n = nvec / A2
        areas[ti] = 0.5 * A2
        normals[ti] = n
        for k in range(3):
            opp = P[(k + 2) % 3] - P[(k + 1) % 3]
            h = np.cross(n, opp) / (2 * areas[ti])
            if h @ (P[k] - P[(k + 1) % 3]) < 0:
                h = -h
            gvecs[ti, k] = h
    cents = pts_o[tris_o].mean(axis=1)
    vnorm = np.zeros((nv, 3))
    for ti, t in enumerate(tris_o):
        for k in range(3):
            vnorm[Tmap[t[k]]] += normals[ti] * areas[ti]
    vnorm /= np.maximum(np.linalg.norm(vnorm, axis=1), 1e-30)[:, None]

    # Theta carrier
    ring = _midwall_ring(pts, pts_o, tris_o, cut, vnorm) if carrier_ring is None else np.asarray(carrier_ring, dtype=float)
    if ring.ndim != 2 or ring.shape[1] != 3 or len(ring) < 3 or not np.all(np.isfinite(ring)):
        raise ValueError("carrier_ring must contain at least three finite 3D points inside the material wall")
    H_ring = _ring_H_factory(ring)
    Theta, theta_jump = _theta_by_path_integration(
        pts_o, tris_o, dup, cut, H_ring)
    Theta = Theta.astype(complex)

    qT = _project_ring_neumann(pts, tris, areas, normals, H_ring, M_inv)
    rK_T = np.zeros(nv, dtype=complex)
    for ti, t in enumerate(tris_o):
        gTh = (gvecs[ti, 0] * Theta[t[0]] + gvecs[ti, 1] * Theta[t[1]]
               + gvecs[ti, 2] * Theta[t[2]])
        for k in range(3):
            rK_T[Tmap[t[k]]] += areas[ti] * (gvecs[ti, k] @ gTh)
    a_col = SL @ (gamma * (M_inv @ rK_T) - qT.astype(complex))

    A_sys = (0.5 * M - DL + gamma * (SL @ M_inv @ K)).astype(complex)
    RHS = (M @ np.asarray(phi_inc_nodal, dtype=complex))

    from .bem_loop_work import _loop_work_row
    fE, e_loop, fP, p_loop, loop_source, work_diagnostics = _loop_work_row(
        bem_solver, pts, tris, tris_o, Tmap, gvecs, normals, areas, Theta,
        Z_s, A_inc_fn, qT, phi_inc_nodal)
    f_alpha = e_loop + 1j * omega * p_loop

    # assemble + solve (Lagrange mean-zero gauge, production style)
    Mrow = M.sum(axis=1).astype(complex)
    N = nv + 2
    A2 = np.zeros((N, N), dtype=complex)
    b2 = np.zeros(N, dtype=complex)
    A2[:nv, :nv] = A_sys
    A2[:nv, nv] = a_col
    A2[:nv, nv + 1] = Mrow
    A2[nv, :nv] = fE + 1j * omega * fP
    A2[nv, nv] = f_alpha
    A2[nv + 1, :nv] = Mrow
    b2[:nv] = RHS
    b2[nv] = -1j * omega * loop_source
    u = np.linalg.solve(A2, b2)
    residual = A2 @ u - b2
    linear_residual_rel = float(np.linalg.norm(residual) / max(np.linalg.norm(b2), np.finfo(float).tiny))
    faraday_residual_rel = float(abs(residual[nv]) / max(abs(b2[nv]), abs((A2 @ u)[nv]), np.finfo(float).tiny))
    if not np.isfinite(linear_residual_rel) or linear_residual_rel > 1e-6 or faraday_residual_rel > 1e-6:
        raise RuntimeError("Loop SIBC true residual / Faraday closure failed")
    phi_u, alpha = u[:nv], complex(u[nv])

    def _P_Ht(phi_u_v, alpha_v):
        phi_o = phi_u_v[Tmap] + alpha_v * Theta
        s2 = 0.0
        for ti, t in enumerate(tris_o):
            g = -(gvecs[ti, 0] * phi_o[t[0]] + gvecs[ti, 1] * phi_o[t[1]]
                  + gvecs[ti, 2] * phi_o[t[2]])
            s2 += areas[ti] * float(np.sum(np.abs(g) ** 2))
        return (0.5 * Z_s.real * s2, math.sqrt(s2 / float(areas.sum())))

    P_total, H_t_rms = _P_Ht(phi_u, alpha)

    # Complete total field, including the multivalued carrier.  A plain
    # single-valued phi.B integral would omit the cut contribution.
    phi_open = phi_u[Tmap] + alpha * Theta
    H_total = -np.einsum('tik,ti->tk', gvecs, phi_open[tris_o])
    H_inc = -np.einsum('tik,ti->tk', gvecs,
                       np.asarray(phi_inc_nodal)[Tmap][tris_o])
    A_inc = np.asarray(A_inc_fn(cents), dtype=complex)
    reaction_magnetic = np.sum(areas * np.einsum(
        'ij,ij->i', np.cross(normals, H_total), A_inc))
    reaction_electric = Z_s / (1j * omega) * np.sum(
        areas * np.einsum('ij,ij->i', H_inc, H_total))

    # frozen sub-solve (diagnostic; == the plain production solve)
    Ng = nv + 1
    A0 = np.zeros((Ng, Ng), dtype=complex)
    b0 = np.zeros(Ng, dtype=complex)
    A0[:nv, :nv] = A_sys
    A0[:nv, nv] = Mrow
    A0[nv, :nv] = Mrow
    b0[:nv] = RHS
    phi_f = np.linalg.solve(A0, b0)[:nv]
    P_frozen, Ht_frozen = _P_Ht(phi_f, 0.0)

    return {
        "loop_work_diagnostics": work_diagnostics,
        "theta_jump_max_deviation": float(max(abs(Theta[dup[v]]-Theta[v]-theta_jump) for v in cut)),
        "linear_residual_rel": linear_residual_rel,
        "faraday_residual_rel": faraday_residual_rel,
        "section_anchor": (rho0, z0),
        "automatic_carrier": automatic_carrier,
        "carrier_diagnostics": carrier_diagnostics,
        "alpha": alpha,
        "P_total": float(P_total),
        "H_t_tri": H_total,
        "q_tri": 0.5 * Z_s.real * np.sum(np.abs(H_total)**2, axis=1),
        "reaction_integral": complex(reaction_magnetic + reaction_electric),
        "P_reaction": float(-0.5 * omega * (reaction_magnetic + reaction_electric).imag),
        "H_t_rms": float(H_t_rms),
        "P_frozen": float(P_frozen),
        "Ht_frozen": float(Ht_frozen),
        "phi_u": phi_u,
        "phi_open": phi_u[Tmap] + alpha * Theta,
        "theta_jump": float(theta_jump),
        "cut_n_vertices": int(n_cut),
        "genus": 1,
    }
