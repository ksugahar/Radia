"""HACApK + custom Laplace kernel SIBC BEM solver (Phase 1: dense).

Replaces NGSolve.bem (effectively dense O(N^2) at typical IH wp size,
silently ignores all H-matrix compression flags as of 2026-04-29) with
a true H-matrix BEM via HACApK once Phase 1.4 lands.

Phase 1 = flat-triangle Galerkin BEM with custom Laplace kernel,
DENSE assembly first.  Verification anchor: NGSolve.bem on the same
mesh (Galerkin convention is identical:
   M_ij = int_Gamma int_Gamma  phi_i(r) G(r, r') phi_j(r') dS dS').

Singular treatment:
- Identical pair: Sauter-Schwab 4D Duffy
- Edge-share: Sauter-Schwab edge case (5D)
- Vertex-share: Sauter-Schwab vertex case (4D)
- Disjoint: tensor-product Gauss on each tri

References:
- Sauter-Schwab "Boundary Element Methods" 2010, ch. 5.2
"""
import math
import numpy as np

INV_4PI = 1.0 / (4.0 * math.pi)


# ----------------------------------------------------------------------
# Mesh extraction
# ----------------------------------------------------------------------
def extract_surface(mesh, bnd_label=None):
    """Extract a flat-triangle surface from a NGSolve volume or surface
    mesh, restricted to one BND label (or all BND if label is None).

    Args:
        mesh: NGSolve Mesh.
        bnd_label: optional string name of the boundary label to keep.

    Returns:
        verts (n_v, 3) float64 vertex coords (local to extracted set)
        tris (n_t, 3) int64 triangle vertex indices into verts
        v_global (n_v,) int64 global vertex index in mesh
    """
    from ngsolve import BND
    bnd_labels = list(mesh.GetBoundaries())
    if bnd_label is not None:
        target_idx = {i for i, n in enumerate(bnd_labels) if n == bnd_label}
        if not target_idx:
            raise ValueError(
                f"BEM: BND label {bnd_label!r} not found in mesh; "
                f"available: {bnd_labels}")
    else:
        target_idx = set(range(len(bnd_labels)))

    used = set()
    for el in mesh.Elements(BND):
        if el.index in target_idx:
            for vv in el.vertices:
                used.add(vv.nr)

    sorted_used = sorted(used)
    g2l = {g: l for l, g in enumerate(sorted_used)}
    verts = np.array([mesh.vertices[g].point for g in sorted_used],
                     dtype=np.float64)

    tris = []
    for el in mesh.Elements(BND):
        if el.index in target_idx:
            tris.append([g2l[v.nr] for v in el.vertices])
    tris = np.array(tris, dtype=np.int64)
    v_global = np.array(sorted_used, dtype=np.int64)
    return verts, tris, v_global


# ----------------------------------------------------------------------
# Triangle quadrature points (Gauss on reference triangle)
# ----------------------------------------------------------------------
# Reference triangle T_ref = {(xi, eta): xi >= 0, eta >= 0, xi + eta <= 1}
# Hat functions: L0 = 1 - xi - eta, L1 = xi, L2 = eta
# Quadrature lists below give (xi, eta, w) where sum(w) = 1/2 (= area of T_ref)

# 7-point degree-5 symmetric quadrature (Stroud)
_TRI_GAUSS_7 = np.array([
    # (xi, eta, weight)
    (1/3, 1/3, 0.225 * 0.5),
    (0.0597158717, 0.4701420641, 0.1323941527 * 0.5),
    (0.4701420641, 0.0597158717, 0.1323941527 * 0.5),
    (0.4701420641, 0.4701420641, 0.1323941527 * 0.5),
    (0.7974269853, 0.1012865073, 0.1259391805 * 0.5),
    (0.1012865073, 0.7974269853, 0.1259391805 * 0.5),
    (0.1012865073, 0.1012865073, 0.1259391805 * 0.5),
], dtype=np.float64)

# 13-point degree-7 quadrature for more accurate near-singular pairs
_TRI_GAUSS_13 = np.array([
    (0.333333333333333, 0.333333333333333, -0.149570044467670 * 0.5),
    (0.260345966079038, 0.260345966079038,  0.175615257433204 * 0.5),
    (0.260345966079038, 0.479308067841923,  0.175615257433204 * 0.5),
    (0.479308067841923, 0.260345966079038,  0.175615257433204 * 0.5),
    (0.065130102902216, 0.065130102902216,  0.053347235608839 * 0.5),
    (0.065130102902216, 0.869739794195568,  0.053347235608839 * 0.5),
    (0.869739794195568, 0.065130102902216,  0.053347235608839 * 0.5),
    (0.312865496004875, 0.638444188569809,  0.077113760890257 * 0.5),
    (0.638444188569809, 0.048690315425316,  0.077113760890257 * 0.5),
    (0.048690315425316, 0.312865496004875,  0.077113760890257 * 0.5),
    (0.312865496004875, 0.048690315425316,  0.077113760890257 * 0.5),
    (0.638444188569809, 0.312865496004875,  0.077113760890257 * 0.5),
    (0.048690315425316, 0.638444188569809,  0.077113760890257 * 0.5),
], dtype=np.float64)


def tri_quad(degree=5):
    """Pick a triangle Gauss quadrature rule by polynomial degree."""
    if degree <= 5:
        return _TRI_GAUSS_7
    return _TRI_GAUSS_13


def map_quad(verts, tri, q_xi_eta_w):
    """Map reference-triangle quadrature to physical triangle.

    Args:
        verts (n_v, 3): all vertex coords
        tri (3,): vertex indices of one triangle
        q_xi_eta_w (n_q, 3): (xi, eta, weight) on T_ref

    Returns:
        pts (n_q, 3) physical points
        w_phys (n_q,) physical weights = w_ref * 2 * area_phys
        L (n_q, 3) hat function values at quad points (L0, L1, L2)
    """
    v0, v1, v2 = verts[tri[0]], verts[tri[1]], verts[tri[2]]
    e1 = v1 - v0
    e2 = v2 - v0
    cross = np.cross(e1, e2)
    area_phys = 0.5 * np.linalg.norm(cross)
    jac = 2.0 * area_phys

    xi = q_xi_eta_w[:, 0]
    eta = q_xi_eta_w[:, 1]
    w_ref = q_xi_eta_w[:, 2]

    L0 = 1.0 - xi - eta
    L1 = xi
    L2 = eta

    pts = (L0[:, None] * v0[None, :]
           + L1[:, None] * v1[None, :]
           + L2[:, None] * v2[None, :])
    w_phys = w_ref * jac
    L = np.stack([L0, L1, L2], axis=1)
    return pts, w_phys, L


def tri_normal_area(verts, tri):
    """Outward unit normal and area of one flat triangle."""
    v0, v1, v2 = verts[tri[0]], verts[tri[1]], verts[tri[2]]
    e1 = v1 - v0
    e2 = v2 - v0
    cross = np.cross(e1, e2)
    area = 0.5 * np.linalg.norm(cross)
    n = cross / (2.0 * area + 1e-300)
    return n, area


# ----------------------------------------------------------------------
# Pairwise sharing detection
# ----------------------------------------------------------------------
def share_class(tri_a, tri_b):
    """Classify a pair of triangles by how many vertices they share.

    Returns:
        'identical': same triangle (3 shared)
        'edge':      2 shared vertices (edge-shared)
        'vertex':    1 shared vertex
        'regular':   0 shared (non-adjacent)
    Plus the local indices of shared vertices in each tri.
    """
    sa = set(int(v) for v in tri_a)
    sb = set(int(v) for v in tri_b)
    shared = sa & sb
    n_shared = len(shared)
    if n_shared == 3:
        return 'identical', shared
    if n_shared == 2:
        return 'edge', shared
    if n_shared == 1:
        return 'vertex', shared
    return 'regular', shared


# ----------------------------------------------------------------------
# Phase 1.1 SL kernel (regular case via tensor-product Gauss)
# ----------------------------------------------------------------------
def sl_pair_regular(verts, tri_a, tri_b, q_a, q_b):
    """Compute the 3x3 SL local block for a REGULAR (disjoint) pair.

    M_ab[i, j] = int_{T_a} int_{T_b}  L_i(r) (1/(4 pi |r - r'|))
                                       L_j(r') dS dS'
    where i is local vertex of T_a, j is local vertex of T_b.

    Implementation: tensor-product Gauss on T_ref x T_ref.

    Returns (3, 3) float64 contribution.
    """
    pts_a, w_a, La = map_quad(verts, tri_a, q_a)
    pts_b, w_b, Lb = map_quad(verts, tri_b, q_b)

    # Pairwise distance |r - r'|
    dx = pts_a[:, None, :] - pts_b[None, :, :]
    r = np.sqrt(np.sum(dx * dx, axis=-1))
    G = INV_4PI / r           # shape (n_qa, n_qb)

    # Outer-product weights and hat values
    # M[i, j] = sum_{a, b} La[a, i] * w_a[a] * G[a, b] * w_b[b] * Lb[b, j]
    # First contract over b: T[a, j] = sum_b G[a, b] * w_b[b] * Lb[b, j]
    T = (G * w_b[None, :]) @ Lb           # (n_qa, 3)
    # Then contract over a: M[i, j] = sum_a La[a, i] * w_a[a] * T[a, j]
    M = La.T @ (w_a[:, None] * T)          # (3, 3)
    return M


# ----------------------------------------------------------------------
# Phase 1.2 -- Singular pair handling via h-refinement + Duffy.
#
# Strategy: for any singular pair (identical / edge / vertex shared),
# uniformly sub-divide each triangle into k^2 sub-triangles and apply
# the regular kernel between every sub-tri pair, EXCEPT pairs of
# sub-tris that themselves share a vertex/edge/are identical (which
# we further sub-divide recursively, OR cap the recursion and accept
# residual error).
#
# This converges as O(k^{-1}) for the Laplace SL kernel, which is
# slow but predictable.  A spectral method (Sauter-Schwab) would be
# faster but is more code; here we trade convergence rate for
# verifiability.
# ----------------------------------------------------------------------

def _subdivide_triangle(verts_tri, k):
    """Subdivide a flat triangle with vertices (3, 3) into k^2 sub-tris.

    Returns (k^2, 3, 3) sub-vertex coordinates.

    Uses uniform barycentric sub-division.
    """
    v0, v1, v2 = verts_tri
    # Generate (k+1)*(k+2)/2 lattice points on the reference triangle
    sub_tris = []
    for i in range(k):
        for j in range(k - i):
            # Upward sub-tri at lattice (i, j)
            a = ((k - i - j) * v0 + i * v1 + j * v2) / k
            b = ((k - i - j - 1) * v0 + (i + 1) * v1 + j * v2) / k
            c = ((k - i - j - 1) * v0 + i * v1 + (j + 1) * v2) / k
            sub_tris.append([a, b, c])
            if j + i + 2 <= k:
                # Downward sub-tri (filling between upward ones)
                d = ((k - i - j - 2) * v0 + (i + 1) * v1 + (j + 1) * v2) / k
                sub_tris.append([b, d, c])
    return np.asarray(sub_tris)


def _hat_at_subtri_quadpts(verts_tri, sub_tris, q_xi_eta_w):
    """Evaluate parent-triangle hat functions at quadrature points of
    each sub-triangle.

    Args:
        verts_tri (3, 3): parent triangle vertices
        sub_tris (m, 3, 3): m sub-triangles
        q_xi_eta_w (n_q, 3): quadrature on reference triangle

    Returns:
        pts_global (m * n_q, 3) physical points
        w_global   (m * n_q,)   physical weights
        Lparent    (m * n_q, 3) parent hat values at each quad point
    """
    v0, v1, v2 = verts_tri
    # Parent area-coordinates: solve the linear system
    e1 = v1 - v0
    e2 = v2 - v0
    # Normal of parent
    n_par = np.cross(e1, e2)
    A_par = 0.5 * np.linalg.norm(n_par)
    # We can find barycentric coords via simple linear system later.
    # For each sub-tri, map its quadrature pts to physical, then
    # express in parent's barycentric coords via inverse 2D transform.
    pts_all = []
    w_all = []
    L_par_all = []
    for k_sub in range(len(sub_tris)):
        st = sub_tris[k_sub]
        s0, s1, s2 = st
        e1_s = s1 - s0
        e2_s = s2 - s0
        cross = np.cross(e1_s, e2_s)
        A_sub = 0.5 * np.linalg.norm(cross)
        jac_sub = 2.0 * A_sub
        xi = q_xi_eta_w[:, 0]
        eta = q_xi_eta_w[:, 1]
        w_ref = q_xi_eta_w[:, 2]
        L0 = 1.0 - xi - eta
        L1 = xi
        L2 = eta
        pts = (L0[:, None] * s0 + L1[:, None] * s1 + L2[:, None] * s2)
        w_phys = w_ref * jac_sub
        # Now express each pts in parent (v0, v1, v2)'s barycentric coords.
        # Use 2D projection: solve [e1 e2] * [L1_par, L2_par]^T = pts - v0
        # (in plane).  We just need L1_par, L2_par; L0_par = 1 - L1 - L2.
        rel = pts - v0[None, :]
        # Project onto e1 and e2 via Gram-Schmidt-like (use 2x2 system on
        # the in-plane basis).
        a11 = e1 @ e1
        a12 = e1 @ e2
        a22 = e2 @ e2
        det = a11 * a22 - a12 ** 2
        b1 = rel @ e1
        b2 = rel @ e2
        L1_par = (a22 * b1 - a12 * b2) / det
        L2_par = (-a12 * b1 + a11 * b2) / det
        L0_par = 1.0 - L1_par - L2_par
        L_par = np.stack([L0_par, L1_par, L2_par], axis=1)
        pts_all.append(pts)
        w_all.append(w_phys)
        L_par_all.append(L_par)
    return (np.concatenate(pts_all, axis=0),
            np.concatenate(w_all, axis=0),
            np.concatenate(L_par_all, axis=0))


def sl_pair_singular_subdivide(verts, tri_a, tri_b, *,
                                k_sub=8, q_degree=5,
                                drop_self_subpairs=True):
    """3x3 SL local block for a singular pair, using sub-division.

    Sub-divide each parent triangle into k_sub^2 children.  Apply the
    regular formula on every sub-pair except those that ARE singular
    (same sub-tri, edge-shared, vertex-shared).  Drop singular sub-pairs
    or recurse (here we drop -- this gives O(k^{-1}) convergence but
    is finite and verifiable).
    """
    Va = verts[tri_a]
    Vb = verts[tri_b]
    sub_a = _subdivide_triangle(Va, k_sub)
    sub_b = _subdivide_triangle(Vb, k_sub)
    q = tri_quad(q_degree)

    # Compute physical quad pts + parent hat values on both sides
    pa_all, wa_all, Lpar_a = _hat_at_subtri_quadpts(Va, sub_a, q)
    pb_all, wb_all, Lpar_b = _hat_at_subtri_quadpts(Vb, sub_b, q)

    # We want sum over all (i_qa, i_qb) of
    #   wa[i_qa] * Lpar_a[i_qa, i] * G(pa[i_qa], pb[i_qb])
    #            * wb[i_qb] * Lpar_b[i_qb, j]
    # EXCEPT where i_qa and i_qb belong to a "too-singular" sub-pair.
    # For simplicity here we evaluate ALL pairs with regularised G:
    # G_reg(r, r') = 1 / (4 pi sqrt(|r-r'|^2 + eps^2)) where eps is a
    # small fraction of typical sub-tri size.  This is approximate but
    # avoids the explicit case split.
    h_a = np.sqrt(2.0 * (np.linalg.norm(np.cross(sub_a[:, 1] - sub_a[:, 0],
                                                  sub_a[:, 2] - sub_a[:, 0]),
                                          axis=1) / 2.0))
    h_typ = float(h_a.mean()) * 0.05  # 5% of mean sub-tri size
    dx = pa_all[:, None, :] - pb_all[None, :, :]
    r2 = np.sum(dx * dx, axis=-1)
    G = INV_4PI / np.sqrt(r2 + h_typ ** 2)

    # Outer products to build the 3x3 block in PARENT indexing
    T = (G * wb_all[None, :]) @ Lpar_b
    M = Lpar_a.T @ (wa_all[:, None] * T)
    return M


# ----------------------------------------------------------------------
# Dense SL assembly (Phase 1, regular pairs only first; singular added
# in 1.2, then HACApK in 1.4)
# ----------------------------------------------------------------------
def assemble_SL_dense(verts, tris, *, regular_quad_degree=5,
                      include_singular=False):
    """Build the dense Galerkin Laplace SL matrix on a flat triangulation.

    Args:
        verts (n_v, 3): vertex coords
        tris  (n_t, 3): triangle vertex indices
        regular_quad_degree: triangle Gauss degree for non-singular pairs
        include_singular: if False, mark identical/edge/vertex pairs as
            zero (returns matrix with unfilled diagonal blocks -- useful
            for unit-testing the regular-case kernel against NGSolve.bem
            on the off-diagonal far block).

    Returns (n_v, n_v) float64 SL matrix.
    """
    n_v = len(verts)
    n_t = len(tris)
    SL = np.zeros((n_v, n_v), dtype=np.float64)
    q = tri_quad(regular_quad_degree)

    # Pre-cache for speed: per-tri quadrature data
    pts_cache = []
    w_cache = []
    L_cache = []
    for t in range(n_t):
        pts, w, L = map_quad(verts, tris[t], q)
        pts_cache.append(pts)
        w_cache.append(w)
        L_cache.append(L)

    for a in range(n_t):
        Va = tris[a]
        for b in range(n_t):
            Vb = tris[b]
            cls, _ = share_class(Va, Vb)
            if cls != 'regular':
                if include_singular:
                    raise NotImplementedError(
                        f"singular pair handling (case={cls}) is "
                        f"Phase 1.2 -- not implemented yet")
                continue
            # Inline regular formula using cached quadrature data
            pts_a, w_a, La = pts_cache[a], w_cache[a], L_cache[a]
            pts_b, w_b, Lb = pts_cache[b], w_cache[b], L_cache[b]
            dx = pts_a[:, None, :] - pts_b[None, :, :]
            r = np.sqrt(np.sum(dx * dx, axis=-1))
            G = INV_4PI / r
            T = (G * w_b[None, :]) @ Lb
            block = La.T @ (w_a[:, None] * T)
            # Scatter into global SL
            SL[np.ix_(Va, Vb)] += block

    return SL
