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

# Higher-order triangle Gauss-Legendre rules generated on demand.
# Use a Duffy-collapsed product of two 1D Gauss-Legendre rules:
#   xi  = s,  eta = s * t   with (s, t) ∈ [0, 1]^2, Jacobian s.
# This integrates exactly any polynomial of total degree (2n-1) in (xi, eta)
# for n-point 1D rules in each direction.  Sum of weights = 1/2 (= area of T_ref).
_TRI_GL_CACHE = {}


def _tri_gl_duffy(n):
    """Duffy-collapsed n*n point Gauss-Legendre quadrature on T_ref.

    Mapping: (xi, eta) = (s*(1-t), s*t), Jacobian = s.  Maps the unit
    square (s, t) ∈ [0,1]^2 onto T_ref = {xi >= 0, eta >= 0, xi+eta<=1}
    (xi+eta = s).  Sum of weights = 1/2 (= area of T_ref).

    Returns (n*n, 3) array of (xi, eta, w).
    Exact for polynomials of total degree 2n - 1 in (xi, eta).
    """
    if n in _TRI_GL_CACHE:
        return _TRI_GL_CACHE[n]
    x, w = np.polynomial.legendre.leggauss(n)
    s_pts = 0.5 * (x + 1.0); s_w = 0.5 * w
    t_pts = 0.5 * (x + 1.0); t_w = 0.5 * w
    out = []
    for i in range(n):
        for j in range(n):
            s = s_pts[i]; t = t_pts[j]
            xi = s * (1.0 - t)
            eta = s * t
            wt = s * s_w[i] * t_w[j]
            out.append((xi, eta, wt))
    arr = np.array(out, dtype=np.float64)
    _TRI_GL_CACHE[n] = arr
    return arr


def tri_quad(degree=5):
    """Pick a triangle Gauss quadrature rule by polynomial degree.

    For degree <= 5 (resp. <= 7), use the Stroud symmetric 7-pt (resp.
    13-pt) rule.  For higher requested degree, fall back to a
    Duffy-collapsed n x n Gauss-Legendre product with n = ceil((degree+1)/2),
    which is exact through total degree 2n - 1.
    """
    if degree <= 5:
        return _TRI_GAUSS_7
    if degree <= 7:
        return _TRI_GAUSS_13
    n = (degree + 1 + 1) // 2  # ceil((degree+1)/2), giving exact 2n-1
    return _tri_gl_duffy(n)


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
# Phase 1.2 -- Sauter-Schwab Duffy 4D for singular pairs.
#
# Direct port of NGSolve.bem's intrules_SauterSchwab.cpp (CommonVertex
# 2 sub-cubes, CommonEdge 5 sub-cubes, IdenticPanel 6 sub-cubes).
# Reference T_ref = [(0,0), (1,0), (0,1)] (xi >= 0, eta >= 0, xi+eta <= 1)
# matches NGSolve.  Hat functions:
#   L_0(xi, eta) = 1 - xi - eta
#   L_1(xi, eta) = xi
#   L_2(xi, eta) = eta
#
# Each rule produces (xi_a, eta_a), (xi_b, eta_b), and a quadrature
# weight w_q.  The pair-wise SL Galerkin entry is then
#   M_ij = sum_q L_i(xi_a^q, eta_a^q) * (1/(4*pi*|r_a^q - r_b^q|))
#                * L_j(xi_b^q, eta_b^q) * (4 * A_a * A_b) * w_q
# where (4*A_a*A_b) is the surface-area Jacobian for both T_ref -> T_phys
# maps, and the singular factor (rho-equivalent) is already absorbed in
# w_q via the Sauter-Schwab Jacobian.
# ----------------------------------------------------------------------

def _gl1d(n):
    """n-point Gauss-Legendre on [0, 1]."""
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w


_SS_CACHE = {}  # (case, n_q) -> (xi_a, eta_a, xi_b, eta_b, weight) all flat


def _ss_quad_pts(case, n_q):
    """Build Sauter-Schwab quadrature nodes for one of the singular cases.

    case ∈ {"common_vertex", "common_edge", "identical"}.  Direct port of
    NGSolve.bem CommonVertexIntegrationRule / CommonEdgeIntegrationRule /
    IdenticPanelIntegrationRule.

    The quadrature is a 4D tensor product of n_q-point 1D Gauss-Legendre,
    multiplied by 2 / 5 / 6 sub-cubes.

    Returns:
        xi_a, eta_a, xi_b, eta_b: (N,) arrays of points on T_ref^a x T_ref^b
        weight: (N,) array of quadrature weights such that
                 sum_q f(xi_a, eta_a, xi_b, eta_b) * weight_q
                 approximates the singular integral
                 int_{T_ref} int_{T_ref} f * (1/|...|) dxi_a deta_a dxi_b deta_b
                 EXCEPT that the Sauter-Schwab xi^3 and similar Jacobian
                 factors are NOT pre-divided by |...| -- we supply weight
                 = J(xi, e1, e2, e3) * w_xi * w_e1 * w_e2 * w_e3 directly.
        Caller multiplies by 1/|r_a - r_b| (with the r_a, r_b computed
        from xi_a, eta_a, xi_b, eta_b) and the physical area Jacobian
        (4 * A_a * A_b) and the kernel prefactor (1/(4*pi)).
    """
    key = (case, n_q)
    if key in _SS_CACHE:
        return _SS_CACHE[key]

    x_q, w_q = _gl1d(n_q)
    # NGSolve: irhex (e1, e2, e3) outer x irsegm (xi) inner.
    # We build all 4D combinations.  Use indexing='ij' to mimic
    # the C++ nested-for ordering, but that doesn't actually matter
    # since we sum over all of them.
    e1_g, e2_g, e3_g, xi_g = np.meshgrid(x_q, x_q, x_q, x_q, indexing='ij')
    we1, we2, we3, wxi = np.meshgrid(w_q, w_q, w_q, w_q, indexing='ij')
    e1 = e1_g.ravel(); e2 = e2_g.ravel(); e3 = e3_g.ravel(); xi = xi_g.ravel()
    w_4d = (we1 * we2 * we3 * wxi).ravel()

    # Each sub-cube contributes (ip0, ip1, ip2, ip3) and a Jacobian J.
    # NGSolve trafo: (xi_a, eta_a) = (ip0 - ip1, ip1),
    #                (xi_b, eta_b) = (ip2 - ip3, ip3).
    sub_cubes = []
    if case == "common_vertex":
        # Duffies[0] = xi*(1, e1, e2, e2*e3), J = xi^3 * e2
        # Duffies[1] = xi*(e2, e2*e3, 1, e1)
        sub_cubes.append((xi*1.0,    xi*e1,    xi*e2,    xi*e2*e3, xi**3 * e2))
        sub_cubes.append((xi*e2,     xi*e2*e3, xi*1.0,   xi*e1,    xi**3 * e2))
    elif case == "common_edge":
        # 5 sub-cubes
        # Duffies[0] = xi*(1, e1*e3, 1-e1*e2, e1*(1-e2)),    J = xi^3 * e1^2
        # Duffies[1] = xi*(1, e1, 1-e1*e2*e3, e1*e2*(1-e3)), J = xi^3 * e1^2 * e2
        # Duffies[2] = xi*(1-e1*e2, e1*(1-e2), 1, e1*e2*e3), J = xi^3 * e1^2 * e2
        # Duffies[3] = xi*(1-e1*e2*e3, e1*e2*(1-e3), 1, e1), J = xi^3 * e1^2 * e2
        # Duffies[4] = xi*(1-e1*e2*e3, e1*(1-e2*e3), 1, e1*e2), J = xi^3 * e1^2 * e2
        J_a = xi**3 * e1**2
        J_b = xi**3 * e1**2 * e2
        sub_cubes.append((xi*1.0,           xi*e1*e3,          xi*(1-e1*e2),     xi*e1*(1-e2),       J_a))
        sub_cubes.append((xi*1.0,           xi*e1,             xi*(1-e1*e2*e3),  xi*e1*e2*(1-e3),    J_b))
        sub_cubes.append((xi*(1-e1*e2),     xi*e1*(1-e2),      xi*1.0,           xi*e1*e2*e3,        J_b))
        sub_cubes.append((xi*(1-e1*e2*e3),  xi*e1*e2*(1-e3),   xi*1.0,           xi*e1,              J_b))
        sub_cubes.append((xi*(1-e1*e2*e3),  xi*e1*(1-e2*e3),   xi*1.0,           xi*e1*e2,           J_b))
    elif case == "identical":
        # 6 sub-cubes, J = xi^3 * e1^2 * e2  for all
        J = xi**3 * e1**2 * e2
        # Duffies[0] = xi*(1, 1-e1+e1*e2, 1-e1*e2*e3, 1-e1)
        # Duffies[1] = xi*(1-e1*e2*e3, 1-e1, 1, 1-e1+e1*e2)
        # Duffies[2] = xi*(1, e1*(1-e2+e2*e3), 1-e1*e2, e1*(1-e2))
        # Duffies[3] = xi*(1-e1*e2, e1*(1-e2), 1, e1*(1-e2+e2*e3))
        # Duffies[4] = xi*(1-e1*e2*e3, e1*(1-e2*e3), 1, e1*(1-e2))
        # Duffies[5] = xi*(1, e1*(1-e2), 1-e1*e2*e3, e1*(1-e2*e3))
        sub_cubes.append((xi*1.0,           xi*(1-e1+e1*e2),    xi*(1-e1*e2*e3),  xi*(1-e1),            J))
        sub_cubes.append((xi*(1-e1*e2*e3),  xi*(1-e1),          xi*1.0,           xi*(1-e1+e1*e2),      J))
        sub_cubes.append((xi*1.0,           xi*e1*(1-e2+e2*e3), xi*(1-e1*e2),     xi*e1*(1-e2),         J))
        sub_cubes.append((xi*(1-e1*e2),     xi*e1*(1-e2),       xi*1.0,           xi*e1*(1-e2+e2*e3),   J))
        sub_cubes.append((xi*(1-e1*e2*e3),  xi*e1*(1-e2*e3),    xi*1.0,           xi*e1*(1-e2),         J))
        sub_cubes.append((xi*1.0,           xi*e1*(1-e2),       xi*(1-e1*e2*e3),  xi*e1*(1-e2*e3),      J))
    else:
        raise ValueError(f"unknown SS case: {case!r}")

    xi_a_all = []
    eta_a_all = []
    xi_b_all = []
    eta_b_all = []
    weight_all = []
    for ip0, ip1, ip2, ip3, J_per in sub_cubes:
        # NGSolve trafo
        xi_a_all.append(ip0 - ip1)
        eta_a_all.append(ip1)
        xi_b_all.append(ip2 - ip3)
        eta_b_all.append(ip3)
        weight_all.append(J_per * w_4d)

    out = (np.concatenate(xi_a_all),
           np.concatenate(eta_a_all),
           np.concatenate(xi_b_all),
           np.concatenate(eta_b_all),
           np.concatenate(weight_all))
    _SS_CACHE[key] = out
    return out


def _ss_block(va, vb, case, n_q):
    """Common-vertex/edge/identical Galerkin SL block in PERMUTED indexing.

    Both va and vb are in NGSolve canonical permutation:
      common_vertex: shared vertex at local 0 of both
      common_edge:   shared edge from local 0 to local 1, with vertex 0 of
                     T_a == vertex 0 of T_b and vertex 1 of T_a == vertex 1
                     of T_b (same edge orientation)
      identical:     T_a == T_b (same vertex order)

    Returns (3, 3) block where rows are L_i^a (canonical local) and cols
    are L_j^b (canonical local).
    """
    e1a = va[1] - va[0]
    e2a = va[2] - va[0]
    e1b = vb[1] - vb[0]
    e2b = vb[2] - vb[0]
    Aa = 0.5 * np.linalg.norm(np.cross(e1a, e2a))
    Ab = 0.5 * np.linalg.norm(np.cross(e1b, e2b))

    xi_a, eta_a, xi_b, eta_b, w_q = _ss_quad_pts(case, n_q)

    # Physical positions
    r_a = va[0] + xi_a[:, None] * e1a + eta_a[:, None] * e2a
    r_b = vb[0] + xi_b[:, None] * e1b + eta_b[:, None] * e2b
    diff = r_a - r_b
    r_dist = np.sqrt(np.sum(diff * diff, axis=1))
    # Guard against exact zero (numerically vanishingly rare on physical pairs)
    r_dist = np.where(r_dist == 0, 1e-300, r_dist)

    # Hat values at quadrature points
    L0_a = 1.0 - xi_a - eta_a
    L1_a = xi_a
    L2_a = eta_a
    L0_b = 1.0 - xi_b - eta_b
    L1_b = xi_b
    L2_b = eta_b
    La = np.stack([L0_a, L1_a, L2_a], axis=1)
    Lb = np.stack([L0_b, L1_b, L2_b], axis=1)

    # Galerkin entry weight per quad point:
    #   1/(4*pi*r_dist) * w_q * (4 * A_a * A_b)   [4 from 2*A_a * 2*A_b]
    # = w_q * A_a * A_b / (pi * r_dist)
    weight = w_q * Aa * Ab / (np.pi * r_dist)
    return La.T @ (weight[:, None] * Lb)


def _ss_common_vertex_block(va, vb, share_a, share_b, n_q=5):
    """SL block for vertex-shared pair, returning (3,3) in ORIGINAL indexing.

    Permutes va, vb so the shared vertex is at local 0, applies
    Sauter-Schwab common-vertex rule, then unpermutes.
    """
    pa = (share_a, (share_a + 1) % 3, (share_a + 2) % 3)
    pb = (share_b, (share_b + 1) % 3, (share_b + 2) % 3)
    va_p = np.array([va[pa[0]], va[pa[1]], va[pa[2]]])
    vb_p = np.array([vb[pb[0]], vb[pb[1]], vb[pb[2]]])
    M_perm = _ss_block(va_p, vb_p, "common_vertex", n_q)
    M_orig = np.zeros((3, 3))
    M_orig[np.ix_(pa, pb)] = M_perm
    return M_orig


def _ss_common_edge_block(va, vb, shared_pairs, n_q=5):
    """SL block for edge-shared pair via Sauter-Schwab common-edge rule.

    shared_pairs: list of 2 (la, lb) tuples; va[la] == vb[lb] for each.
    Returns (3, 3) block in ORIGINAL indexing.

    Permutation: place both shared vertices at local 0 and 1 of the
    permuted triangles, with consistent edge orientation
    (va_p[0] == vb_p[0], va_p[1] == vb_p[1]).
    """
    (la0, lb0), (la1, lb1) = shared_pairs
    apex_a = ({0, 1, 2} - {la0, la1}).pop()
    apex_b = ({0, 1, 2} - {lb0, lb1}).pop()
    pa = (la0, la1, apex_a)
    pb = (lb0, lb1, apex_b)
    va_p = np.array([va[pa[0]], va[pa[1]], va[pa[2]]])
    vb_p = np.array([vb[pb[0]], vb[pb[1]], vb[pb[2]]])
    M_perm = _ss_block(va_p, vb_p, "common_edge", n_q)
    M_orig = np.zeros((3, 3))
    M_orig[np.ix_(pa, pb)] = M_perm
    return M_orig


def _ss_identical_block(va, vb, shared_pairs, n_q=5):
    """SL block for identical pair via Sauter-Schwab IdenticPanel rule.

    Vertex sets of T_a and T_b are equal but possibly permuted.
    Returns (3, 3) block in ORIGINAL indexing of va (rows), vb (cols).
    """
    perm_b = [None, None, None]
    for la, lb in shared_pairs:
        perm_b[la] = lb
    # vb_in_va_order[la] = vb[perm_b[la]] == va[la]
    vb_in_va_order = np.array([vb[perm_b[0]], vb[perm_b[1]], vb[perm_b[2]]])
    M_va = _ss_block(va, vb_in_va_order, "identical", n_q)
    # M_va[i, k] is in (va index, va_order_b index).  Map va_order_b -> vb_orig.
    M_orig = np.zeros((3, 3))
    for la in range(3):
        M_orig[:, perm_b[la]] = M_va[:, la]
    return M_orig


def sl_pair_singular(verts, tri_a, tri_b, n_q=5):
    """Public Phase 1.2 entry: SL block for any singular pair.

    Dispatches to vertex / edge / identical handler.

    Args:
        verts (n_v, 3): mesh vertex coords.
        tri_a, tri_b: (3,) global vertex indices of T_a, T_b.
        n_q: 1D Gauss-Legendre order (n_q^4 nodes per sub-cube; 2 / 5 / 6
            sub-cubes for vertex / edge / identical).

    Returns:
        (3, 3) SL block in tri_a (rows) and tri_b (cols) local indexing.
    """
    va = verts[tri_a]
    vb = verts[tri_b]
    cls, _shared = share_class(tri_a, tri_b)
    shared_pairs = []
    for la in range(3):
        for lb in range(3):
            if int(tri_a[la]) == int(tri_b[lb]):
                shared_pairs.append((la, lb))

    if cls == 'vertex':
        la, lb = shared_pairs[0]
        return _ss_common_vertex_block(va, vb, share_a=la, share_b=lb, n_q=n_q)
    if cls == 'edge':
        return _ss_common_edge_block(va, vb, shared_pairs, n_q=n_q)
    if cls == 'identical':
        return _ss_identical_block(va, vb, shared_pairs, n_q=n_q)
    raise ValueError(f"sl_pair_singular called on non-singular pair: cls={cls}")


# ----------------------------------------------------------------------
# Dense SL assembly (Phase 1, regular pairs only first; singular added
# in 1.2, then HACApK in 1.4)
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Phase 1.3 -- Laplace double-layer (DL) kernel.
#
# Galerkin DL convention (matches NGSolve.bem LaplaceDL):
#   M_ij = int_{T_a} int_{T_b} L_i^a(r)
#                              [ -(r - r')·n_y / (4*pi |r-r'|^3) ]
#                              L_j^b(r')   dS dS'
# where n_y is the OUTWARD normal of T_b at r'.
#
# Same Sauter-Schwab Duffy 4D quadrature as SL; the kernel substitution
# is the only difference.  The (r-r') · n_y / |r-r'|^3 = O(1/r^2)
# singularity is absorbed by the same xi^3 Jacobian (since (r-r')·n_y =
# O(rho * r^something) typically vanishes at the shared singularity).
# ----------------------------------------------------------------------


def _ss_block_dl(va, vb, case, n_q, n_b):
    """Common-vertex/edge/identical Galerkin DL block in PERMUTED indexing.

    n_b: outward unit normal of T_b (constant on a flat triangle).

    DL Galerkin entry:
        M[i, j] = sum_q L_i^a * [-(r_a - r_b)·n_b / (4*pi*r_dist^3)]
                                  * L_j^b * (4*A_a*A_b) * w_q
    """
    e1a = va[1] - va[0]
    e2a = va[2] - va[0]
    e1b = vb[1] - vb[0]
    e2b = vb[2] - vb[0]
    Aa = 0.5 * np.linalg.norm(np.cross(e1a, e2a))
    Ab = 0.5 * np.linalg.norm(np.cross(e1b, e2b))

    xi_a, eta_a, xi_b, eta_b, w_q = _ss_quad_pts(case, n_q)

    r_a = va[0] + xi_a[:, None] * e1a + eta_a[:, None] * e2a
    r_b = vb[0] + xi_b[:, None] * e1b + eta_b[:, None] * e2b
    diff = r_a - r_b
    r2 = np.sum(diff * diff, axis=1)
    r_dist = np.sqrt(r2)
    r_dist = np.where(r_dist == 0, 1e-300, r_dist)

    L0_a = 1.0 - xi_a - eta_a
    L1_a = xi_a
    L2_a = eta_a
    L0_b = 1.0 - xi_b - eta_b
    L1_b = xi_b
    L2_b = eta_b
    La = np.stack([L0_a, L1_a, L2_a], axis=1)
    Lb = np.stack([L0_b, L1_b, L2_b], axis=1)

    # DL kernel ∂G/∂n_y, NGSolve.bem convention (positive sign):
    # ∂G/∂n_y = (r - r') · n_y / (4*pi * |r-r'|^3).
    # Combined with the 4*A_a*A_b = 2*A_a * 2*A_b area Jacobian:
    #   weight per quad pt = (r-r')·n_y / (4pi r^3) * 4*A_a*A_b * w_q
    #                      = (r-r')·n_y * A_a * A_b / (pi * r^3) * w_q
    dot_n = diff @ n_b
    weight = dot_n / (np.pi * r_dist**3) * w_q * Aa * Ab
    return La.T @ (weight[:, None] * Lb)


def _outward_normal(verts, tri):
    """Unit normal of a flat triangle, oriented per the vertex ordering."""
    v0, v1, v2 = verts[tri[0]], verts[tri[1]], verts[tri[2]]
    cross = np.cross(v1 - v0, v2 - v0)
    n = cross / (np.linalg.norm(cross) + 1e-300)
    return n


def dl_pair_singular(verts, tri_a, tri_b, n_q=8):
    """DL block for a singular pair (any sharing class).

    Args:
        verts (n_v, 3): vertex coords.
        tri_a, tri_b: (3,) global vertex indices.
        n_q: 1D Gauss-Legendre order (n_q^4 nodes per sub-cube).

    Returns:
        (3, 3) DL block in tri_a (rows) x tri_b (cols) local indexing.
    """
    va = verts[tri_a]
    vb = verts[tri_b]
    n_b = _outward_normal(verts, tri_b)
    cls, _ = share_class(tri_a, tri_b)
    shared_pairs = [(la, lb) for la in range(3) for lb in range(3)
                    if int(tri_a[la]) == int(tri_b[lb])]

    if cls == 'vertex':
        la, lb = shared_pairs[0]
        pa = (la, (la + 1) % 3, (la + 2) % 3)
        pb = (lb, (lb + 1) % 3, (lb + 2) % 3)
        va_p = np.array([va[pa[0]], va[pa[1]], va[pa[2]]])
        vb_p = np.array([vb[pb[0]], vb[pb[1]], vb[pb[2]]])
        M_perm = _ss_block_dl(va_p, vb_p, "common_vertex", n_q, n_b)
        M_orig = np.zeros((3, 3))
        M_orig[np.ix_(pa, pb)] = M_perm
        return M_orig
    if cls == 'edge':
        (la0, lb0), (la1, lb1) = shared_pairs
        apex_a = ({0, 1, 2} - {la0, la1}).pop()
        apex_b = ({0, 1, 2} - {lb0, lb1}).pop()
        pa = (la0, la1, apex_a)
        pb = (lb0, lb1, apex_b)
        va_p = np.array([va[pa[0]], va[pa[1]], va[pa[2]]])
        vb_p = np.array([vb[pb[0]], vb[pb[1]], vb[pb[2]]])
        M_perm = _ss_block_dl(va_p, vb_p, "common_edge", n_q, n_b)
        M_orig = np.zeros((3, 3))
        M_orig[np.ix_(pa, pb)] = M_perm
        return M_orig
    if cls == 'identical':
        # For identical pair, (r - r') · n_b = 0 since both r and r' lie on
        # the same plane.  So DL kernel is identically zero.  But the limit
        # from one side gives the well-known +/- 1/2 jump.  The Cauchy
        # principal value of the DL integral (excluding the diagonal jump)
        # is zero for a flat triangle on itself.
        # NGSolve.bem returns the Cauchy PV (consistent with the
        # 1/2 - DL identity on closed surfaces).
        return np.zeros((3, 3))
    raise ValueError(f"dl_pair_singular: unexpected class {cls}")


def assemble_DL_dense(verts, tris, *, regular_quad_degree=11,
                       include_singular=True, singular_n_q=8):
    """Build the dense Galerkin Laplace DL matrix on a flat triangulation.

    Convention matches NGSolve.bem LaplaceDL:
        M_ij = int_{T_a} int_{T_b} L_i^a(r) [-(r-r')·n_y / (4*pi*|r-r'|^3)]
                                  L_j^b(r') dS dS'

    For the diagonal (identical) pair on a flat triangle, the kernel
    (r - r') · n_y vanishes identically, so the Cauchy PV is zero.

    Args: as for assemble_SL_dense.
    """
    n_v = len(verts)
    n_t = len(tris)
    DL = np.zeros((n_v, n_v), dtype=np.float64)
    q = tri_quad(regular_quad_degree)

    # Pre-cache per-tri quadrature data + outward normals
    pts_cache = []
    w_cache = []
    L_cache = []
    n_cache = []
    for t in range(n_t):
        pts, w, L = map_quad(verts, tris[t], q)
        pts_cache.append(pts)
        w_cache.append(w)
        L_cache.append(L)
        n_cache.append(_outward_normal(verts, tris[t]))

    for a in range(n_t):
        Va = tris[a]
        for b in range(n_t):
            Vb = tris[b]
            cls, _ = share_class(Va, Vb)
            n_b = n_cache[b]
            if cls == 'regular':
                pts_a, w_a, La = pts_cache[a], w_cache[a], L_cache[a]
                pts_b, w_b, Lb = pts_cache[b], w_cache[b], L_cache[b]
                dx = pts_a[:, None, :] - pts_b[None, :, :]
                r2 = np.sum(dx * dx, axis=-1)
                r = np.sqrt(r2)
                # ∂G/∂n_y = (r-r')·n_y / (4*pi*|r-r'|^3) (NGSolve convention)
                dot_n = dx @ n_b   # (n_qa, n_qb)
                K = dot_n / (4.0 * np.pi * r2 * r)
                # Weights: w_a (n_qa,), w_b (n_qb,)
                # block[i, j] = sum_a sum_b La[a, i] * w_a[a] * K[a, b] * w_b[b] * Lb[b, j]
                Tmat = (K * w_b[None, :]) @ Lb
                block = La.T @ (w_a[:, None] * Tmat)
            else:
                if not include_singular:
                    continue
                block = dl_pair_singular(verts, Va, Vb, n_q=singular_n_q)
            DL[np.ix_(Va, Vb)] += block

    return DL


# ----------------------------------------------------------------------
# Phase 1.5 -- end-to-end pipeline: build SL, DL, M, K from a NGSolve mesh
# without invoking ngsolve.bem.  Wires into bem_sibc_solver.py as a
# drop-in replacement.
# ----------------------------------------------------------------------


def assemble_bem_matrices(mesh, *, bnd_label=None,
                            regular_quad_degree=11, singular_n_q=8):
    """Build SL, DL, M (surface mass), K (Laplace-Beltrami) for one BND
    label of an NGSolve mesh, using our in-tree Sauter-Schwab BEM.

    Args:
        mesh: NGSolve Mesh.  May be a volume mesh (will use its BND elements)
              or a surface mesh (BND elements are the surface tris).
        bnd_label: optional string name of the boundary label to keep.
                   If None, all BND elements are used.
        regular_quad_degree: outer Gauss degree for non-singular tri pairs.
        singular_n_q: 1D Gauss-Legendre order for the Sauter-Schwab quadrature.

    Returns:
        dict with keys:
            SL  (ndof, ndof) Galerkin Laplace single-layer matrix
            DL  (ndof, ndof) Galerkin Laplace double-layer matrix
            M   (ndof, ndof) surface mass matrix (P1 hat)
            K   (ndof, ndof) Laplace-Beltrami stiffness matrix
            verts (n_v, 3) vertex coords (P1 DOF coordinates)
            tris  (n_t, 3) triangle vertex indices into verts
            v_global (n_v,) the original mesh.vertices index for each
                local vertex, so the caller can map to/from full-mesh
                vector layouts when ``bnd_label`` was used.
    """
    verts, tris, v_global = extract_surface(mesh, bnd_label=bnd_label)
    n_v = len(verts)
    n_t = len(tris)

    SL = assemble_SL_dense(verts, tris,
                            regular_quad_degree=regular_quad_degree,
                            include_singular=True,
                            singular_n_q=singular_n_q)
    DL = assemble_DL_dense(verts, tris,
                            regular_quad_degree=regular_quad_degree,
                            include_singular=True,
                            singular_n_q=singular_n_q)

    # Surface mass and Laplace-Beltrami for P1 hat functions on flat tris.
    # Both have closed-form per-element 3x3 blocks.
    #   M_local = (Area / 12) * [[2,1,1],[1,2,1],[1,1,2]]
    #   K_local = (1 / (4*Area)) * (B^T B), where B[k,:] = (v[k+1] - v[k+2])
    #             rotated by 90deg in the triangle's plane.  Equivalent:
    #   K_local[i,j] = -(e_i . e_j) / (4 Area)  where e_i is the OPPOSITE
    #             edge of vertex i, AS A 3D VECTOR.  Also equals
    #             cot(theta_k) / 2 between vertices i, j -- the standard
    #             cotangent formula.
    M = np.zeros((n_v, n_v))
    K = np.zeros((n_v, n_v))
    for t in range(n_t):
        v0, v1, v2 = verts[tris[t]]
        e0 = v2 - v1   # opposite vertex 0
        e1 = v0 - v2   # opposite vertex 1
        e2 = v1 - v0   # opposite vertex 2
        edges = (e0, e1, e2)
        Acr = np.cross(v1 - v0, v2 - v0)
        A = 0.5 * np.linalg.norm(Acr)

        # Mass
        M_loc = (A / 12.0) * np.array([[2, 1, 1],
                                        [1, 2, 1],
                                        [1, 1, 2]], dtype=np.float64)
        # Laplace-Beltrami: K_ij = (e_i . e_j) / (4 A)
        # (e_i is the edge OPPOSITE vertex i.  Off-diag ij sign as in the
        # standard cotan-Laplace formula; diagonal is sum of -off-diags.)
        K_loc = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                K_loc[i, j] = (edges[i] @ edges[j]) / (4.0 * A)

        idx = tris[t]
        M[np.ix_(idx, idx)] += M_loc
        K[np.ix_(idx, idx)] += K_loc

    return {
        "SL": SL, "DL": DL, "M": M, "K": K,
        "verts": verts, "tris": tris, "v_global": v_global,
    }


def assemble_SL_dense(verts, tris, *, regular_quad_degree=11,
                      include_singular=True, singular_n_q=8):
    """Build the dense Galerkin Laplace SL matrix on a flat triangulation.

    Args:
        verts (n_v, 3): vertex coords
        tris  (n_t, 3): triangle vertex indices
        regular_quad_degree: triangle Gauss degree for non-singular pairs
            (5 = 7-pt; 7 = 13-pt)
        include_singular: if True (default), use Sauter-Schwab Duffy on
            identical/edge/vertex pairs.  If False, set them to zero
            (used by Phase 1.1 verification harness).
        singular_n_q: 1D Gauss-Legendre order for the Sauter-Schwab 4D
            tensor product (default 5 -- adequate for ~1e-8 accuracy on
            smooth post-Duffy integrands).

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
            if cls == 'regular':
                pts_a, w_a, La = pts_cache[a], w_cache[a], L_cache[a]
                pts_b, w_b, Lb = pts_cache[b], w_cache[b], L_cache[b]
                dx = pts_a[:, None, :] - pts_b[None, :, :]
                r = np.sqrt(np.sum(dx * dx, axis=-1))
                G = INV_4PI / r
                Tmat = (G * w_b[None, :]) @ Lb
                block = La.T @ (w_a[:, None] * Tmat)
            else:
                if not include_singular:
                    continue
                block = sl_pair_singular(verts, Va, Vb, n_q=singular_n_q)
            SL[np.ix_(Va, Vb)] += block

    return SL
