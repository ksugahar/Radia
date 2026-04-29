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
def extract_surface_curved(mesh, bnd_label=None, geom_order=2):
    """Extract a CURVED P2 (or higher) surface representation from a NGSolve
    mesh with mesh.Curve() already applied.

    For each BND triangle in the chosen label set, this returns the 6 P2
    Lagrange node coordinates (3 corners + 3 mid-edge nodes).  Higher
    geom_order is treated as P2 here -- the SS kernel only needs the
    surface position and its first-derivative Jacobian, and a P2 nodal
    representation is sufficient for the geometry up to second-order
    accuracy.  The Sauter-Schwab quadrature itself is reference-frame
    based and is geom-order independent.

    Args:
        mesh: NGSolve Mesh; mesh.Curve(p) MUST have been called for p>=2.
        bnd_label: optional string name of the boundary label to keep.
        geom_order: nominal mesh curving order (>=2 expected for curved
            geometry; geom_order=1 is also accepted but degenerate -- mid-
            edge nodes coincide with corner midpoints).

    Returns:
        verts (n_v, 3) corner vertex coords, indexed as in the P1 mesh.
        tris (n_t, 3) triangle corner-vertex indices into verts.
        v_global (n_v,) global vertex index in mesh.
        tri_p2_nodes (n_t, 6, 3) physical positions of the 6 P2 nodes per
            tri, ordered [v0, v1, v2, mid01, mid12, mid20] in the same
            local indexing as tris.
    """
    from ngsolve import BND, IntegrationRule

    bnd_labels = list(mesh.GetBoundaries())
    if bnd_label is not None:
        target_idx = {i for i, n in enumerate(bnd_labels) if n == bnd_label}
        if not target_idx:
            raise ValueError(
                f"BEM: BND label {bnd_label!r} not found in mesh; "
                f"available: {bnd_labels}")
    else:
        target_idx = set(range(len(bnd_labels)))

    # Reuse extract_surface for the corner topology
    verts, tris, v_global = extract_surface(mesh, bnd_label=bnd_label)
    n_t = len(tris)
    tri_p2_nodes = np.zeros((n_t, 6, 3), dtype=np.float64)

    # NGSolve T_ref = [(0,0), (1,0), (0,1)].  P2 reference nodes:
    #   k=0  -> (0, 0)        corner v0
    #   k=1  -> (1, 0)        corner v1
    #   k=2  -> (0, 1)        corner v2
    #   k=3  -> (0.5, 0)      mid-edge v0-v1
    #   k=4  -> (0.5, 0.5)    mid-edge v1-v2
    #   k=5  -> (0, 0.5)      mid-edge v2-v0
    P2_REF = np.array([
        (0.0, 0.0),
        (1.0, 0.0),
        (0.0, 1.0),
        (0.5, 0.0),
        (0.5, 0.5),
        (0.0, 0.5),
    ], dtype=np.float64)
    ir6 = IntegrationRule(points=[(p[0], p[1], 0) for p in P2_REF],
                          weights=[1.0/6.0]*6)

    # Build a global-vertex -> local-vertex mapping
    g2l = {int(g): l for l, g in enumerate(v_global)}

    # Iterate BND elements; locate each tri by its corner-vertex set
    # and write the 6 P2 nodes into tri_p2_nodes[t].
    # NGSolve does NOT guarantee that el iteration order matches our
    # local tri ordering (which came from sorted vertex iteration in
    # extract_surface).  Reconstruct the (sorted-vertex-set) -> tri-row
    # mapping so we can write into the right row.
    tris_set_to_row = {}
    for t_row in range(n_t):
        key = tuple(sorted(int(v) for v in tris[t_row]))
        tris_set_to_row[key] = t_row

    seen = 0
    for el in mesh.Elements(BND):
        if el.index not in target_idx:
            continue
        local_corners = [g2l[v.nr] for v in el.vertices]
        key = tuple(sorted(local_corners))
        t_row = tris_set_to_row.get(key)
        if t_row is None:
            continue

        # The corner ordering in tri_row may NOT match el.vertices order.
        # Build a permutation: for each ref P2 node corner k=0..2, find
        # which physical corner of THIS element it corresponds to (by
        # matching the trafo image to the corner vertex coord).
        trafo = mesh.GetTrafo(el)
        # Evaluate trafo at all 6 P2 ref points in one IR call
        physical = np.zeros((6, 3), dtype=np.float64)
        for k, ip in enumerate(ir6):
            mip = trafo(ip)
            physical[k, 0] = mip.point[0]
            physical[k, 1] = mip.point[1]
            physical[k, 2] = mip.point[2]

        # Match the 3 corners of THIS element (physical[0..2]) against
        # the 3 stored corners (verts[tris[t_row, 0..2]]).  Build perm
        # so that perm[i] = local corner index of THIS element matching
        # tris[t_row, i].
        perm = [-1, -1, -1]
        for i in range(3):
            target = verts[tris[t_row, i]]
            best = -1
            best_d = 1e30
            for j in range(3):
                d = np.linalg.norm(physical[j] - target)
                if d < best_d:
                    best_d = d; best = j
            perm[i] = best

        # Write into tri_p2_nodes:
        # corner k of OUR storage = el's corner perm[k]
        # mid-edge k+3 (between OUR corners k and (k+1)%3) = el's mid-edge
        #     between perm[k] and perm[(k+1)%3]
        # GMSH-like edge index lookup: edge {a, b} -> el ref index
        # mid-edge map: el's edge (i, j) with i<j corresponds to:
        #   (0, 1) -> ref node 3
        #   (1, 2) -> ref node 4
        #   (0, 2) -> ref node 5
        EDGE_REF_NODE = {(0, 1): 3, (1, 2): 4, (0, 2): 5}
        for k in range(3):
            tri_p2_nodes[t_row, k, :] = physical[perm[k]]
        for k in range(3):
            a = perm[k]; b = perm[(k+1) % 3]
            edge_key = (min(a, b), max(a, b))
            mid_ref = EDGE_REF_NODE[edge_key]
            tri_p2_nodes[t_row, 3 + k, :] = physical[mid_ref]
        seen += 1

    if seen != n_t:
        raise RuntimeError(
            f"extract_surface_curved: matched {seen} BND elements but "
            f"expected {n_t}")
    return verts, tris, v_global, tri_p2_nodes


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
# Phase 1.6 -- P2 Lagrange surface geometry (curved tris).
#
# Evaluate physical position r(xi, eta), tangent vectors dr/dxi, dr/deta,
# and the surface area Jacobian |dr/dxi x dr/deta| at quadrature points
# (xi, eta) on the reference triangle T_ref.  The 6 P2 nodes per tri are
# stored in tri_p2_nodes[t, k, :] with k=0..2 corners, k=3..5 mid-edges
# in the order [mid01, mid12, mid20].  P2 shape functions:
#   N_0 = (1-xi-eta) * (1 - 2*xi - 2*eta)
#   N_1 = xi * (2*xi - 1)
#   N_2 = eta * (2*eta - 1)
#   N_3 = 4 * (1-xi-eta) * xi
#   N_4 = 4 * xi * eta
#   N_5 = 4 * (1-xi-eta) * eta
# ----------------------------------------------------------------------


def _p2_basis_and_grad(xi, eta):
    """Evaluate P2 Lagrange shape functions and their (xi, eta) derivatives.

    Args:
        xi, eta: (Q,) arrays of reference coords in T_ref.

    Returns:
        N (Q, 6) shape function values
        dN_dxi (Q, 6) ∂N/∂xi
        dN_deta (Q, 6) ∂N/∂eta
    """
    L0 = 1.0 - xi - eta
    L1 = xi
    L2 = eta
    Q = xi.shape[0]
    N = np.empty((Q, 6), dtype=np.float64)
    dN_dxi = np.empty((Q, 6), dtype=np.float64)
    dN_deta = np.empty((Q, 6), dtype=np.float64)

    # Corners
    N[:, 0] = L0 * (2.0 * L0 - 1.0)
    N[:, 1] = L1 * (2.0 * L1 - 1.0)
    N[:, 2] = L2 * (2.0 * L2 - 1.0)
    # Mid-edges
    N[:, 3] = 4.0 * L0 * L1   # edge v0-v1 (xi-axis)
    N[:, 4] = 4.0 * L1 * L2   # edge v1-v2 (hypotenuse)
    N[:, 5] = 4.0 * L0 * L2   # edge v2-v0 (eta-axis)

    # ∂L0/∂xi = -1, ∂L0/∂eta = -1
    # ∂L1/∂xi = +1, ∂L1/∂eta =  0
    # ∂L2/∂xi =  0, ∂L2/∂eta = +1
    # Corners: dN0 = (4*L0 - 1) * (-1), etc.
    dN_dxi[:, 0]  = -(4.0 * L0 - 1.0)
    dN_deta[:, 0] = -(4.0 * L0 - 1.0)
    dN_dxi[:, 1]  =  (4.0 * L1 - 1.0)
    dN_deta[:, 1] = 0.0
    dN_dxi[:, 2]  = 0.0
    dN_deta[:, 2] =  (4.0 * L2 - 1.0)
    # Mid-edges
    dN_dxi[:, 3]  = 4.0 * (L0 - L1)              # 4 (∂L0/∂xi * L1 + L0 * ∂L1/∂xi) = 4(-L1 + L0)
    dN_deta[:, 3] = -4.0 * L1                    # 4 (∂L0/∂eta * L1) = -4 L1
    dN_dxi[:, 4]  = 4.0 * L2                     # 4 (∂L1/∂xi * L2) = 4 L2
    dN_deta[:, 4] = 4.0 * L1                     # 4 (L1 * ∂L2/∂eta) = 4 L1
    dN_dxi[:, 5]  = -4.0 * L2                    # 4 (∂L0/∂xi * L2) = -4 L2
    dN_deta[:, 5] = 4.0 * (L0 - L2)              # 4 (∂L0/∂eta L2 + L0 ∂L2/∂eta) = 4(-L2 + L0)
    return N, dN_dxi, dN_deta


def _eval_p2_geom(p2_nodes, xi, eta):
    """Evaluate physical position r(xi, eta), tangent vectors, and the
    surface Jacobian on a P2-curved triangle.

    Args:
        p2_nodes (6, 3) the 6 P2 node coordinates of one tri.
        xi, eta (Q,) arrays of reference coords.

    Returns:
        r       (Q, 3) physical positions
        dr_dxi  (Q, 3) tangent ∂r/∂xi
        dr_deta (Q, 3) tangent ∂r/∂eta
        J       (Q,)   surface element |dr/dxi x dr/deta|
    """
    N, dN_dxi, dN_deta = _p2_basis_and_grad(xi, eta)
    r       = N       @ p2_nodes        # (Q, 6) @ (6, 3) -> (Q, 3)
    dr_dxi  = dN_dxi  @ p2_nodes
    dr_deta = dN_deta @ p2_nodes
    cross = np.cross(dr_dxi, dr_deta)
    J = np.sqrt(np.sum(cross * cross, axis=1))
    return r, dr_dxi, dr_deta, J


def _eval_p2_geom_with_normal(p2_nodes, xi, eta):
    """Same as _eval_p2_geom but also returns the unit outward normal at
    each quadrature point (constructed from the cross product, sign
    inherited from the (v0, v1, v2) winding order)."""
    N, dN_dxi, dN_deta = _p2_basis_and_grad(xi, eta)
    r       = N       @ p2_nodes
    dr_dxi  = dN_dxi  @ p2_nodes
    dr_deta = dN_deta @ p2_nodes
    cross = np.cross(dr_dxi, dr_deta)
    J = np.sqrt(np.sum(cross * cross, axis=1))
    n_hat = cross / (J[:, None] + 1e-300)
    return r, dr_dxi, dr_deta, J, n_hat


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


def _ss_block_curved(p2a, p2b, case, n_q):
    """Curved (P2 geometry) Galerkin SL block in PERMUTED indexing.

    Same NGSolve-canonical permutation as _ss_block:
      common_vertex: shared vertex at local 0 of both
      common_edge:   shared edge from local 0 to local 1
      identical:     T_a == T_b

    Both p2a and p2b are (6, 3) P2 node arrays in the canonical
    permutation: corners 0..2 in canonical order, mid-edges 3..5 in
    [mid01, mid12, mid20].

    Returns (3, 3) block where rows are P1 hat L_i^a (canonical local)
    and cols are L_j^b (canonical local).
    """
    xi_a, eta_a, xi_b, eta_b, w_q = _ss_quad_pts(case, n_q)

    # Curved physical positions and per-quadpt Jacobians
    r_a, _, _, J_a = _eval_p2_geom(p2a, xi_a, eta_a)
    r_b, _, _, J_b = _eval_p2_geom(p2b, xi_b, eta_b)

    diff = r_a - r_b
    r_dist = np.sqrt(np.sum(diff * diff, axis=1))
    r_dist = np.where(r_dist == 0, 1e-300, r_dist)

    # P1 hat values in reference coords (geom-order independent)
    L0_a = 1.0 - xi_a - eta_a
    L1_a = xi_a
    L2_a = eta_a
    L0_b = 1.0 - xi_b - eta_b
    L1_b = xi_b
    L2_b = eta_b
    La = np.stack([L0_a, L1_a, L2_a], axis=1)
    Lb = np.stack([L0_b, L1_b, L2_b], axis=1)

    # Galerkin entry: 1/(4 pi r) * J_a(xi_a) * J_b(xi_b) * w_q
    # (Sauter-Schwab w_q already absorbs the singular Jacobian factor;
    # the per-quadpt J_a, J_b replace the constant 4*A_a*A_b of the flat
    # version.)
    weight = (J_a * J_b) / (4.0 * np.pi * r_dist) * w_q
    return La.T @ (weight[:, None] * Lb)


def _ss_block_dl_curved(p2a, p2b, case, n_q):
    """Curved (P2 geometry) Galerkin DL block in PERMUTED indexing.

    DL kernel ∂G/∂n_y = +(r-r')·n_y / (4*pi*|r-r'|^3) with n_y the unit
    outward normal of T_b at r' (= curved-tri normal at quadpt).
    """
    xi_a, eta_a, xi_b, eta_b, w_q = _ss_quad_pts(case, n_q)

    r_a, _, _, J_a            = _eval_p2_geom(p2a, xi_a, eta_a)
    r_b, _, _, J_b, n_hat_b   = _eval_p2_geom_with_normal(p2b, xi_b, eta_b)

    diff = r_a - r_b
    r2 = np.sum(diff * diff, axis=1)
    r_dist = np.sqrt(r2)
    r_dist = np.where(r_dist == 0, 1e-300, r_dist)
    dot_n = np.sum(diff * n_hat_b, axis=1)   # (Q,)

    L0_a = 1.0 - xi_a - eta_a
    L1_a = xi_a
    L2_a = eta_a
    L0_b = 1.0 - xi_b - eta_b
    L1_b = xi_b
    L2_b = eta_b
    La = np.stack([L0_a, L1_a, L2_a], axis=1)
    Lb = np.stack([L0_b, L1_b, L2_b], axis=1)

    weight = dot_n / (4.0 * np.pi * r2 * r_dist) * J_a * J_b * w_q
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


def _permute_p2_corners(p2_nodes, perm):
    """Permute the 6 P2 nodes so corners follow `perm` (a 3-tuple of
    new corner indices in old ordering).  Mid-edge order follows the
    canonical [mid01, mid12, mid20] in the new local indexing.

    Args:
        p2_nodes (6, 3) original P2 nodes [v0, v1, v2, mid01, mid12, mid20]
        perm     (3,)   new_corner_k = old_corner[perm[k]]

    Returns:
        (6, 3) permuted P2 nodes.
    """
    out = np.zeros_like(p2_nodes)
    for k in range(3):
        out[k] = p2_nodes[perm[k]]
    # mid-edge mapping: old mid-edge between old corners (i, j) stored as
    # k=3 (i=0,j=1), k=4 (i=1,j=2), k=5 (i=0,j=2 NOTE: stored as v2-v0 wrap).
    # Actually our convention: idx 3 = mid01, 4 = mid12, 5 = mid20.
    # mid between old corners (a, b) (a<b):
    #   (0,1) -> 3, (1,2) -> 4, (0,2) -> 5
    OLD_MID = {(0, 1): 3, (1, 2): 4, (0, 2): 5}
    for k in range(3):
        a = perm[k]; b = perm[(k+1) % 3]
        edge = (min(a, b), max(a, b))
        out[3 + k] = p2_nodes[OLD_MID[edge]]
    return out


def sl_pair_singular_curved(p2_a, p2_b, tri_a, tri_b, n_q=8):
    """SL block for any singular pair using P2-curved geometry.

    Returns (3, 3) block in tri_a (rows) x tri_b (cols) ORIGINAL local
    indexing (matching the corner ordering passed in via tri_a/tri_b).
    """
    cls, _ = share_class(tri_a, tri_b)
    shared_pairs = [(la, lb) for la in range(3) for lb in range(3)
                    if int(tri_a[la]) == int(tri_b[lb])]

    if cls == 'vertex':
        la, lb = shared_pairs[0]
        pa = (la, (la + 1) % 3, (la + 2) % 3)
        pb = (lb, (lb + 1) % 3, (lb + 2) % 3)
        p2a = _permute_p2_corners(p2_a, pa)
        p2b = _permute_p2_corners(p2_b, pb)
        M_perm = _ss_block_curved(p2a, p2b, "common_vertex", n_q)
        M_orig = np.zeros((3, 3))
        M_orig[np.ix_(pa, pb)] = M_perm
        return M_orig
    if cls == 'edge':
        (la0, lb0), (la1, lb1) = shared_pairs
        apex_a = ({0, 1, 2} - {la0, la1}).pop()
        apex_b = ({0, 1, 2} - {lb0, lb1}).pop()
        pa = (la0, la1, apex_a)
        pb = (lb0, lb1, apex_b)
        p2a = _permute_p2_corners(p2_a, pa)
        p2b = _permute_p2_corners(p2_b, pb)
        M_perm = _ss_block_curved(p2a, p2b, "common_edge", n_q)
        M_orig = np.zeros((3, 3))
        M_orig[np.ix_(pa, pb)] = M_perm
        return M_orig
    if cls == 'identical':
        # Resolve perm: vb_in_va_order[la] = vb[perm_b[la]] == va[la]
        perm_b = [None, None, None]
        for la, lb in shared_pairs:
            perm_b[la] = lb
        p2b_in_va = _permute_p2_corners(p2_b, tuple(perm_b))
        M_va = _ss_block_curved(p2_a, p2b_in_va, "identical", n_q)
        # Permute cols back to vb original indexing
        M_orig = np.zeros((3, 3))
        for la in range(3):
            M_orig[:, perm_b[la]] = M_va[:, la]
        return M_orig
    raise ValueError(f"sl_pair_singular_curved: non-singular cls={cls}")


def _perm_parity_sign(perm):
    """+1 if perm is a cyclic (= even) permutation of (0, 1, 2),
    -1 if odd (one swap from cyclic).  For triangle permutations,
    even = winding preserved (= natural outward normal preserved).
    """
    a, b, c = int(perm[0]), int(perm[1]), int(perm[2])
    cyclic_orderings = {(0, 1, 2), (1, 2, 0), (2, 0, 1)}
    return +1 if (a, b, c) in cyclic_orderings else -1


def dl_pair_singular_curved(p2_a, p2_b, tri_a, tri_b, n_q=8):
    """DL block for any singular pair using P2-curved geometry.

    NOTE on orientation: when we permute corners to bring the SS rule's
    canonical "shared vertex/edge at local 0/0-1" form, the resulting
    canonical winding may differ from the element's NATURAL winding.
    For SL (kernel = 1/r) this is irrelevant.  For DL (kernel
    (r-r')*n_y/r^3) the n_y direction depends on the cross-product of
    edge vectors, which flips sign if the canonical permutation is a
    non-cyclic permutation of (0, 1, 2).  We detect non-cyclic perm
    via _perm_parity_sign and negate the result accordingly so n_y
    always points outward.
    """
    cls, _ = share_class(tri_a, tri_b)
    shared_pairs = [(la, lb) for la in range(3) for lb in range(3)
                    if int(tri_a[la]) == int(tri_b[lb])]

    if cls == 'vertex':
        la, lb = shared_pairs[0]
        pa = (la, (la + 1) % 3, (la + 2) % 3)
        pb = (lb, (lb + 1) % 3, (lb + 2) % 3)
        # pa, pb are always cyclic perms (constructed via cyclic shift),
        # so winding is preserved and no sign correction needed.
        p2a = _permute_p2_corners(p2_a, pa)
        p2b = _permute_p2_corners(p2_b, pb)
        M_perm = _ss_block_dl_curved(p2a, p2b, "common_vertex", n_q)
        # Sign: pa and pb both cyclic -> normals point outward -> +1
        M_orig = np.zeros((3, 3))
        M_orig[np.ix_(pa, pb)] = M_perm
        return M_orig
    if cls == 'edge':
        (la0, lb0), (la1, lb1) = shared_pairs
        apex_a = ({0, 1, 2} - {la0, la1}).pop()
        apex_b = ({0, 1, 2} - {lb0, lb1}).pop()
        pa = (la0, la1, apex_a)
        pb = (lb0, lb1, apex_b)
        # If pa or pb is non-cyclic, the canonical T_a/T_b winding is
        # reversed from natural (= cross gives INWARD normal).  Flip the
        # contribution sign for each non-cyclic perm.
        sign_a = _perm_parity_sign(pa)
        sign_b = _perm_parity_sign(pb)
        p2a = _permute_p2_corners(p2_a, pa)
        p2b = _permute_p2_corners(p2_b, pb)
        M_perm = _ss_block_dl_curved(p2a, p2b, "common_edge", n_q)
        # n_b only enters via dot_n in the kernel -> sign of n_b flips
        # the kernel sign.  pa winding doesn't enter the kernel
        # (dr_a/dxi etc. only set positions, not n).  So only sign_b
        # affects DL.
        M_perm *= sign_b
        M_orig = np.zeros((3, 3))
        M_orig[np.ix_(pa, pb)] = M_perm
        return M_orig
    if cls == 'identical':
        perm_b = [None, None, None]
        for la, lb in shared_pairs:
            perm_b[la] = lb
        sign_b = _perm_parity_sign(tuple(perm_b))
        p2b_in_va = _permute_p2_corners(p2_b, tuple(perm_b))
        M_va = _ss_block_dl_curved(p2_a, p2b_in_va, "identical", n_q)
        M_va *= sign_b   # flip if T_b winding got reversed
        M_orig = np.zeros((3, 3))
        for la in range(3):
            M_orig[:, perm_b[la]] = M_va[:, la]
        return M_orig
    raise ValueError(f"dl_pair_singular_curved: non-singular cls={cls}")


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


# ----------------------------------------------------------------------
# Phase 1.8 (REJECTED) -- mesh.GetTrafo-based geometry experiment.
#
# Hypothesis was that NGSolve's mesh.Curve(2) uses higher-order curving
# than P2 Lagrange, and Phase 1.6's P2 Lagrange interpolation would
# leave residual geometry error.
#
# Empirically REFUTED on N=99 sphere: Phase 1.6 SL matches NGSolve.bem
# (bonus_intorder=10 anchor) to rel 4.3e-10 already.  The earlier
# "0.1-0.2 geometry mismatch" diagnostic was a buggy comparison test
# that confused corner permutations between NGSolve's reference frame
# (corner 0 at (1,0)) and our local frame (corner 0 at (0,0)).
#
# Phase 1.8 trafo path has been tested and is BOTH 70x slower (per-
# quadpt Python trafo calls) AND less accurate (rel 0.24 vs reference)
# due to permutation logic complexity.  Code below is kept as a
# reference for future debugging only; it is NOT used by any public
# entry point.
# ----------------------------------------------------------------------


def _trafo_eval(mesh, el_index_in_bnd, xi, eta):
    """Evaluate physical position, dr/dxi, dr/deta, surface Jacobian, and
    outward unit normal at points (xi, eta) on a single BND element.

    Args:
        mesh: NGSolve Mesh.
        el_index_in_bnd: int, index in the BND iteration order (= ngsolve
            ElementId(BND, idx) used by mesh.GetTrafo).  Caller is
            responsible for ensuring this matches their tri row.
        xi, eta: (Q,) arrays.

    Returns:
        r       (Q, 3) physical positions
        J       (Q,)   surface element |dr/dxi x dr/deta|
        n_hat   (Q, 3) unit outward normal (cross / J)
    """
    from ngsolve import IntegrationRule, ElementId, BND
    Q = xi.shape[0]
    pts_with_w = [(xi[k], eta[k], 0) for k in range(Q)]
    ws = [1.0] * Q  # weights are dummy; we extract pos/Jac/normal only
    ir = IntegrationRule(points=pts_with_w, weights=ws)
    eid = ElementId(BND, el_index_in_bnd)
    trafo = mesh.GetTrafo(eid)

    r = np.zeros((Q, 3), dtype=np.float64)
    J = np.zeros((Q,), dtype=np.float64)
    n_hat = np.zeros((Q, 3), dtype=np.float64)
    for k, ip in enumerate(ir):
        mip = trafo(ip)
        r[k, 0] = mip.point[0]; r[k, 1] = mip.point[1]; r[k, 2] = mip.point[2]
        J[k] = mip.measure
        jac = np.asarray(mip.jacobi)   # shape (3, 2): ∂r/∂xi, ∂r/∂eta
        cross = np.cross(jac[:, 0], jac[:, 1])
        nrm = np.linalg.norm(cross) + 1e-300
        n_hat[k] = cross / nrm
    return r, J, n_hat


def _ss_block_curved_trafo(mesh, bnd_el_a, bnd_el_b, case, n_q,
                            perm_a, perm_b):
    """Curved Galerkin SL block evaluated via mesh.GetTrafo at SS quad pts.

    Args:
        mesh: NGSolve Mesh.
        bnd_el_a, bnd_el_b: BND element indices (in mesh.Elements(BND) order).
        case: "common_vertex" / "common_edge" / "identical".
        n_q: 1D Gauss-Legendre order.
        perm_a, perm_b: 3-tuples mapping canonical-perm corner k to
            ORIGINAL local corner index.  We need to inverse-permute the
            SS quad reference coords (xi, eta) before passing to the
            mesh trafo.

    Returns:
        (3, 3) SL block in ORIGINAL local indexing of bnd_el_a/_b.
    """
    xi_a_canon, eta_a_canon, xi_b_canon, eta_b_canon, w_q = _ss_quad_pts(case, n_q)

    # The SS quad coords are in CANONICAL permutation (shared corner at local 0).
    # The trafo evaluates in ORIGINAL element ordering.  We need to map
    # canonical (xi, eta) -> original (xi, eta) using the inverse of perm_a/_b.
    # Coords transformation: new_corner_k = old_corner[perm[k]].  In barycentric:
    #   L_canon = (L_canon_0, L_canon_1, L_canon_2)
    #   L_orig[perm[k]] = L_canon[k]
    # i.e. L_orig[perm[0]] = L_canon[0] = 1-xi-eta, L_orig[perm[1]] = xi, L_orig[perm[2]] = eta.
    # Since perm is a permutation of (0,1,2), inverse perm sends each L back.
    # In (xi, eta) coords: L_0 = 1-xi-eta, L_1 = xi, L_2 = eta.
    # We want L_orig = inverse_perm(L_canon).  Then xi_orig = L_orig[1], eta_orig = L_orig[2].

    def map_canonical_to_orig(xi_c, eta_c, perm):
        # L_canon = (1 - xi_c - eta_c, xi_c, eta_c)
        L_canon = np.stack([1.0 - xi_c - eta_c, xi_c, eta_c], axis=0)
        # L_orig[perm[k]] = L_canon[k]
        L_orig = np.zeros_like(L_canon)
        for k in range(3):
            L_orig[perm[k]] = L_canon[k]
        xi_o = L_orig[1]
        eta_o = L_orig[2]
        return xi_o, eta_o

    xi_a, eta_a = map_canonical_to_orig(xi_a_canon, eta_a_canon, perm_a)
    xi_b, eta_b = map_canonical_to_orig(xi_b_canon, eta_b_canon, perm_b)

    r_a, J_a, _      = _trafo_eval(mesh, bnd_el_a, xi_a, eta_a)
    r_b, J_b, n_hat_b = _trafo_eval(mesh, bnd_el_b, xi_b, eta_b)

    diff = r_a - r_b
    r2 = np.sum(diff * diff, axis=1)
    r_dist = np.sqrt(r2)
    r_dist = np.where(r_dist == 0, 1e-300, r_dist)

    # Hat values: in original local indexing, L_orig_i evaluated at the
    # ORIGINAL (xi, eta) coords.  We have those.
    L0_a = 1.0 - xi_a - eta_a
    L1_a = xi_a
    L2_a = eta_a
    L0_b = 1.0 - xi_b - eta_b
    L1_b = xi_b
    L2_b = eta_b
    La = np.stack([L0_a, L1_a, L2_a], axis=1)
    Lb = np.stack([L0_b, L1_b, L2_b], axis=1)

    weight = (J_a * J_b) / (4.0 * np.pi * r_dist) * w_q
    return La.T @ (weight[:, None] * Lb)


def _build_bnd_to_row_map(mesh, verts, tris, v_global, bnd_label):
    """Build the bnd_to_row map: index in mesh.Elements(BND) iteration ->
    row index in our `tris` array.  Used by trafo-based assembly to look
    up each BND element's row.
    """
    from ngsolve import BND
    g2l = {int(g): l for l, g in enumerate(v_global)}
    bnd_labels = list(mesh.GetBoundaries())
    if bnd_label is not None:
        target_idx = {i for i, n in enumerate(bnd_labels) if n == bnd_label}
    else:
        target_idx = set(range(len(bnd_labels)))
    tri_set_to_row = {}
    for r in range(len(tris)):
        key = tuple(sorted(int(c) for c in tris[r]))
        tri_set_to_row[key] = r

    bnd_to_row = []      # bnd iter index -> our row (or -1)
    bnd_corner_perm = []  # bnd iter idx -> 3-tuple perm such that
                          # verts[tris[r, perm[k]]] == el.vertices[k]
    for el in mesh.Elements(BND):
        if el.index not in target_idx:
            bnd_to_row.append(-1)
            bnd_corner_perm.append(None)
            continue
        local_corners = [g2l[v.nr] for v in el.vertices]
        key = tuple(sorted(local_corners))
        r = tri_set_to_row.get(key, -1)
        bnd_to_row.append(r)
        if r >= 0:
            # Build perm: perm[k] = local index in tris[r] of el.vertices[k]
            tris_row = list(int(c) for c in tris[r])
            perm = [tris_row.index(c) for c in local_corners]
            bnd_corner_perm.append(tuple(perm))
        else:
            bnd_corner_perm.append(None)
    return bnd_to_row, bnd_corner_perm


def _assemble_SL_dense_trafo(mesh, verts, tris, v_global, *,
                              bnd_label=None,
                              regular_quad_degree=11,
                              singular_n_q=8):
    """Build dense Galerkin SL using mesh.GetTrafo for curved geometry.

    Slower than P2 Lagrange (multiple trafo calls per pair) but uses
    exactly the same curving as NGSolve.bem.

    Returns:
        SL (n_v, n_v) dense matrix, indexed by surface-local vertex
        (matches extract_surface conventions).
    """
    from ngsolve import BND
    n_v = len(verts)
    n_t = len(tris)
    SL = np.zeros((n_v, n_v), dtype=np.float64)
    bnd_to_row, bnd_corner_perm = _build_bnd_to_row_map(
        mesh, verts, tris, v_global, bnd_label)
    # Inverse: row -> bnd iter index
    row_to_bnd = [-1] * n_t
    for bi, r in enumerate(bnd_to_row):
        if r >= 0:
            row_to_bnd[r] = bi

    # Pre-cache regular Gauss reference quadrature points.
    q = tri_quad(regular_quad_degree)
    qxi  = q[:, 0].astype(np.float64)
    qeta = q[:, 1].astype(np.float64)
    qw_ref = q[:, 2].astype(np.float64)

    # Pre-evaluate per-tri regular-pair physical positions and Jacobians
    # via trafo (one call per tri, n_q points each).
    pts_cache = []
    w_cache = []
    L_cache = []
    for t in range(n_t):
        bi = row_to_bnd[t]
        # Apply inverse perm so we evaluate at the ORIGINAL element local
        # frame (where mesh.GetTrafo expects).  Our `tris[t, k]` is in
        # surface-local order; we need to know which local index in the
        # NGSolve element corresponds to our local k.
        # bnd_corner_perm[bi][k] = our_local_idx of el.vertices[k].
        # So inv_perm[our_k] = el_k.  We want xi, eta in EL frame given
        # they're set in OUR frame (with corners 0,1,2 of tris[t]).
        perm_el_to_our = bnd_corner_perm[bi]
        perm_our_to_el = [perm_el_to_our.index(k) for k in range(3)]
        # qxi, qeta are in OUR frame.  Convert to EL frame: L_our_i evaluated
        # at our (xi, eta).  L_el_k = L_our_{perm_el_to_our[k]}.
        # Then xi_el = L_el_1, eta_el = L_el_2.
        L_our = np.stack([1.0 - qxi - qeta, qxi, qeta], axis=0)
        L_el = np.zeros_like(L_our)
        for k in range(3):
            L_el[k] = L_our[perm_el_to_our[k]]
        xi_el = L_el[1]
        eta_el = L_el[2]
        r_t, J_t, _ = _trafo_eval(mesh, bi, xi_el, eta_el)
        pts_cache.append(r_t)
        w_cache.append(qw_ref * J_t)
        L_cache.append(np.stack([1.0 - qxi - qeta, qxi, qeta], axis=1))

    for a in range(n_t):
        Va = tris[a]
        bi_a = row_to_bnd[a]
        for b in range(n_t):
            Vb = tris[b]
            bi_b = row_to_bnd[b]
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
                # Build canonical-perm permutations for SS rule.
                shared_pairs = [(la, lb) for la in range(3) for lb in range(3)
                                if int(Va[la]) == int(Vb[lb])]
                if cls == 'vertex':
                    la, lb = shared_pairs[0]
                    perm_a_canon = (la, (la + 1) % 3, (la + 2) % 3)
                    perm_b_canon = (lb, (lb + 1) % 3, (lb + 2) % 3)
                elif cls == 'edge':
                    (la0, lb0), (la1, lb1) = shared_pairs
                    apex_a = ({0, 1, 2} - {la0, la1}).pop()
                    apex_b = ({0, 1, 2} - {lb0, lb1}).pop()
                    perm_a_canon = (la0, la1, apex_a)
                    perm_b_canon = (lb0, lb1, apex_b)
                else:  # identical
                    perm_a_canon = (0, 1, 2)
                    perm_b_canon = [None, None, None]
                    for la, lb in shared_pairs:
                        perm_b_canon[la] = lb
                    perm_b_canon = tuple(perm_b_canon)

                # Compose: canonical -> our_local -> el_local
                # canonical perm_a maps canonical k -> our local
                # our_local -> el_local via bnd_corner_perm (which is el_k -> our_k)
                # So total perm: canonical k -> el_k where:
                #   our_local = perm_a_canon[k]
                #   el_local = bnd_corner_perm[bi_a].index(our_local)
                perm_a_full = tuple(
                    bnd_corner_perm[bi_a].index(perm_a_canon[k]) for k in range(3))
                perm_b_full = tuple(
                    bnd_corner_perm[bi_b].index(perm_b_canon[k]) for k in range(3))

                cls_quad = {"vertex": "common_vertex",
                             "edge": "common_edge",
                             "identical": "identical"}[cls]
                block_perm = _ss_block_curved_trafo(
                    mesh, bi_a, bi_b, cls_quad, singular_n_q,
                    perm_a_full, perm_b_full)
                # block_perm is in EL local indexing (rows: el_a, cols: el_b).
                # Map back to OUR local (rows: our_a, cols: our_b).
                # bnd_corner_perm[bi_a][k] = our_local_idx of el.vertices[k].
                # block_perm[el_i, el_j] in el local.  We want
                #   block_our[our_i, our_j] = block_perm[el_i, el_j] where
                #   el_i = bnd_corner_perm[bi_a].index(our_i), etc.
                block = np.zeros((3, 3))
                for our_i in range(3):
                    el_i = bnd_corner_perm[bi_a].index(our_i)
                    for our_j in range(3):
                        el_j = bnd_corner_perm[bi_b].index(our_j)
                        block[our_i, our_j] = block_perm[el_i, el_j]
            SL[np.ix_(Va, Vb)] += block

    return SL


def assemble_bem_matrices_trafo(mesh, *, bnd_label=None,
                                  regular_quad_degree=11, singular_n_q=8):
    """Phase 1.8 entry point: build SL, M, K via mesh.GetTrafo (no P2 Lagrange).

    DL is omitted in this initial Phase 1.8 (its diagonal-jump convention
    needs separate reconciliation; the SIBC formula folds the (1/2)M
    jump into the system matrix anyway).

    Returns dict with keys SL, M, K, verts, tris, v_global.
    """
    verts, tris, v_global = extract_surface(mesh, bnd_label=bnd_label)
    SL = _assemble_SL_dense_trafo(mesh, verts, tris, v_global,
                                    bnd_label=bnd_label,
                                    regular_quad_degree=regular_quad_degree,
                                    singular_n_q=singular_n_q)

    # M and K: closed-form per-tri 3x3 blocks (P1 hat).  Use the
    # CURVED Jacobian via trafo so M and K match NGSolve exactly on
    # curved meshes.
    n_v = len(verts)
    n_t = len(tris)
    M = np.zeros((n_v, n_v))
    K = np.zeros((n_v, n_v))
    bnd_to_row, bnd_corner_perm = _build_bnd_to_row_map(
        mesh, verts, tris, v_global, bnd_label)
    row_to_bnd = [-1] * n_t
    for bi, r in enumerate(bnd_to_row):
        if r >= 0:
            row_to_bnd[r] = bi

    # Use a moderate quadrature for M, K; their integrand is polynomial
    # times surface-element factor.
    q = tri_quad(regular_quad_degree)
    qxi  = q[:, 0]; qeta = q[:, 1]; qw_ref = q[:, 2]
    # Hat values on T_ref
    L_ref = np.stack([1.0 - qxi - qeta, qxi, qeta], axis=1)   # (n_q, 3)

    for t in range(n_t):
        bi = row_to_bnd[t]
        perm_el_to_our = bnd_corner_perm[bi]
        L_our = np.stack([1.0 - qxi - qeta, qxi, qeta], axis=0)
        L_el = np.zeros_like(L_our)
        for k in range(3):
            L_el[k] = L_our[perm_el_to_our[k]]
        xi_el = L_el[1]; eta_el = L_el[2]
        # r evaluated at trafo, plus Jacobian
        r_t, J_t, _ = _trafo_eval(mesh, bi, xi_el, eta_el)
        # Surface element factor: w_q_phys = qw_ref * J_t
        w_q_phys = qw_ref * J_t
        # M_ij = sum_q L_ref_i(q) * L_ref_j(q) * w_q_phys
        M_loc = L_ref.T @ (w_q_phys[:, None] * L_ref)
        # For K (Laplace-Beltrami): use the cotangent formula on the
        # CURVED triangle (approx by flat tri at corners).  For Phase 1.8
        # production we use the same formula as Phase 1.5.  TODO: use
        # surface-grad of P1 hats in (xi, eta) and evaluate on curved
        # geometry for a more accurate K.
        v0, v1, v2 = verts[tris[t]]
        e0_v = v2 - v1; e1_v = v0 - v2; e2_v = v1 - v0
        cross = np.cross(v1 - v0, v2 - v0)
        A_flat = 0.5 * np.linalg.norm(cross)
        K_loc = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                ei = (e0_v, e1_v, e2_v)[i]
                ej = (e0_v, e1_v, e2_v)[j]
                K_loc[i, j] = (ei @ ej) / (4.0 * A_flat)
        idx = tris[t]
        M[np.ix_(idx, idx)] += M_loc
        K[np.ix_(idx, idx)] += K_loc

    return {"SL": SL, "M": M, "K": K,
            "verts": verts, "tris": tris, "v_global": v_global}


def assemble_SL_dense_curved(verts, tris, tri_p2_nodes, *,
                                regular_quad_degree=11,
                                singular_n_q=8):
    """Build dense Galerkin Laplace SL on a curved-P2 surface mesh.

    Args:
        verts (n_v, 3): corner vertex coords.
        tris  (n_t, 3): triangle corner-vertex indices.
        tri_p2_nodes (n_t, 6, 3): the 6 P2 node coords per tri (corners
            + mid-edges).  See extract_surface_curved.
        regular_quad_degree, singular_n_q: as for assemble_SL_dense.

    Returns (n_v, n_v) float64 SL matrix.
    """
    n_v = len(verts)
    n_t = len(tris)
    SL = np.zeros((n_v, n_v), dtype=np.float64)

    # Per-tri pre-computed regular quadrature: (xi, eta, w_ref) on T_ref,
    # plus their physical-position evaluations on the curved geometry.
    q = tri_quad(regular_quad_degree)
    qxi  = q[:, 0].astype(np.float64)
    qeta = q[:, 1].astype(np.float64)
    qw_ref = q[:, 2].astype(np.float64)

    pts_cache = []   # (n_q, 3) physical positions per tri
    w_cache = []     # (n_q,)   physical weights = w_ref * J(xi, eta)
    L_cache = []     # (n_q, 3) hat values (geom-order independent)
    for t in range(n_t):
        r_t, _, _, J_t = _eval_p2_geom(tri_p2_nodes[t], qxi, qeta)
        pts_cache.append(r_t)
        w_cache.append(qw_ref * J_t)
        L_cache.append(np.stack([1.0 - qxi - qeta, qxi, qeta], axis=1))

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
                block = sl_pair_singular_curved(
                    tri_p2_nodes[a], tri_p2_nodes[b],
                    Va, Vb, n_q=singular_n_q)
            SL[np.ix_(Va, Vb)] += block

    return SL


def assemble_DL_dense_curved(verts, tris, tri_p2_nodes, *,
                               regular_quad_degree=11,
                               singular_n_q=8):
    """Build dense Galerkin Laplace DL on a curved-P2 surface mesh.

    Convention matches NGSolve.bem LaplaceDL with outward normal of T_b.
    For the diagonal (identical) pair on a CURVED triangle, n_y is no
    longer constant, so (r-r')·n_y does NOT vanish identically -- the
    contribution is small (O(curvature)) but non-zero.  We compute it
    via Sauter-Schwab as for any other singular pair.
    """
    n_v = len(verts)
    n_t = len(tris)
    DL = np.zeros((n_v, n_v), dtype=np.float64)

    q = tri_quad(regular_quad_degree)
    qxi  = q[:, 0].astype(np.float64)
    qeta = q[:, 1].astype(np.float64)
    qw_ref = q[:, 2].astype(np.float64)

    pts_cache = []
    w_cache = []
    L_cache = []
    n_cache = []     # (n_q, 3) per-quadpt unit normal per tri
    for t in range(n_t):
        r_t, _, _, J_t, n_t_arr = _eval_p2_geom_with_normal(
            tri_p2_nodes[t], qxi, qeta)
        pts_cache.append(r_t)
        w_cache.append(qw_ref * J_t)
        L_cache.append(np.stack([1.0 - qxi - qeta, qxi, qeta], axis=1))
        n_cache.append(n_t_arr)

    for a in range(n_t):
        Va = tris[a]
        for b in range(n_t):
            Vb = tris[b]
            cls, _ = share_class(Va, Vb)
            if cls == 'regular':
                pts_a, w_a, La = pts_cache[a], w_cache[a], L_cache[a]
                pts_b, w_b, Lb = pts_cache[b], w_cache[b], L_cache[b]
                n_b_q = n_cache[b]   # (n_qb, 3)
                dx = pts_a[:, None, :] - pts_b[None, :, :]
                r2 = np.sum(dx * dx, axis=-1)
                r = np.sqrt(r2)
                # ∂G/∂n_y per quad pair: dx · n_b[b_q] / (4 pi r^3)
                # Use einsum to broadcast over (a_q, b_q, 3)
                dot = np.einsum('abc,bc->ab', dx, n_b_q)
                K = dot / (4.0 * np.pi * r2 * r)
                Tmat = (K * w_b[None, :]) @ Lb
                block = La.T @ (w_a[:, None] * Tmat)
            else:
                block = dl_pair_singular_curved(
                    tri_p2_nodes[a], tri_p2_nodes[b],
                    Va, Vb, n_q=singular_n_q)
            DL[np.ix_(Va, Vb)] += block

    return DL


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
