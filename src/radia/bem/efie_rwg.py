"""Intree BEM-A Galerkin assembler for the Laplace single-layer operator
on RWG (HDivSurface order=0) basis on flat triangular surface meshes.

Companion to ``radia.bem.sibc_hacapk`` (scalar SS Galerkin on H1 basis).
This module provides the vector / surface-current analog needed for
coil self-inductance via the Weggler stabilized EFIE saddle-point::

    [ SL   D^T ] [ J ]   [ 0 ]
    [ D    0   ] [ p ] = [ g ]

The Galerkin entry is::

    L_ij = (1 / 4 pi) * int_{T_a} int_{T_b}  (J_i(r) . J_j(r')) / |r - r'|
           dS dS'

The Sauter-Schwab Duffy 4D quadrature points are imported from the
**already validated** scalar pipeline in ``radia.bem.sibc_hacapk``
(Phase 1.2-1.5, agrees to 1e-10 against ``ngsolve.bem.LaplaceSL`` on
H1).  Only the basis evaluation differs.

RWG orientation is derived **automatically** from global vertex
indices (sigma_k = +1 if the CCW traversal of local edge k goes from
the smaller to the larger global vertex ID).  Permutation of corners
to the SS canonical form does NOT require manual sign tracking.
"""
from __future__ import annotations

import math
import numpy as np

from radia.bem.sibc_hacapk import _ss_quad_pts as _scalar_ss_quad_pts

INV_4PI = 1.0 / (4.0 * math.pi)


# ----------------------------------------------------------------------
# Edge / DOF topology
# ----------------------------------------------------------------------
# Local edge indexing: edge k is opposite corner k.
#   edge 0 spans corners (1, 2)
#   edge 1 spans corners (2, 0)
#   edge 2 spans corners (0, 1)
# CCW traversal of edge k goes from corner (k+1)%3 to corner (k+2)%3.
LOCAL_EDGE_ENDPOINTS = ((1, 2), (2, 0), (0, 1))


def build_edge_topology(tris):
    """Build the edge-DOF data structures used by RWG assembly.

    Args:
        tris: (n_t, 3) int array of triangle vertex indices (CCW).

    Returns:
        edges: (n_e, 2) int array of (vmin, vmax) per global edge.
        tri_edges: (n_t, 3) int array, tri_edges[t, k] = global edge
            index for local edge k of triangle t.
        edge_tris: (n_e, 2) int array, two adjacent triangles per edge
            (-1 in column 1 if boundary).
        edge_tri_local: (n_e, 2) int array, local edge index k inside
            each adjacent triangle.
    """
    tris = np.asarray(tris, dtype=np.int64)
    n_t = len(tris)
    edge_dict = {}
    tri_edges = np.zeros((n_t, 3), dtype=np.int64)
    for t in range(n_t):
        v = tris[t]
        for k, (a, b) in enumerate(LOCAL_EDGE_ENDPOINTS):
            va, vb = int(v[a]), int(v[b])
            key = (min(va, vb), max(va, vb))
            if key not in edge_dict:
                edge_dict[key] = len(edge_dict)
            tri_edges[t, k] = edge_dict[key]

    n_e = len(edge_dict)
    edges = np.zeros((n_e, 2), dtype=np.int64)
    for key, idx in edge_dict.items():
        edges[idx] = key

    edge_tris = -np.ones((n_e, 2), dtype=np.int64)
    edge_tri_local = -np.ones((n_e, 2), dtype=np.int64)
    counts = np.zeros(n_e, dtype=np.int64)
    for t in range(n_t):
        for k in range(3):
            e = tri_edges[t, k]
            c = counts[e]
            if c < 2:
                edge_tris[e, c] = t
                edge_tri_local[e, c] = k
                counts[e] = c + 1
    return edges, tri_edges, edge_tris, edge_tri_local


# ----------------------------------------------------------------------
# Triangle reference quadrature  (degree-5 hard-coded; higher via Duffy)
# ----------------------------------------------------------------------
_TRI_DEG5 = np.array([
    [0.33333333333333333, 0.33333333333333333, 0.22500000000000000 / 2],
    [0.79742698535308720, 0.10128650732345633, 0.12593918054482717 / 2],
    [0.10128650732345633, 0.79742698535308720, 0.12593918054482717 / 2],
    [0.10128650732345633, 0.10128650732345633, 0.12593918054482717 / 2],
    [0.05971587178976981, 0.47014206410511505, 0.13239415278850618 / 2],
    [0.47014206410511505, 0.05971587178976981, 0.13239415278850618 / 2],
    [0.47014206410511505, 0.47014206410511505, 0.13239415278850618 / 2],
], dtype=np.float64)


def _tri_gl_duffy(n):
    """Triangle quadrature via Duffy product of two 1D Gauss-Legendre."""
    x_gl, w_gl = np.polynomial.legendre.leggauss(n)
    s = 0.5 * (x_gl + 1.0)
    t = 0.5 * (x_gl + 1.0)
    ws = 0.5 * w_gl
    wt = 0.5 * w_gl
    pts = np.empty((n * n, 3), dtype=np.float64)
    k = 0
    for i in range(n):
        for j in range(n):
            xi = s[i]
            eta = s[i] * t[j]
            w = ws[i] * wt[j] * s[i]
            pts[k, 0] = xi
            pts[k, 1] = eta
            pts[k, 2] = w
            k += 1
    return pts


def tri_quad(degree=5):
    if degree <= 5:
        return _TRI_DEG5.copy()
    n = max(2, (degree + 2) // 2)
    return _tri_gl_duffy(n)


# ----------------------------------------------------------------------
# RWG basis evaluation (orientation derived from global vertex indices)
# ----------------------------------------------------------------------
def _rwg_signs_from_global(tri_verts):
    """Compute the 3 RWG orientation signs for a triangle (CCW corner
    order) from its global vertex indices.

    Sigma_k = +1 when the CCW traversal of local edge k (from corner
    (k+1)%3 to corner (k+2)%3) goes from the SMALLER to the LARGER
    global vertex ID, -1 otherwise.
    """
    s = np.empty(3, dtype=np.int64)
    for k, (a, b) in enumerate(LOCAL_EDGE_ENDPOINTS):
        va, vb = int(tri_verts[a]), int(tri_verts[b])
        s[k] = 1 if va < vb else -1
    return s


def rwg_eval(verts, tri_verts, xi, eta, *, signs_override=None):
    """Evaluate the 3 RWG basis vectors on a flat triangle at given
    reference-coordinate points.

    Args:
        verts: (n_v, 3)
        tri_verts: (3,) global vertex indices to be used as local
            corners (in any order; need not be the natural CCW order).
        xi, eta: (n_q,) reference-triangle coords.
        signs_override: optional (3,) sign array.  If None, signs are
            computed from the local CCW direction of ``tri_verts``;
            this is correct ONLY if ``tri_verts`` is the natural CCW
            order.  When the caller has applied a non-CCW corner
            permutation (e.g. to bring shared corners to local 0, 1
            for the Sauter-Schwab common_edge canonical layout), the
            caller MUST pass the GLOBAL signs (computed from the
            ORIGINAL CCW ordering and permuted to the new local edge
            indexing) here -- otherwise the basis function sign flips
            on the shared edge (smoking gun: SL has negative
            eigenvalues).

    Returns:
        pts: (n_q, 3) physical positions.
        Js:  (n_q, 3, 3) RWG basis vectors -- Js[q, k, :] = J_local_k(r_q).
        A:   scalar physical-area of the triangle.
    """
    v0 = verts[tri_verts[0]]
    v1 = verts[tri_verts[1]]
    v2 = verts[tri_verts[2]]
    A = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))

    L0 = 1.0 - xi - eta
    L1 = xi
    L2 = eta
    pts = (L0[:, None] * v0
           + L1[:, None] * v1
           + L2[:, None] * v2)

    # NGSolve's HDivSurface(0) basis is RT_0-normalised (NO L_e factor):
    #     J_k(r) = sigma_k * (r - v_corner_k) / (2 A_T)
    # so that  div(J_k) = sigma_k / A_T  and  int_T div(J_k) dS = sigma_k.
    v_opp = np.stack([v0, v1, v2], axis=0)
    if signs_override is None:
        signs = _rwg_signs_from_global(tri_verts)
    else:
        signs = signs_override

    Js = np.empty((xi.shape[0], 3, 3), dtype=np.float64)
    for k in range(3):
        scale = signs[k] / (2.0 * A)
        Js[:, k, :] = scale * (pts - v_opp[k])
    return pts, Js, A


# ----------------------------------------------------------------------
# Pair classification
# ----------------------------------------------------------------------
def share_class(tri_a, tri_b):
    """Return ('identical' | 'common_edge' | 'common_vertex' | 'regular',
    sorted shared vertex list)."""
    sa = set(int(v) for v in tri_a)
    sb = set(int(v) for v in tri_b)
    shared = sorted(sa & sb)
    n = len(shared)
    if n == 3:
        return 'identical', shared
    if n == 2:
        return 'common_edge', shared
    if n == 1:
        return 'common_vertex', shared
    return 'regular', shared


# ----------------------------------------------------------------------
# Regular (non-singular) pair: tensor-product triangle Gauss
# ----------------------------------------------------------------------
def l_pair_regular(verts, tri_a, tri_b, q):
    """3x3 local Galerkin block for a regular triangle pair (T_a, T_b).

    Returns block[i, j] = (1/4 pi) ∫∫ (J_i^a · J_j^b) / |r - r'| dS dS'
    (positive Laplace kernel; caller scatters via global tri_edges).
    """
    xi = q[:, 0]
    eta = q[:, 1]
    w_ref = q[:, 2]
    pts_a, Js_a, A_a = rwg_eval(verts, tri_a, xi, eta)
    pts_b, Js_b, A_b = rwg_eval(verts, tri_b, xi, eta)

    # Physical-area weights: w_phys = 2 A * w_ref (per triangle).
    w_a = (2.0 * A_a) * w_ref
    w_b = (2.0 * A_b) * w_ref

    dx = pts_a[:, None, :] - pts_b[None, :, :]
    r = np.sqrt(np.sum(dx * dx, axis=-1))
    G = INV_4PI / r

    JJ = np.einsum('aix,bjx->abij', Js_a, Js_b)
    Wab = w_a[:, None] * G * w_b[None, :]
    block = np.einsum('ab,abij->ij', Wab, JJ)
    return block


# ----------------------------------------------------------------------
# Sauter-Schwab canonical permutations
# ----------------------------------------------------------------------
def _canonical_perm(case, tri_a, tri_b):
    """Return (perm_a, perm_b) bringing the triangles into Sauter-Schwab
    canonical orientation.

    SS canonical layout:

    - ``identical``:     T_a == T_b, same vertex order.
    - ``common_vertex``: shared corner at local 0 of both, CCW preserved.
    - ``common_edge``:   shared corners at local (0, 1) of BOTH, with
                         ``perm_a[0]_global == perm_b[0]_global`` and
                         ``perm_a[1]_global == perm_b[1]_global``.  In a
                         closed manifold with consistent outward
                         normals, T_b's natural CCW traverses the shared
                         edge OPPOSITE to T_a, so ``perm_b`` is in
                         general NOT a CCW rotation -- it locally
                         "flips" T_b.  The Sauter-Schwab singular
                         integrand for SL only depends on |r - r'|
                         (no normal term), so a non-CCW perm_b is fine.
                         RWG basis evaluation also tolerates this
                         because basis-function orientation signs are
                         derived from GLOBAL vertex IDs (not from the
                         local CCW order of the permuted triangle).
    """
    sa = [int(v) for v in tri_a]
    sb = [int(v) for v in tri_b]
    shared = set(sa) & set(sb)

    if case == 'identical':
        return [0, 1, 2], [0, 1, 2]

    if case == 'common_vertex':
        v_shared = next(iter(shared))
        ka = sa.index(v_shared)
        kb = sb.index(v_shared)
        perm_a = [(ka + i) % 3 for i in range(3)]
        perm_b = [(kb + i) % 3 for i in range(3)]
        return perm_a, perm_b

    if case == 'common_edge':
        # Find the position p_a in T_a where two CCW-consecutive corners
        # (sa[p_a], sa[(p_a+1)%3]) are both shared.
        p_a = None
        for p in range(3):
            if sa[p] in shared and sa[(p + 1) % 3] in shared:
                p_a = p
                break
        if p_a is None:
            raise RuntimeError(
                f"common_edge: shared corners not CCW-consecutive "
                f"in T_a = {sa}")
        v_first = sa[p_a]
        v_second = sa[(p_a + 1) % 3]
        # T_a CCW preserved: [first, second, third].
        perm_a = [p_a, (p_a + 1) % 3, (p_a + 2) % 3]
        # T_b: place SAME global vertex order at (0, 1) -- may flip CCW.
        q_first = sb.index(v_first)
        q_second = sb.index(v_second)
        q_third = 3 - q_first - q_second   # only remaining position
        perm_b = [q_first, q_second, q_third]
        return perm_a, perm_b

    raise ValueError(f"unknown case: {case!r}")


def _apply_corner_perm(tri_verts, perm):
    return np.array([int(tri_verts[p]) for p in perm], dtype=np.int64)


def _edge_perm_inverse(corner_perm):
    """Edge i is opposite corner i, so the edge permutation equals the
    corner permutation.  The inverse mapping (used to scatter the
    permuted block back to original local edge indices) is argsort.
    """
    return np.argsort(corner_perm)


# ----------------------------------------------------------------------
# Singular pair: SS Duffy via the validated scalar quadrature
# ----------------------------------------------------------------------
def l_pair_singular(verts, tri_a, tri_b, case, n_q=6):
    """3x3 local Galerkin block for a SINGULAR triangle pair via
    Sauter-Schwab Duffy.

    Re-uses the validated scalar quadrature from
    ``radia.bem.sibc_hacapk._ss_quad_pts`` -- only the basis evaluation
    differs (vector RWG instead of scalar hat).

    Sign-handling note (2026-05-02): when ``perm_b`` is non-CCW (the
    common_edge case in a closed manifold with consistent normals),
    the basis-function orientation signs MUST be computed from the
    original triangle CCW order and then permuted along with the
    edges -- recomputing from the permuted (non-CCW) order gives the
    OPPOSITE sign on the shared edge (because the local CCW
    traversal direction has flipped relative to the global vmin->vmax
    convention).  Smoking gun before this fix: SL had negative
    eigenvalues at order=0 on a closed surface.
    """
    perm_a, perm_b = _canonical_perm(case, tri_a, tri_b)
    tri_a_p = _apply_corner_perm(tri_a, perm_a)
    tri_b_p = _apply_corner_perm(tri_b, perm_b)

    # GLOBAL basis signs computed from the ORIGINAL CCW order,
    # then permuted to the new local edge indexing.  Edge index
    # equals corner index (edge k is opposite corner k), so the
    # edge permutation == corner permutation.
    signs_a_orig = _rwg_signs_from_global(tri_a)
    signs_b_orig = _rwg_signs_from_global(tri_b)
    signs_a_p = signs_a_orig[np.asarray(perm_a, dtype=np.int64)]
    signs_b_p = signs_b_orig[np.asarray(perm_b, dtype=np.int64)]

    xi_a, eta_a, xi_b, eta_b, w_q = _scalar_ss_quad_pts(case, n_q)

    pts_a, Js_a, A_a = rwg_eval(verts, tri_a_p, xi_a, eta_a,
                                 signs_override=signs_a_p)
    pts_b, Js_b, A_b = rwg_eval(verts, tri_b_p, xi_b, eta_b,
                                 signs_override=signs_b_p)

    diff = pts_a - pts_b
    r = np.sqrt(np.sum(diff * diff, axis=1))
    r = np.where(r == 0, 1e-300, r)

    # Per-quadpt scalar Galerkin weight (Sauter-Schwab w_q already
    # absorbs the singular Jacobian factor; the (4 A_a A_b) combines
    # the two reference-to-physical area maps).
    weight = INV_4PI * w_q * (4.0 * A_a * A_b) / r

    JJ = np.einsum('qix,qjx->qij', Js_a, Js_b)
    block_perm = np.einsum('q,qij->ij', weight, JJ)

    # Un-permute basis indices back to ORIGINAL local edge order.
    inv_a = _edge_perm_inverse(perm_a)
    inv_b = _edge_perm_inverse(perm_b)
    block = block_perm[np.ix_(inv_a, inv_b)]
    return block


# ----------------------------------------------------------------------
# Top-level dense assembly
# ----------------------------------------------------------------------
def assemble_M_rwg_dense(verts, tris, *, quad_degree=4):
    """Build the dense (n_e, n_e) mass matrix for the RWG basis.

    M_ij = int_S J_i(r) . J_j(r) dS

    Both basis functions live on the SAME surface (no kernel; this is
    a local assembly).  For RT_0-normalised RWG (J_k = sigma_k * (r -
    v_k) / (2 A_T) on T), the per-triangle local 3x3 block is
    integrated by Gauss quadrature; degree=4 is exact for the
    quadratic (r-v_k).(r-v_l) integrand.

    Used for SIBC dissipation (Z_s ∫ J . J dS gives Joule loss) and
    AC stabilised EFIE (Z_s/(j omega mu_0) Mass added to SL).
    """
    edges, tri_edges, _, _ = build_edge_topology(tris)
    n_e = len(edges)
    n_t = len(tris)
    M = np.zeros((n_e, n_e), dtype=np.float64)
    q = tri_quad(quad_degree)

    for t in range(n_t):
        Vt = tris[t]
        pts, Js, A = rwg_eval(verts, Vt, q[:, 0], q[:, 1])
        w_phys = (2.0 * A) * q[:, 2]
        # local M[k, l] = sum_q Js[q, k, :] . Js[q, l, :] * w_phys[q]
        local = np.einsum('q,qkx,qlx->kl', w_phys, Js, Js)
        Et = tri_edges[t]
        for k in range(3):
            for l in range(3):
                M[Et[k], Et[l]] += local[k, l]
    return M, edges, tri_edges


def build_div_matrix(tris, edges, tri_edges, verts):
    """Build the divergence matrix ``D[t, i] = int_T_t div(J_i) dS``.

    For RWG basis at order=0 the per-triangle divergence is constant::

        div(J_i) on T_t = sigma_(i, t) * L_i / A_T_t      (T_t adjacent to edge i)
                       = 0                                (else)

    integrated over T_t::

        D[t, i] = sigma_(i, t) * L_i

    where ``sigma_(i, t) = +1`` when the local CCW traversal of triangle
    t along its local edge corresponding to global edge i goes from the
    smaller-ID corner to the larger-ID corner, ``-1`` otherwise.

    Returns a (n_t, n_e) sparse-as-dense float64 matrix.
    """
    n_t = len(tris)
    n_e = len(edges)
    D = np.zeros((n_t, n_e), dtype=np.float64)
    for t in range(n_t):
        signs = _rwg_signs_from_global(tris[t])
        for k in range(3):
            ei = tri_edges[t, k]
            # RT_0 normalisation -> int_T div(J_i) dS = sigma_(i,t)
            D[t, ei] = signs[k]
    return D


def solve_inductance_source_sink_intree(verts, tris,
                                        source_tri_mask,
                                        sink_tri_mask, *,
                                        regular_quad_degree=11,
                                        singular_n_q=8,
                                        Z_s=0.0):
    """End-to-end self-inductance via the intree saddle-point EFIE.

    Replaces the ngsolve.bem path of
    ``examples/induction_heating/bem_reference/bem_inductance.py
    ::compute_inductance_source_sink`` with the intree assembler.  The
    saddle-point system is the Weggler stabilized EFIE::

        [ SL    D^T ] [ J ]   [ 0 ]
        [ D     0   ] [ p ] = [ g ]

    with ``g[t] = A_t / A_src`` for source tris and ``-A_t / A_snk``
    for sink tris (so the unit total current 1 A enters through source
    and exits through sink), and ``L = mu_0 * J^T SL J`` is the
    self-inductance for that unit current.

    Args:
        verts: (n_v, 3) float64.
        tris:  (n_t, 3) int.
        source_tri_mask: (n_t,) bool array, True for triangles on the
            source cap.
        sink_tri_mask:   (n_t,) bool array, True for triangles on the
            sink cap.
        regular_quad_degree, singular_n_q: forwarded to the SS Galerkin
            assembler.

    Returns:
        dict (mirrors compute_inductance_source_sink) with keys
        ``L``, ``J``, ``SL``, ``D``, ``A_source``, ``A_sink``,
        ``residual``, ``n_J``, ``n_t``.
    """
    from scipy.linalg import solve as scipy_solve

    MU_0 = 4.0e-7 * math.pi

    SL, edges, tri_edges = assemble_L_rwg_dense(
        verts, tris,
        regular_quad_degree=regular_quad_degree,
        include_singular=True,
        singular_n_q=singular_n_q)
    n_e = SL.shape[0]
    n_t = len(tris)

    D = build_div_matrix(tris, edges, tri_edges, verts)

    # Per-triangle areas
    A_t = np.empty(n_t)
    for t in range(n_t):
        v0 = verts[tris[t, 0]]
        v1 = verts[tris[t, 1]]
        v2 = verts[tris[t, 2]]
        A_t[t] = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))

    A_source = float(np.sum(A_t[source_tri_mask]))
    A_sink = float(np.sum(A_t[sink_tri_mask]))
    if A_source < 1e-30 or A_sink < 1e-30:
        return {"error": f"degenerate ports: A_src={A_source}, A_snk={A_sink}"}

    g = np.zeros(n_t, dtype=np.float64)
    g[source_tri_mask] = A_t[source_tri_mask] / A_source
    g[sink_tri_mask] = -A_t[sink_tri_mask] / A_sink

    # Saddle point with one constraint dropped (SurfaceL2 has a 1D
    # null space when the surface is closed: the constant function).
    D_red = D[:-1, :]
    g_red = g[:-1]
    n_constraint = D_red.shape[0]

    K = np.block([
        [SL,                       D_red.T],
        [D_red, np.zeros((n_constraint, n_constraint))]
    ])
    rhs = np.zeros(n_e + n_constraint, dtype=np.float64)
    rhs[n_e:] = g_red

    sol = scipy_solve(K, rhs)
    J = sol[:n_e]
    residual = float(np.max(np.abs(D @ J - g)))

    L_self = float(MU_0 * (J @ SL @ J))

    # DC SIBC: optional R extraction.  Z_s may be a scalar (uniform
    # sheet impedance, real for DC, complex for AC) or an (n_t,) array
    # (per-triangle).  R = J^T Re(Z_s) Mass J for unit prescribed
    # current (g normalised so that 1 A enters through source).
    R_dc = 0.0
    if not (np.isscalar(Z_s) and Z_s == 0.0):
        Mmat, _, _ = assemble_M_rwg_dense(verts, tris, quad_degree=4)
        Z_s_arr = np.asarray(Z_s)
        if Z_s_arr.ndim == 0:
            R_dc = float(np.real(Z_s_arr) * (J @ Mmat @ J))
        else:
            # Per-triangle Z_s -- requires per-triangle Mass blocks.
            # For DC scalar uniform, the simple form above suffices.
            raise NotImplementedError(
                "Per-triangle Z_s is not yet wired through the "
                "assembled Mass matrix; pass a scalar Z_s for now.")
    else:
        Mmat = None

    return {
        "L": L_self,
        "R": R_dc,
        "J": J,
        "SL": SL,
        "D": D,
        "Mass": Mmat,
        "A_source": A_source,
        "A_sink": A_sink,
        "residual": residual,
        "n_J": n_e,
        "n_t": n_t,
        "edges": edges,
        "tri_edges": tri_edges,
    }


def compute_dc_resistance_intree(verts, tris, source_tri_mask, sink_tri_mask,
                                  *, Z_s, mass_quad_degree=4):
    """End-to-end DC resistance via the (Z_s * Mass)-based saddle point.

    For a thin conductor with surface impedance Z_s = R_sheet
    (Ohm/square, real for DC; complex for AC ohmic loss), solves the
    DC current-flow problem on the boundary mesh::

        [ Z_s Mass    D^T ] [ J ]   [ 0 ]
        [ D            0  ] [ p ] = [ g ]

    The unknown ``J`` is the surface current that *minimises ohmic
    loss* (Z_s J^T Mass J / 2) subject to current conservation
    ``D J = g`` (unit total current entering through ``source`` and
    exiting through ``sink``).  This is fundamentally different from
    the PEC inductance solution: at DC the conductor sees no inductive
    impedance, only resistance, and the J distribution looks like the
    DC current flow on a sheet.

    For PCB trace resistance (R_sheet known per material/thickness)
    and arbitrary 3D conductor DC R, this method gives the DC R
    directly.  The companion ``solve_inductance_source_sink_intree``
    gives the high-frequency PEC inductance (geometric L).

    Args:
        verts, tris: surface mesh.
        source_tri_mask, sink_tri_mask: (n_t,) bool masks.
        Z_s: scalar surface impedance [Ohm/square] (real for DC).
        mass_quad_degree: Gauss quadrature degree for Mass assembly.

    Returns:
        dict with keys::
            R          : DC resistance [Ohm]
            J          : surface current vector
            Mass       : (n_e, n_e) mass matrix
            D          : (n_t, n_e) divergence matrix
            A_source, A_sink : port areas
            residual   : ||D J - g||_inf
            n_J, n_t   : DOF / triangle counts
    """
    from scipy.linalg import solve as scipy_solve

    edges, tri_edges, _, _ = build_edge_topology(tris)
    n_e = len(edges)
    n_t = len(tris)

    Mmat, _, _ = assemble_M_rwg_dense(verts, tris,
                                       quad_degree=mass_quad_degree)
    D = build_div_matrix(tris, edges, tri_edges, verts)

    # Per-triangle areas for port g.
    A_t = np.empty(n_t)
    for t in range(n_t):
        v0 = verts[tris[t, 0]]
        v1 = verts[tris[t, 1]]
        v2 = verts[tris[t, 2]]
        A_t[t] = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
    A_source = float(np.sum(A_t[source_tri_mask]))
    A_sink = float(np.sum(A_t[sink_tri_mask]))
    if A_source < 1e-30 or A_sink < 1e-30:
        return {"error": f"degenerate ports: A_src={A_source}, A_snk={A_sink}"}
    g = np.zeros(n_t)
    g[source_tri_mask] = A_t[source_tri_mask] / A_source
    g[sink_tri_mask] = -A_t[sink_tri_mask] / A_sink

    # Drop one constraint row (SurfaceL2 has 1D null space on closed surface).
    D_red = D[:-1, :]
    g_red = g[:-1]
    n_constraint = D_red.shape[0]

    Z_real = float(np.real(Z_s))
    K = np.block([
        [Z_real * Mmat,           D_red.T],
        [D_red, np.zeros((n_constraint, n_constraint))]
    ])
    rhs = np.zeros(n_e + n_constraint)
    rhs[n_e:] = g_red

    sol = scipy_solve(K, rhs)
    J = sol[:n_e]
    residual = float(np.max(np.abs(D @ J - g)))

    R_dc = float(Z_real * (J @ Mmat @ J))

    return {
        "R": R_dc,
        "J": J,
        "Mass": Mmat,
        "D": D,
        "A_source": A_source,
        "A_sink": A_sink,
        "residual": residual,
        "n_J": n_e,
        "n_t": n_t,
        "edges": edges,
        "tri_edges": tri_edges,
    }


def solve_ac_efie_sibc_intree(verts, tris,
                                source_tri_mask, sink_tri_mask, *,
                                omega, Z_s,
                                regular_quad_degree=11,
                                singular_n_q=8,
                                mass_quad_degree=4):
    """AC EFIE+SIBC saddle-point solve at frequency ``omega`` with
    surface impedance ``Z_s`` (complex scalar).

    Solves the complex saddle::

        [ j omega mu_0 SL + Z_s Mass    D^T ] [ J ]   [ 0 ]
        [ D                              0   ] [ p ] = [ g ]

    The SL term gives the magnetic stored energy, the Z_s Mass term
    gives the surface-impedance loss (skin effect via Leontovich SIBC).
    For Z_s = (1 + j) / (sigma * delta) at frequency omega in a
    non-magnetic conductor, the system captures the thin-skin
    distribution self-consistently -- distinct from the PEC-limit
    `solve_inductance_source_sink_intree` which assumes Z_s = 0 and
    extracts only the geometric L.

    Args:
        verts, tris: surface mesh.
        source_tri_mask, sink_tri_mask: (n_t,) bool masks for the
            two ports.
        omega: angular frequency [rad/s], > 0.
        Z_s: complex surface impedance [Ohm/sq].
        regular_quad_degree, singular_n_q: SS Galerkin SL parameters.
        mass_quad_degree: triangle Gauss degree for Mass assembly.

    Returns:
        dict with::
            L_port  : Im(Z_port) / omega  -- self inductance [H]
            R_port  : Re(Z_port)          -- AC resistance [Ohm]
            Z_port  : complex impedance R + j*omega*L
            J       : complex (n_e,) surface current vector
            SL, Mass, D : assembled matrices
            residual : ||D Re(J) - g||_inf
            n_J, n_t  : DOF counts
    """
    from scipy.linalg import solve as scipy_solve

    if omega <= 0:
        raise ValueError("solve_ac_efie_sibc_intree requires omega > 0; "
                         "use compute_dc_resistance_intree for DC.")

    SL, edges, tri_edges = assemble_L_rwg_dense(
        verts, tris,
        regular_quad_degree=regular_quad_degree,
        include_singular=True,
        singular_n_q=singular_n_q)
    n_e = SL.shape[0]
    n_t = len(tris)

    Mmat, _, _ = assemble_M_rwg_dense(verts, tris,
                                       quad_degree=mass_quad_degree)
    D = build_div_matrix(tris, edges, tri_edges, verts)

    # Per-triangle areas + port g
    A_t = np.empty(n_t)
    for t in range(n_t):
        v0 = verts[tris[t, 0]]
        v1 = verts[tris[t, 1]]
        v2 = verts[tris[t, 2]]
        A_t[t] = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
    A_source = float(np.sum(A_t[source_tri_mask]))
    A_sink = float(np.sum(A_t[sink_tri_mask]))
    if A_source < 1e-30 or A_sink < 1e-30:
        return {"error": f"degenerate ports: A_src={A_source}, A_snk={A_sink}"}

    g = np.zeros(n_t)
    g[source_tri_mask] = A_t[source_tri_mask] / A_source
    g[sink_tri_mask] = -A_t[sink_tri_mask] / A_sink

    D_red = D[:-1, :]
    g_red = g[:-1]
    n_constraint = D_red.shape[0]

    MU_0 = 4.0e-7 * math.pi
    A_complex = (1j * omega * MU_0) * SL + complex(Z_s) * Mmat
    K = np.zeros((n_e + n_constraint, n_e + n_constraint), dtype=complex)
    K[:n_e, :n_e] = A_complex
    K[:n_e, n_e:] = D_red.T
    K[n_e:, :n_e] = D_red

    rhs = np.zeros(n_e + n_constraint, dtype=complex)
    rhs[n_e:] = g_red

    sol = scipy_solve(K, rhs)
    J = sol[:n_e]                          # complex
    p = sol[n_e:]
    residual = float(np.max(np.abs(D @ J.real - g)))

    # Energy / dissipation extraction (real symmetric SL, Mass:
    # J^H SL J and J^H Mass J are real positive).
    J_H_SL_J = float(np.real(np.vdot(J, SL @ J)))
    J_H_Mass_J = float(np.real(np.vdot(J, Mmat @ J)))

    L_port = MU_0 * J_H_SL_J
    R_port = float(np.real(Z_s)) * J_H_Mass_J
    Z_port = R_port + 1j * omega * L_port

    return {
        "L_port": float(L_port),
        "R_port": float(R_port),
        "Z_port": complex(Z_port),
        "J": J,
        "SL": SL,
        "Mass": Mmat,
        "D": D,
        "A_source": A_source,
        "A_sink": A_sink,
        "residual": residual,
        "J_H_SL_J": J_H_SL_J,
        "J_H_Mass_J": J_H_Mass_J,
        "n_J": n_e,
        "n_t": n_t,
        "edges": edges,
        "tri_edges": tri_edges,
    }


def assemble_L_rwg_dense(verts, tris, *,
                         regular_quad_degree=11,
                         include_singular=True,
                         singular_n_q=6):
    """Build the dense (n_e, n_e) Galerkin Laplace L matrix on RWG basis.

    Args:
        verts: (n_v, 3) float64 vertex coords.
        tris:  (n_t, 3) int triangle vertex indices (CCW).
        regular_quad_degree: triangle Gauss degree for non-adjacent
            pairs (default 11 -> ~36-pt Duffy product per tri).
        include_singular: when False, singular pairs are skipped.
        singular_n_q: 1D Gauss-Legendre order in each Duffy direction
            for the singular cases (n_q=6 -> 1296 pts per sub-cube,
            6/5/2 sub-cubes for identical/common_edge/common_vertex).

    Returns:
        L: (n_e, n_e) float64 dense matrix.
        edges: (n_e, 2) global edge endpoint indices.
        tri_edges: (n_t, 3) global edge index per local edge.
    """
    edges, tri_edges, _, _ = build_edge_topology(tris)
    n_e = len(edges)
    n_t = len(tris)
    L = np.zeros((n_e, n_e), dtype=np.float64)
    q = tri_quad(regular_quad_degree)

    for a in range(n_t):
        Va = tris[a]
        for b in range(n_t):
            Vb = tris[b]
            cls, _shared = share_class(Va, Vb)
            if cls == 'regular':
                block = l_pair_regular(verts, Va, Vb, q)
            else:
                if not include_singular:
                    continue
                block = l_pair_singular(
                    verts, Va, Vb, cls, n_q=singular_n_q)
            for i in range(3):
                ei = tri_edges[a, i]
                for j in range(3):
                    ej = tri_edges[b, j]
                    L[ei, ej] += block[i, j]
    return L, edges, tri_edges
