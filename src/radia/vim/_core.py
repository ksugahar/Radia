"""radia.vim._core -- HDiv-type VIM demag SPARSE assembly on UNSTRUCTURED meshes (tet/hex/wedge).

The NGSolve HDiv(order=0) extraction is element-AGNOSTIC:
  B = [ -div(u) tested on L2 (volume charge rho) ;  u.Trace().n tested on SurfaceL2 (surface charge sigma) ]
  M_mass = the HDiv(0) mass.
`build_demag` returns the SPARSE charge map B + the SPARSE HDiv mass M_mass + the per-cell / per-face
geometry that the C++ analytic charge-Gram H-matrix (`_ChargeGramHMatrix`) consumes.  The C++ kernel IS
the demag operator (N v = B^T (H.matvec(B v))) and folds IMA (image charges) in via image_masks /
image_signs -- the dense Python charge-Gram path (the O(N^2) Gram G, the SVD loop basis, the dense
N = B^T G B, and the analytic / image / Wilton dense Gram builders) was REMOVED.

Only the Gram self-energy geometry (the preconditioner diagonal) is element-specific -> TET (barycentric
sub-points) and TRIANGLE faces, with the self constants c_tet / c_tri, plus the analytic polytope
self-energy for hex/wedge cells + quad faces.
"""
from math import pi, sqrt

import numpy as np
import scipy.sparse as sp
from scipy.spatial import ConvexHull

import ngsolve as ng

EPS = 1e-6

# ---- self constants c = INT INT_{unit-measure shape} 1/|x-y| (sub-cell self correction) ----
# cube/square via exact quadrature (as the hex prototype); tet/tri extrapolated (compute_self_constants.py)
from scipy import integrate  # noqa: E402
_c_sq, _ = integrate.dblquad(lambda v, u: (1 - u) * (1 - v) / sqrt(u * u + v * v + 1e-300), 0, 1, 0, 1)
C_SQ = 4.0 * _c_sq
C_TET = 1.77   # INT INT 1/|x-y| over a unit-volume regular tet (Richardson, ~0.3% uncertainty)
C_TRI = 2.89   # INT INT 1/|x-y| over a unit-area equilateral triangle


def _bary_tet(nsub):
    """Barycentric sub-point lattice of a reference tet: (i,j,k,l) lattice nodes with i+j+k+l=nsub-1,
    point (i+1/4, j+1/4, k+1/4, l+1/4)/nsub.  The 4 coords sum to EXACTLY 1 (inside the tet) and
    nsub=1 -> the centroid (1/4,1/4,1/4,1/4).  Equal sub-weights.  (The earlier [.+0.5,...]/nsub form
    summed to 1+1/nsub -> sub-points OUTSIDE the tet -- a quadrature bug; fixed 2026-06-07.)"""
    lam = []
    for i in range(nsub):
        for j in range(nsub - i):
            for k in range(nsub - i - j):
                l = nsub - 1 - i - j - k
                lam.append([(i + 0.25) / nsub, (j + 0.25) / nsub, (k + 0.25) / nsub, (l + 0.25) / nsub])
    return np.array(lam)   # (m, 4) barycentric, each row sums to 1


def _bary_tri(nsub):
    """Barycentric sub-point lattice of a reference triangle: (i,j,k) with i+j+k=nsub-1, point
    (i+1/3, j+1/3, k+1/3)/nsub (sums to EXACTLY 1; nsub=1 -> centroid).  Equal sub-weights."""
    lam = []
    for i in range(nsub):
        for j in range(nsub - i):
            k = nsub - 1 - i - j
            lam.append([(i + 1.0 / 3) / nsub, (j + 1.0 / 3) / nsub, (k + 1.0 / 3) / nsub])
    return np.array(lam)   # (m, 3), each row sums to 1


def _gauss_duffy_tet(o):
    """Gauss-Duffy collapsed-cube tet rule: o Gauss-Legendre pts/dim (o^3 nodes), ref-tet coords
    (lam1,lam2,lam3) + weights summing to 1/6.  For the SMOOTH outer integral of an EXACT analytic inner
    (tri_potential) this converges to ~machine precision.  Used by the polytope cell self-energy
    quadrature (_cell_outer_quad).  Same rule as radia.vim._vim._tet_ref(o)."""
    x, w = np.polynomial.legendre.leggauss(o)
    s, ws = 0.5 * (x + 1.0), 0.5 * w
    P, W = [], []
    for a, wa in zip(s, ws):
        for b, wb in zip(s, ws):
            for c, wc in zip(s, ws):
                P.append((a, b * (1 - a), c * (1 - a) * (1 - b)))
                W.append(wa * wb * wc * (1 - a) ** 2 * (1 - b))
    return np.array(P), np.array(W)   # (m,3) ref pts (lam1,lam2,lam3), weights sum 1/6


def tet_self_energy(V, vol, nsub):
    """G_aa for a tet (4 verts V) of volume `vol`: cross sub-point sum + sub-cell self (c_tet).

    NOTE (2026-06-12): this empirical cross-subpoint + C_TET*w^(5/3) self is ~10-15% LOW -- the same
    volume-self quadrature defect that was fixed in the C++ analytic ctor (the production order-0 Gram).
    It is NOT fixed here because it feeds only the PRECONDITIONER diagonal, never the production demag
    operator (which uses the C++ ctor's own QuadDot(a,a))."""
    lam = _bary_tet(nsub)
    C = lam @ V                                   # (m,3) sub-point positions
    w = np.full(len(C), vol / len(C))             # equal sub-weights, sum = vol
    D = np.linalg.norm(C[:, None] - C[None, :], axis=2)
    np.fill_diagonal(D, np.inf)
    cross = np.sum(np.outer(w, w) / (4 * pi * D))
    selfsub = np.sum(C_TET * w ** (5.0 / 3.0) / (4 * pi))
    return cross + selfsub


def tri_self_energy(V, area, nsub):
    lam = _bary_tri(nsub)
    C = lam @ V
    w = np.full(len(C), area / len(C))
    D = np.linalg.norm(C[:, None] - C[None, :], axis=2)
    np.fill_diagonal(D, np.inf)
    cross = np.sum(np.outer(w, w) / (4 * pi * D))
    selfsub = np.sum(C_TRI * w ** 1.5 / (4 * pi))
    return cross + selfsub


# --- exact analytic flat-triangle potential (tri_potential) + Dunavant outer quadrature ---
# tri_potential is the exact Newtonian potential of a uniformly-charged flat triangle (Wilton et al.,
# IEEE TAP 32(3):276 (1984); Graglia, IEEE TAP 41(10):1448 (1993)).  It feeds the polytope / face
# self-energy diagonal (the preconditioner) via _polytope_potential / _poly_face_self_energy.

# Dunavant degree-5 symmetric triangle rule (7 points, barycentric; weights sum to 1) for the OUTER
# integral over the observation triangle (the inner integral is the exact triangle potential).
_DUN5 = [(1.0 / 3, 1.0 / 3, 1.0 / 3, 0.225)]
for _a, _w in [(0.0597158717, 0.1323941527), (0.7974269853, 0.1259391805)]:
    _b = (1.0 - _a) / 2.0
    _DUN5 += [(_a, _b, _b, _w), (_b, _a, _b, _w), (_b, _b, _a, _w)]
_DUN5 = np.array(_DUN5)


def tri_potential(V, r):
    """Exact INT_T 1/|r - r'| dA' over a flat triangle T (vertices V, 3x3) at observation point(s) r
    (constant unit density) -- the Wilton/Graglia analytic potential integral.  Verified vs fine
    numerical quadrature to the reference's own accuracy (~1e-4, limited by the reference near the
    plane).  Handles r on either side of / on the triangle plane.

    Vectorized: r may be a single point (3,) -> float, or a batch (M,3) -> (M,) -- the batch form is
    what makes the O(n_bf^2) surface Gram practical."""
    r = np.asarray(r, float)
    single = (r.ndim == 1)
    R = r.reshape(1, 3) if single else r            # (M,3)
    v0, v1, v2 = V[0], V[1], V[2]
    n = np.cross(v1 - v0, v2 - v0)
    n = n / np.linalg.norm(n)
    d = (R - v0) @ n                                # (M,) signed heights above the plane
    p = R - d[:, None] * n                          # (M,3) projections onto the plane
    ad = np.abs(d)
    I = np.zeros(len(R))
    vs = (v0, v1, v2)
    for i in range(3):
        a = vs[i]
        b = vs[(i + 1) % 3]
        lh = (b - a) / np.linalg.norm(b - a)
        uh = np.cross(lh, n)                         # in-plane unit normal to the edge
        P0 = (a - p) @ uh                            # (M,) signed perp distance projection -> edge
        sm = (a - p) @ lh                            # (M,) endpoint projections along the edge
        sp = (b - p) @ lh
        Rm = np.linalg.norm(R - a, axis=1)           # (M,) = |r - a|
        Rp = np.linalg.norm(R - b, axis=1)
        R0sq = P0 * P0 + d * d
        dm, dp = Rm + sm, Rp + sp
        safe = (dp > 1e-300) & (dm > 1e-300)
        ratio = np.where(safe, dp, 1.0) / np.where(safe, dm, 1.0)
        f = np.where(safe, np.log(ratio), 0.0)
        beta = np.arctan2(P0 * sp, R0sq + ad * Rp) - np.arctan2(P0 * sm, R0sq + ad * Rm)
        I += P0 * f - ad * beta
    return float(I[0]) if single else I


# --- polytope (hex / wedge cell, quad face) charge-Gram geometry ---------------------------------
# The exact triangle potential (tri_potential) generalizes to ANY flat-faced convex cell + quad face
# with NO new singular quadrature: the cell Newtonian potential is the divergence-theorem sum over the
# cell's (convex-hull) triangular faces of the SAME exact triangle potential, and a quad face is two
# flat triangles.  These helpers build the polytope self-energy diagonal (the preconditioner) and the
# flat triangle-soup geometry (polytope_flat_geom) that the C++ POLYTOPE _ChargeGramHMatrix consumes.


def _cell_hull_tris(V):
    """(triangle (3,3), outward unit normal) faces of a convex cell (hex/wedge/...) via ConvexHull.
    Feeds the polytope volume-charge Newtonian potential (the self-energy diagonal + the C++ polytope
    charge-Gram triangle soup)."""
    V = np.asarray(V, float)
    hull = ConvexHull(V)
    cen = V.mean(0)
    tris = []
    for simplex in hull.simplices:
        P = V[simplex]
        nrm = np.cross(P[1] - P[0], P[2] - P[0])
        ln = np.linalg.norm(nrm)
        if ln < 1e-300:
            continue
        nrm = nrm / ln
        if np.dot(nrm, P.mean(0) - cen) < 0:
            nrm = -nrm                                       # outward
        tris.append((P, nrm))
    return tris


def _polytope_potential(tris, R):
    """VECTORIZED INT_cell 1/|R-r'| dV' over a flat-faced convex cell (hull faces `tris`) at points R,
    via the divergence theorem on the (vectorized) exact triangle potential (tri_potential).
    R: (3,)->float, (M,3)->(M,)."""
    R = np.asarray(R, float)
    single = (R.ndim == 1)
    RR = R.reshape(1, 3) if single else R
    tot = np.zeros(len(RR))
    for (P, nrm) in tris:
        d = (P[0] - RR) @ nrm                                # (M,) signed distance R -> face plane
        tot += d * tri_potential(P, RR)
    tot *= 0.5
    return float(tot[0]) if single else tot


def _face_subtris(V):
    """Sub-triangles of a boundary face: triangle (3 verts) -> [V]; quad (4 verts) -> 2-triangle fan."""
    V = np.asarray(V, float)
    if len(V) == 3:
        return [V]
    if len(V) == 4:
        return [V[[0, 1, 2]], V[[0, 2, 3]]]
    raise ValueError("HDiv-VIM: boundary face with %d vertices is not supported (triangle or quad only)"
                     % len(V))


def _cell_outer_quad(V, vol, lam_t, w_t):
    """Interior quadrature (points (m,3), weights (m,) summing to `vol`) over a convex cell: centroid-fan
    tets (centroid + each hull face), each filled by the Gauss-Duffy tet rule (lam_t, w_t)."""
    V = np.asarray(V, float)
    cen = V.mean(0)
    Ps, Ws = [], []
    for (T, _n) in _cell_hull_tris(V):
        tet = np.vstack([cen, T])
        tvol = abs(np.linalg.det(tet[1:] - tet[0])) / 6.0
        Ps.append(tet[0] + lam_t @ (tet[1:] - tet[0]))
        Ws.append(w_t * 6.0 * tvol)
    return np.vstack(Ps), np.concatenate(Ws)


def _face_outer_quad(V, area):
    """Outer quadrature (points, weights summing to `area`) over a boundary face: Dunavant-5 per sub-tri."""
    Ps, Ws = [], []
    for T in _face_subtris(V):
        a = 0.5 * np.linalg.norm(np.cross(T[1] - T[0], T[2] - T[0]))
        Ps.append(_DUN5[:, :3] @ T)
        Ws.append(_DUN5[:, 3] * a)
    return np.vstack(Ps), np.concatenate(Ws)


def _poly_cell_self_energy(V, vol):
    """Analytic volume self-energy G_aa = (1/4pi) INT_cell INT_cell 1/r for a convex polytope cell (its
    own outer quadrature against its own Newtonian potential) -- the hex/wedge analog of tet_self_energy,
    but EXACT-quadrature rather than the C_TET extrapolation.  Used only for the monopole-path
    preconditioner diagonal (the analytic Gram computes its own diagonal)."""
    lam_t, w_t = _gauss_duffy_tet(4)
    P, W = _cell_outer_quad(V, vol, lam_t, w_t)
    return float(np.sum(W * _polytope_potential(_cell_hull_tris(V), P)) / (4.0 * pi))


def _poly_face_self_energy(V, area):
    """Analytic surface self-energy G_aa for a quad/triangle boundary face (Dunavant outer against its
    sub-triangle Wilton potentials).  Monopole-path preconditioner diagonal only."""
    P, W = _face_outer_quad(V, area)
    pot = sum(tri_potential(T, P) for T in _face_subtris(V))
    return float(np.sum(W * pot) / (4.0 * pi))


def parse_image_string(image):
    """Parse a Radia IMA image string ('+x-z', '-x+y-z', ...) into a list of (axis, sign): axis in
    {0,1,2} for x/y/z, sign +1 (symmetric across that plane) / -1 (antisymmetric).  The string is a
    concatenation of (+|-)(x|y|z) tokens.  Matches the rad.Solve(image=...) convention (CLAUDE.md
    IMA Sign Selection: field PARALLEL to the mirror -> '+', PERPENDICULAR -> '-')."""
    s = image.strip().lower().replace(" ", "")
    planes, i = [], 0
    axis_of = {"x": 0, "y": 1, "z": 2}
    while i < len(s):
        if s[i] not in "+-" or i + 1 >= len(s) or s[i + 1] not in axis_of:
            raise ValueError("bad IMA image string %r (expected tokens like '+x','-z')" % image)
        planes.append((axis_of[s[i + 1]], 1 if s[i] == "+" else -1))
        i += 2
    axes = [a for a, _ in planes]
    if len(set(axes)) != len(axes):
        raise ValueError("IMA image string %r repeats an axis" % image)
    return planes


def image_group(planes):
    """Image-method reflection group from parsed planes [(axis,sign),...]: every NON-EMPTY subset of
    the mirror planes gives one image, whose reflection negates the coords on those axes and whose
    scalar sign is the PRODUCT of the per-plane signs (standard image method; the paper eq 23
    G^IMA = G + sum_m s_m G_{i,m(j)}).  Returns [(axes_tuple, sign), ...] (2^P - 1 images)."""
    out = []
    P = len(planes)
    for mask in range(1, 1 << P):
        axes, sign = [], 1
        for k in range(P):
            if mask & (1 << k):
                axes.append(planes[k][0]); sign *= planes[k][1]
        out.append((tuple(sorted(axes)), sign))
    return out


def polytope_flat_geom(el_V, bf_V, el_vol, bf_area):
    """Flatten the per-cell convex-hull triangles + per-face sub-triangles (+ vertex-mean centroids and
    measures) into the flat triangle-soup arrays the C++ POLYTOPE _ChargeGramHMatrix constructor consumes
    (hex/wedge scalable path).  Cells -> _cell_hull_tris, boundary faces -> _face_subtris; the C++ rebuilds
    the centroid-fan (cells) / Dunavant (faces) outer quadrature and the divergence-theorem source potential
    from these.  Returns flat float arrays + int CSR offsets (triangle soup = 9 doubles/tri)."""
    cell_tris, cell_troff, cell_cent = [], [0], []
    for V in el_V:
        V = np.asarray(V, float)
        cell_cent.append(V.mean(0))
        tris = _cell_hull_tris(V)
        for (P, _n) in tris:
            cell_tris.append(np.asarray(P, float).ravel())
        cell_troff.append(cell_troff[-1] + len(tris))
    face_tris, face_troff, face_cent = [], [0], []
    for V in bf_V:
        V = np.asarray(V, float)
        face_cent.append(V.mean(0))
        subs = _face_subtris(V)
        for T in subs:
            face_tris.append(np.asarray(T, float).ravel())
        face_troff.append(face_troff[-1] + len(subs))
    return dict(
        cell_tris=(np.concatenate(cell_tris) if cell_tris else np.zeros(0)),
        cell_troff=np.asarray(cell_troff, np.int32),
        cell_cent=(np.asarray(cell_cent, float).ravel() if cell_cent else np.zeros(0)),
        cell_meas=np.asarray(el_vol, float),
        face_tris=(np.concatenate(face_tris) if face_tris else np.zeros(0)),
        face_troff=np.asarray(face_troff, np.int32),
        face_cent=(np.asarray(face_cent, float).ravel() if face_cent else np.zeros(0)),
        face_meas=np.asarray(bf_area, float),
    )


def _csr_sp(bf):
    """NGSolve BilinearForm.mat -> scipy CSR (SPARSE). The build keeps the charge maps + HDiv mass
    sparse (the C++ charge-Gram H-matrix is the demag operator; no dense n_charge^2 object)."""
    m = bf.mat
    r, c, v = m.COO()
    return sp.csr_matrix((np.array(v), (np.array(r), np.array(c))), shape=(m.height, m.width))


def build_demag(mesh, nsub=4):
    """Assemble the SPARSE pieces of the HDiv-type VIM demag operator + the C++ charge-Gram geometry.

    Returns the sparse charge map B (= the scaled charge map, each row / its measure), the sparse HDiv
    mass M_mass, the diagonal charge self-energies, the charge centroids/measures, and the per-cell /
    per-face vertex geometry that the C++ analytic charge-Gram H-matrix (`_ChargeGramHMatrix`) consumes.
    There is NO dense N: the dense Python charge-Gram path (the O(N^2) Gram G, the SVD loop basis, the
    dense N = B^T G B, and the analytic / image / Wilton dense Gram builders) was REMOVED -- the C++
    `_ChargeGramHMatrix` kernel is the demag operator, and IMA (image charges) folds INTO it via the
    image_masks / image_signs arguments (see _solve.py).  The demag apply is N v = B^T (H.matvec(B v))
    with H = the C++ charge-Gram H-matrix.

    nsub controls the sub-point lattice for the diagonal charge self-energies (the preconditioner
    diagonal only); the C++ charge Gram computes its own diagonal.

    Per the caller-wraps policy this HELPER does NOT open a TaskManager; the CALLER wraps build_demag in
    `with ng.TaskManager():` (CLAUDE.md "TaskManager Wrap Policy: Caller Wraps, Helper Does NOT")."""
    # NGSolve assembly (FES, charge-map BilinearForms, HDiv mass).
    fes = ng.HDiv(mesh, order=0)
    ndof = fes.ndof
    nn = ng.specialcf.normal(mesh.dim)
    L2v, L2b = ng.L2(mesh, order=0), ng.SurfaceL2(mesh, order=0)
    u = fes.TrialFunction()
    bv = ng.BilinearForm(trialspace=fes, testspace=L2v)
    bv += (-ng.div(u)) * L2v.TestFunction() * ng.dx
    bv.Assemble()
    bb = ng.BilinearForm(trialspace=fes, testspace=L2b)
    bb += (u.Trace() * nn) * L2b.TestFunction() * ng.ds
    bb.Assemble()
    Bv_sp, Bb_sp = _csr_sp(bv), _csr_sp(bb)   # charge maps as sparse CSR
    massv = ng.BilinearForm(L2v); massv += L2v.TrialFunction() * L2v.TestFunction() * ng.dx; massv.Assemble()
    massb = ng.BilinearForm(L2b); massb += L2b.TrialFunction() * L2b.TestFunction() * ng.ds; massb.Assemble()
    el_vol = _csr_sp(massv).diagonal(); bf_area = _csr_sp(massb).diagonal()   # L2-0 mass is diagonal == measures
    # HDiv mass (the physical demag-factor metric) -- kept SPARSE (no dense n_charge^2 object anywhere)
    vh = fes.TestFunction()
    mh = ng.BilinearForm(fes); mh += u * vh * ng.dx; mh.Assemble()
    M_mass_sp = _csr_sp(mh)
    el_V = [np.array([mesh[v].point for v in el.vertices]) for el in mesh.Elements(ng.VOL)]
    bf_V = [np.array([mesh[v].point for v in el.vertices]) for el in mesh.Elements(ng.BND)]
    # uniform M_z, L2-projected onto RT0 (for the demag Rayleigh quotient)
    gfu = ng.GridFunction(fes); gfu.Set(ng.CoefficientFunction((0, 0, 1))); m_unit = np.array(gfu.vec)
    el_c = np.array([V.mean(0) for V in el_V]); bf_c = np.array([V.mean(0) for V in bf_V])
    n_el, n_bf = len(el_c), len(bf_c)
    cent = np.vstack([el_c, bf_c]); meas = np.concatenate([el_vol, bf_area])
    # scaled charge map B (each row / its measure) as sparse CSR (each charge cell touches only its own
    # faces) -- this is both B and B_csr in the returned dict.
    _Brows = [sp.diags(1.0 / el_vol) @ Bv_sp]
    if n_bf:
        _Brows.append(sp.diags(1.0 / bf_area) @ Bb_sp)
    B_sp = sp.vstack(_Brows).tocsr()

    # diagonal charge self-energies (O(N); the preconditioner diagonal -- the C++ charge Gram computes
    # its own diagonal).  tet/triangle use the validated tet_self_energy/tri_self_energy; hex/wedge cells
    # + quad faces use the analytic polytope self-energy (EXACT-quadrature) so a non-tet mesh does not
    # crash here and gets a physical diagonal.
    diagG = np.empty(n_el + n_bf)
    for k, V in enumerate(el_V):
        diagG[k] = tet_self_energy(V, el_vol[k], nsub) if len(V) == 4 else _poly_cell_self_energy(V, el_vol[k])
    for k, V in enumerate(bf_V):
        diagG[n_el + k] = (tri_self_energy(V, bf_area[k], nsub) if len(V) == 3
                           else _poly_face_self_energy(V, bf_area[k]))
    # SPARSE-ONLY: the C++ analytic _ChargeGramHMatrix is the demag operator (N v = B^T (H.matvec(B v))),
    # so there is NO dense n_charge^2 object here -- M_mass + B stay SPARSE and the charge-Gram BUILD needs
    # only the verts (straight from the mesh).  The dense Python Gram path (G, the loop SVD, N = B^T G B,
    # the analytic / image / Wilton dense Gram builders) was REMOVED.
    M_mass = M_mass_sp            # sparse CSR (RT0 HDiv mass is local -> sparse)
    B = B_sp                      # sparse CSR scaled charge map
    # charge geometry (for the C++ HACApK charge-Gram H-matrix path, #1b): centroids, measures,
    # diagonal self-energies, and the sparse charge map B as scipy CSR.
    # cell_verts [n_el*12] (tets, 4 verts) + face_verts [n_bf*9] (tris, 3 verts) feed the ANALYTIC
    # C++ charge Gram (the TET cells + triangle faces path); ravel only when every cell/face has the same
    # vertex count (uniform tet/triangle mesh).  A hex/wedge or mixed mesh leaves these empty -- it takes
    # the C++ polytope triangle-soup path (poly, below) instead.
    cell_verts = (np.asarray(el_V, float).ravel()
                  if n_el and len({len(V) for V in el_V}) == 1 else np.zeros(0))
    face_verts = (np.asarray(bf_V, float).ravel()
                  if n_bf and len({len(V) for V in bf_V}) == 1 else np.zeros(0))
    # POLYTOPE flat geometry for the C++ HEX/WEDGE charge-Gram H-matrix (the polytope _ChargeGramHMatrix
    # ctor): computed only when the mesh has non-tet cells / non-triangle faces (a pure tet mesh takes the
    # cell_verts/face_verts path, so the per-cell ConvexHull is skipped).
    need_poly = any(len(V) != 4 for V in el_V) or any(len(V) != 3 for V in bf_V)
    poly = (polytope_flat_geom(el_V, bf_V, el_vol, bf_area) if need_poly else None)
    return dict(M_mass=M_mass, B=B, ndof=ndof, m_unit=m_unit,
                cent=cent, meas=meas, self_energy=diagG, B_csr=B_sp, n_charge=n_el + n_bf,
                cell_verts=cell_verts, face_verts=face_verts, n_el=n_el, poly=poly)


def _charge_subpoints(mesh, nsub):
    """Sub-points + equal weights per charge (cells then boundary faces), element-type aware:
    tet cells use the (fixed) barycentric tet lattice, boundary tri faces the triangle lattice."""
    el_V = [np.array([mesh[v].point for v in el.vertices]) for el in mesh.Elements(ng.VOL)]
    bf_V = [np.array([mesh[v].point for v in el.vertices]) for el in mesh.Elements(ng.BND)]
    lam_t, lam_f = _bary_tet(nsub), _bary_tri(nsub)
    SP, SW, kind = [], [], []
    for V in el_V:
        vol = abs(np.linalg.det(V[1:] - V[0])) / 6.0
        SP.append(lam_t @ V); SW.append(np.full(len(lam_t), vol / len(lam_t))); kind.append(0)
    for V in bf_V:
        ar = 0.5 * np.linalg.norm(np.cross(V[1] - V[0], V[2] - V[0]))
        SP.append(lam_f @ V); SW.append(np.full(len(lam_f), ar / len(lam_f))); kind.append(1)
    return SP, SW, kind


def build_near_correction(mesh, d, nsub=4, near_factor=2.0, skip_surface_surface=False):
    """Sparse near-field Gram correction (exact sub-point MINUS centroid-monopole) for NEAR charge
    pairs -- the standard H-matrix near-field correction.  The scalable Gram is then the compressed
    monopole H-matrix (far, cheap) PLUS this sparse local correction (near, exact): it lifts the
    sphere demag from the monopole ~0.31 to the Gram-exact ~0.33 (-> analytic 1/3) at O(N) extra cost.
    Cells contribute ~0 for uniform M (div M = 0), but the correction is computed for all near pairs
    for generality.  Returns a scipy CSR (n_charge x n_charge).

    skip_surface_surface=True omits surface-surface pairs (for callers whose surface-surface block is
    already exact); the near-correction then only fixes the VOLUME-involving (cell-cell, cell-face) near
    pairs that the per-element NONLINEAR Newton needs (without it the volume near-field is un-corrected
    and Newton finds a wrong root)."""
    import scipy.sparse as sp
    inv4pi = 1.0 / (4.0 * pi)
    cent, meas = d["cent"], d["meas"]
    n_cell = d["n_charge"] - sum(1 for _ in mesh.Elements(ng.BND))
    SP, SW, _ = _charge_subpoints(mesh, nsub)
    n = len(cent)
    size = np.array([meas[a] ** (1.0 / 3.0) if a < n_cell else meas[a] ** 0.5 for a in range(n)])
    rows, cols, vals = [], [], []
    for a in range(n):
        ca, sa, wa = cent[a], SP[a], SW[a]
        for b in range(a + 1, n):
            if skip_surface_surface and a >= n_cell and b >= n_cell:
                continue                         # caller handles the surface-surface block exactly
            dx = ca - cent[b]
            r = float(np.sqrt(dx @ dx))
            if r < near_factor * (size[a] + size[b]):
                D = np.linalg.norm(sa[:, None, :] - SP[b][None, :, :], axis=2)
                exact = float(np.sum(np.outer(wa, SW[b]) * inv4pi / D))
                mono = meas[a] * meas[b] * inv4pi / r
                delta = exact - mono
                rows += [a, b]; cols += [b, a]; vals += [delta, delta]
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
