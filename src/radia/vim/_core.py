"""hdiv_demag_tet.py -- HDiv-type VIM demag on UNSTRUCTURED TET meshes (production step #1a).

The structured-hex prototype (hdiv_demag_quad_self.py) is generalized to real tet RT0 meshes from
NGSolve.  The NGSolve HDiv(order=0) extraction is element-AGNOSTIC:
  B = [ -div(u) tested on L2 (volume charge rho) ;  u.Trace().n tested on SurfaceL2 (surface charge sigma) ]
  M_mass = the HDiv(0) mass (the physical demag factors are eig(N, M_mass), basis-invariant).
Only the Gram self-energy geometry is element-specific -> here TET (barycentric sub-points) and
TRIANGLE faces, with the self constants c_tet / c_tri (the tet/tri analogs of c_cube / c_sq).

Verification (the physics gate): a tet-meshed SPHERE has demag factor EXACTLY 1/3 (isotropic);
a CUBE has the three axis demag factors equal to 1/3 each (sum to 1 by symmetry).  We compute the
demag factor of a UNIFORM M_z as the Rayleigh quotient (m^T N m)/(m^T M_mass m) with m = the RT0
L2-projection of (0,0,1).  Sphere -> ~1/3 confirms the unstructured tet operator is physical.

This is the Python reference; production step #1b moves the dense Gram G to the C++ HACApK charge
H-matrix (scalable) and checks it reproduces this.
"""
import json
import os
from math import pi, sqrt, log, atan2

import numpy as np
import scipy.sparse as sp

import ngsolve as ng
from netgen.csg import CSGeometry, Sphere, Pnt, OrthoBrick

ng.SetNumThreads(4)
HERE = os.path.dirname(os.path.abspath(__file__))
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
    (phi_tet/tri_potential) this converges to ~machine precision -- it REPLACES the crude equal-weight
    _bary_tet outer rule, which under-integrated the volume self-energy by ~6.5% (golden-invisible because
    every uniform-M demag has div M = 0).  Same rule as radia.vim._vim._tet_ref(o)."""
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
    volume-self quadrature defect that was fixed in the C++ analytic ctor (the production order-0 Gram) and
    in analytic_charge_gram (the dense reference).  It is NOT fixed here because it feeds only the dense-path
    PRECONDITIONER diagonal + the monopole-path Gram diagonal, never the production scalable operator (which
    uses the C++ ctor's own QuadDot(a,a)).  Correcting it makes the wilton_surface=True non-uniform solve
    CONVERGE (to ~0.6% of analytic_gram) instead of diverging, which retires the surface-only fail-loud
    guard (test_hdiv_vim_cyoke_nonlinear) -- a behaviour change deferred for a separate decision."""
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


# --- Wilton analytic surface Gram (the #1 production accuracy fix, 2026-06-07) ---
# For a UNIFORM M the demag is ENTIRELY surface charge (div M = 0 -> zero volume charge), so the demag
# factor is governed by the surface-surface (boundary-triangle) Gram block.  The centroid-MONOPOLE
# off-diagonal under-resolves adjacent/near boundary triangles -> the demag factor comes out ~5-6% low
# (cube 0.311, sphere 0.314 vs the exact 1/3), and the sub-point near-correction does NOT fix the cube.
# Wilton's exact analytic potential of a uniformly-charged flat triangle closes this: cube/sphere demag
# -> 1/3 to <0.15%.  Reference: Wilton et al., IEEE TAP 32(3):276 (1984); Graglia, IEEE TAP 41(10):1448
# (1993).  Off-diagonal only -- the validated tri_self_energy keeps the diagonal.

# Dunavant degree-5 symmetric triangle rule (7 points, barycentric; weights sum to 1) for the OUTER
# integral over the observation triangle (the inner integral is the exact Wilton potential).
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


def wilton_surface_block(bf_V):
    """Surface-surface Gram block (n_bf x n_bf) by the Wilton analytic inner integral + a Dunavant
    outer rule: G[a][b] = (1/4pi) INT_{tri_a} (INT_{tri_b} 1/|r-r'| dA') dA.  Symmetric; the diagonal
    here is the Dunavant-outer self (shape-exact for any triangle).  Vectorized over the observation
    points: for each SOURCE triangle b, evaluate its Wilton potential at ALL outer quad points at once
    -> the Python loop is O(n_bf), not O(n_bf^2)."""
    nb = len(bf_V)
    area = np.array([0.5 * np.linalg.norm(np.cross(V[1] - V[0], V[2] - V[0])) for V in bf_V])
    QP = np.array([_DUN5[:, :3] @ V for V in bf_V])  # (nb, 7, 3) outer quad points per triangle
    allP = QP.reshape(-1, 3)                          # (nb*7, 3)
    wq = _DUN5[:, 3]                                  # (7,)
    G = np.zeros((nb, nb))
    for b in range(nb):
        phi = tri_potential(bf_V[b], allP).reshape(nb, 7)   # potential of source b at every outer point
        G[:, b] = (phi * wq[None, :]).sum(axis=1) * area / (4 * pi)
    return 0.5 * (G + G.T)


def phi_tet(V, P):
    """Newtonian potential INT_tet 1/|P - r'| dV' of a uniform tetrahedron (vertices V, 4x3) at P --
    the VOLUME analog of tri_potential, via the divergence theorem (nabla'^2 R = 2/R):
        INT_V 1/R dV = (1/2) sum_{4 faces} d_face * INT_face 1/R dA'
    with d_face = (r' - P).n_hat (constant on a flat face) and the inner face integral = tri_potential.
    Reuses the exact Wilton triangle potential -> the volume Gram needs NO new singular quadrature.
    P: (3,) -> float, or (M,3) -> (M,).  Verified vs fine volume quadrature to the reference's accuracy."""
    P = np.asarray(P, float)
    single = (P.ndim == 1)
    R = P.reshape(1, 3) if single else P
    cen = V.mean(0)
    tot = np.zeros(len(R))
    for f in ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)):
        Fv = V[list(f)]
        n = np.cross(Fv[1] - Fv[0], Fv[2] - Fv[0])
        n = n / np.linalg.norm(n)
        if np.dot(Fv.mean(0) - cen, n) < 0:
            n = -n                                          # outward face normal
        d = (Fv[0] - R) @ n                                 # (M,) signed distance P -> face plane
        tot += d * tri_potential(Fv, R)
    tot *= 0.5
    return float(tot[0]) if single else tot


def analytic_charge_gram(el_V, bf_V, el_vol, bf_area, nsub_tet=3):
    """FULL analytic charge Gram (n_charge x n_charge, cells then boundary faces): every block via the
    EXACT analytic potential -- cell sources via phi_tet, face sources via tri_potential -- integrated
    against an outer quadrature on the target cell (tet sub-points) / face (Dunavant).  Replaces the
    centroid-monopole G everywhere, so the vol-vol and vol-surf (cell-cell, cell-face) blocks -- which
    matter for NON-uniform M (nonlinear, div M != 0) -- are analytic, not monopole+sub-point.  This is
    the #B 'volume Gram' that closes the residual sharp-body non-uniform gap (cube ~8.7%, C-yoke ~4%).
    Symmetric; O(n_charge) Python loop (vectorized over target points)."""
    n_el, n_bf = len(el_V), len(bf_V)
    n = n_el + n_bf
    lam_t, w_t = _gauss_duffy_tet(4)                         # Gauss-Duffy tet outer rule (nsub_tet now unused:
    QP, QW = [], []                                          # the old _bary_tet outer under-integrated vol-self ~6.5%)
    for k, V in enumerate(el_V):
        P = V[0] + lam_t @ (V[1:] - V[0])                    # physical outer pts (V0 + l1 e1 + l2 e2 + l3 e3)
        QP.append(P); QW.append(w_t * 6.0 * el_vol[k])       # phys weight = w_ref * |J|, |J| = 6*vol
    for k, V in enumerate(bf_V):
        QP.append(_DUN5[:, :3] @ V); QW.append(_DUN5[:, 3] * bf_area[k])
    allP = np.vstack(QP)
    allW = np.concatenate(QW)
    tidx = np.concatenate([np.full(len(QP[i]), i) for i in range(n)])
    inv4pi = 1.0 / (4.0 * pi)
    G = np.zeros((n, n))
    for s in range(n):
        pot = phi_tet(el_V[s], allP) if s < n_el else tri_potential(bf_V[s - n_el], allP)
        G[:, s] = np.bincount(tidx, weights=allW * pot * inv4pi, minlength=n)
    return 0.5 * (G + G.T)


def _csr_sp(bf):
    """NGSolve BilinearForm.mat -> scipy CSR (SPARSE). The scalable (skip_dense_gram) build keeps the
    charge maps + HDiv mass sparse; only the dense reference path densifies via _csr()."""
    m = bf.mat
    r, c, v = m.COO()
    return sp.csr_matrix((np.array(v), (np.array(r), np.array(c))), shape=(m.height, m.width))


def _csr(bf):
    return _csr_sp(bf).toarray()


def build_demag(mesh, nsub=4, wilton_surface=False, analytic_gram=False, skip_dense_gram=False):
    """Assemble the HDiv-type VIM demag operator N = B^T G B on a tet mesh + the HDiv mass M_mass.
    Returns N, M_mass, the charge map B, the loop basis, and diagnostics.

    wilton_surface=True replaces the surface-surface (boundary-triangle) OFF-diagonal Gram block with
    the exact Wilton analytic integral (the diagonal keeps the validated tri_self_energy).  This makes
    the demag factor exact to <0.15% on cube AND sphere (vs the ~5-6% monopole error that the sub-point
    near-correction cannot fix on the cube) -- the #1 production accuracy fix.  O(n_bf^2) dense; the
    scalable C++ HACApK charge-Gram path is the next step."""
    # NGSolve assembly (FES, charge-map BilinearForms, HDiv mass).  Per the caller-wraps policy this HELPER
    # does NOT open a TaskManager; the CALLER wraps build_demag in `with ng.TaskManager():` (CLAUDE.md
    # "TaskManager Wrap Policy: Caller Wraps, Helper Does NOT").
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
    Bv_sp, Bb_sp = _csr_sp(bv), _csr_sp(bb)   # charge maps as sparse CSR (densified only on the ref path)
    massv = ng.BilinearForm(L2v); massv += L2v.TrialFunction() * L2v.TestFunction() * ng.dx; massv.Assemble()
    massb = ng.BilinearForm(L2b); massb += L2b.TrialFunction() * L2b.TestFunction() * ng.ds; massb.Assemble()
    el_vol = _csr_sp(massv).diagonal(); bf_area = _csr_sp(massb).diagonal()   # L2-0 mass is diagonal == measures
    # HDiv mass (the physical demag-factor metric) -- kept SPARSE; the skip_dense_gram path never densifies it
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
    # scaled charge map B (each row / its measure) as sparse CSR -- the single source for both B_csr and,
    # on the dense reference path, the dense B.  (sparse: each charge cell touches only its own faces.)
    _Brows = [sp.diags(1.0 / el_vol) @ Bv_sp]
    if n_bf:
        _Brows.append(sp.diags(1.0 / bf_area) @ Bb_sp)
    B_sp = sp.vstack(_Brows).tocsr()

    # diagonal charge self-energies (O(N); kept even on the scalable path, for the preconditioner)
    diagG = np.empty(n_el + n_bf)
    for k, V in enumerate(el_V):
        diagG[k] = tet_self_energy(V, el_vol[k], nsub)
    for k, V in enumerate(bf_V):
        diagG[n_el + k] = tri_self_energy(V, bf_area[k], nsub)
    if skip_dense_gram:
        # SCALABLE path: SKIP the O(N^2) dense G/N and the O(N^2) loop SVD (the biggest dense costs;
        # G/N are n_charge^2).  The C++ analytic _ChargeGramHMatrix + SolveMaterialMINRES use
        # cent / meas / cell_verts / face_verts / diagG + the sparse B_csr only -- never the dense G.
        # M_mass + B are kept SPARSE here (no n_charge^2 dense), so this path has NO dense N^2 object at
        # all -- the analytic-Gram BUILD needs only the verts (straight from the mesh).
        G = N = loops = None
        n_loop = None
        M_mass = M_mass_sp        # sparse CSR (RT0 HDiv mass is local -> sparse)
        B = B_sp                  # sparse CSR scaled charge map
    else:
        M_mass = M_mass_sp.toarray()                       # dense reference path: small N
        B = B_sp.toarray()
        Dd = np.linalg.norm(cent[:, None, :] - cent[None, :, :], axis=2); np.fill_diagonal(Dd, np.inf)
        G = (meas[:, None] * meas[None, :]) / (4 * pi * Dd)
        np.fill_diagonal(G, diagG)
        if analytic_gram:
            # FULL analytic charge Gram (#B volume Gram): every block (vol-vol, vol-surf, surf-surf) via
            # the exact phi_tet / tri_potential potentials -- closes the NON-uniform-M (nonlinear) gap.
            G = analytic_charge_gram(el_V, bf_V, el_vol, bf_area)
        elif wilton_surface and n_bf > 0:
            # exact Wilton surface-surface block (shape-exact for ANY triangle; demag 1/3 to <0.15%).
            G[n_el:, n_el:] = wilton_surface_block(bf_V)
        N = B.T @ G @ B
        Q = np.vstack([Bv_sp.toarray(), Bb_sp.toarray()]) if n_bf else Bv_sp.toarray()
        sv = np.linalg.svd(Q, compute_uv=False)
        rankQ = int(np.sum(sv > 1e-9 * sv.max()))
        n_loop = ndof - rankQ
        _, _, Vt = np.linalg.svd(Q)
        loops = Vt[rankQ:, :]
    # charge geometry (for the C++ HACApK charge-Gram H-matrix path, #1b): centroids, measures,
    # diagonal self-energies, the dense Gram G, and the sparse charge map B as scipy CSR.
    # cell_verts [n_el*12] (tets, 4 verts) + face_verts [n_bf*9] (tris, 3 verts) feed the ANALYTIC
    # C++ charge Gram (M2b: _ChargeGramHMatrix analytic mode == this dense analytic_gram path).
    cell_verts = np.asarray(el_V, float).ravel() if n_el else np.zeros(0)
    face_verts = np.asarray(bf_V, float).ravel() if n_bf else np.zeros(0)
    return dict(N=N, M_mass=M_mass, B=B, ndof=ndof, n_loop=n_loop, loops=loops, m_unit=m_unit,
                cent=cent, meas=meas, self_energy=diagG, G=G,
                B_csr=B_sp, n_charge=n_el + n_bf,
                cell_verts=cell_verts, face_verts=face_verts, n_el=n_el)


def demag_factor(d):
    m = d["m_unit"]
    return float((m @ d["N"] @ m) / (m @ d["M_mass"] @ m))


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

    skip_surface_surface=True omits surface-surface pairs -- use this together with the Wilton surface
    Gram (build_demag(wilton_surface=True)), which already makes the surface-surface block exact; the
    near-correction then only fixes the VOLUME-involving (cell-cell, cell-face) near pairs that the
    per-element NONLINEAR Newton needs (without it the volume near-field is un-corrected and Newton
    finds a wrong root)."""
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
                continue                         # surface-surface handled exactly by the Wilton Gram
            dx = ca - cent[b]
            r = float(np.sqrt(dx @ dx))
            if r < near_factor * (size[a] + size[b]):
                D = np.linalg.norm(sa[:, None, :] - SP[b][None, :, :], axis=2)
                exact = float(np.sum(np.outer(wa, SW[b]) * inv4pi / D))
                mono = meas[a] * meas[b] * inv4pi / r
                delta = exact - mono
                rows += [a, b]; cols += [b, a]; vals += [delta, delta]
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))


def report(mesh, tag, nsub=4):
    d = build_demag(mesh, nsub)
    Nn = np.linalg.norm(d["N"], 2)
    asym = np.linalg.norm(d["N"] - d["N"].T) / Nn
    loop_res = (max(np.linalg.norm(d["N"] @ d["loops"][k]) for k in range(d["n_loop"])) / Nn
                if d["n_loop"] else 0.0)
    Dz = demag_factor(d)
    print(f"[{tag}] ndof={d['ndof']} n_loop={d['n_loop']} asym={asym:.1e} loop_res={loop_res:.1e} "
          f"demag_z={Dz:.4f}")
    return {"tag": tag, "ndof": d["ndof"], "n_loop": d["n_loop"], "asym": asym,
            "loop_res": loop_res, "demag_z": Dz}


if __name__ == "__main__":
    res = {}
    # cube (tet-meshed): each axis demag factor ~ 1/3
    geo = CSGeometry(); geo.Add(OrthoBrick(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5)))
    res["cube"] = report(ng.Mesh(geo.GenerateMesh(maxh=0.25)), "cube tet", nsub=4)
    # sphere (tet-meshed): demag factor EXACTLY 1/3
    for h in (0.5, 0.35, 0.25):
        geo = CSGeometry(); geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
        res[f"sphere_h{h}"] = report(ng.Mesh(geo.GenerateMesh(maxh=h)), f"sphere tet h={h}", nsub=4)
    with open(os.path.join(HERE, "hdiv_demag_tet.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("saved", os.path.join(HERE, "hdiv_demag_tet.json"))
