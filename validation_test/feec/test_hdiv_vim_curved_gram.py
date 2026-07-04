"""Golden: the PRODUCTION C++ CURVED charge Gram (RadHACApKChargeGram curved ctor, exposed via
radia.vim.ChargeGram(curve_order=2) -> _ChargeGramHMatrix(cell_nodes=.., curve_order=2, ..)) and the
end-to-end curved demag solve Solve(curve_order=2).

This is the piece-3 WIRING of the curved-geometry HDiv-VIM charge Gram (the singular-quadrature CORE is the
already-golden-locked rad_hdiv::CurvedTri/TetPotential, test_hdiv_vim_curved_{tri,tet}.py).  The production
H-matrix curved Gram is validated against an INDEPENDENT dense build of the SAME curved math (curved outer
quad via the P2 map + curved measure x the curved inner potential), on representative charge pairs.

ROLE (de-risked 2026-06-28, memory hdiv-vim-sauter-schwab-cg): curved helps NEAR-SURFACE FIELD / FLUX
accuracy (sigma=M.n on the TRUE curved surface), NOT the volume-averaged demag FACTOR (curving-insensitive
~3e-5 on a sphere -- see also test_hdiv_vim_curved.py which measures the win vs the exact dipole FIELD).  So
the GRAM-correctness + curved-solve-runs goldens lock the math, NOT a demag-factor improvement.

UPDATE 2026-06-29 (Sugahara, memory hdiv-vim-tet-rt1-only METRIC FLIP): test_curved_moment_beats_flat ALSO
locks curving's HEADLINE engineering win -- the MAGNETIC MOMENT m = INT M dV (the dipole far-field quantity),
which the demag factor (volume-insensitive), per-element leak (RT1+Piola noise), and M_avg (volume divided
out) ALL hide.  A faceted P1 sphere under-counts the volume ~9% -> moment ~9% low; P2 curving recovers it to
<0.2% -> moment ~0.13% (flat-vs-curved >10x vs the analytic linear-sphere dipole).
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
rp = pytest.importorskip("radia._radia_pybind")
from netgen.occ import Sphere, Pnt, OCCGeometry          # noqa: E402
from radia.vim._vim import _charge_basis_curved, _tet_ref, _tri_ref, _g01  # noqa: E402
from radia.vim import Solve, ChargeGram                    # noqa: E402

INV4PI = 1.0 / (4.0 * np.pi)


def _Ntet(xi, eta, ze):
    L1, L2, L3, L4 = 1-xi-eta-ze, xi, eta, ze
    return np.array([L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1), L4*(2*L4-1),
                     4*L1*L2, 4*L2*L3, 4*L3*L1, 4*L1*L4, 4*L2*L4, 4*L3*L4])


def _dNtet(xi, eta, ze):
    L1, L2, L3, L4 = 1-xi-eta-ze, xi, eta, ze
    dL = [np.array([-1., -1, -1]), np.array([1., 0, 0]), np.array([0., 1, 0]), np.array([0., 0, 1])]
    d = np.zeros((10, 3))
    for k in range(3):
        d[0, k] = (4*L1-1)*dL[0][k]; d[1, k] = (4*L2-1)*dL[1][k]; d[2, k] = (4*L3-1)*dL[2][k]; d[3, k] = (4*L4-1)*dL[3][k]
        d[4, k] = 4*(dL[0][k]*L2+L1*dL[1][k]); d[5, k] = 4*(dL[1][k]*L3+L2*dL[2][k]); d[6, k] = 4*(dL[2][k]*L1+L3*dL[0][k])
        d[7, k] = 4*(dL[0][k]*L4+L1*dL[3][k]); d[8, k] = 4*(dL[1][k]*L4+L2*dL[3][k]); d[9, k] = 4*(dL[2][k]*L4+L3*dL[3][k])
    return d


def _Ntri(xi, eta):
    L1, L2, L3 = 1-xi-eta, xi, eta
    return np.array([L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1), 4*L1*L2, 4*L2*L3, 4*L3*L1])


def _dNtri(xi, eta):
    L1 = 1-xi-eta
    dxi = np.array([1-4*L1, 4*xi-1, 0.0, 4*(L1-xi), 4*eta, -4*eta])
    det = np.array([1-4*L1, 0.0, 4*eta-1, -4*xi, 4*xi, 4*(L1-eta)])
    return dxi, det


def _outer_quad(cb, a, rtp, rtw, rsp, rsw, cn, fn):
    """curved outer quad of charge a: physical points X(xi_q) + monomial-folded curved-measure weights."""
    host, kind, expo = cb["host"][a], cb["kind"][a], np.array(cb["expo"][3*a:3*a+3])
    Xs, Ws = [], []
    if kind == 0:
        nd = cn[host]
        for (xi, eta, ze), w in zip(rtp, rtw):
            Xs.append(_Ntet(xi, eta, ze) @ nd)
            Jv = abs(np.linalg.det((_dNtet(xi, eta, ze).T @ nd).T))
            Ws.append(w * Jv * xi**expo[0] * eta**expo[1] * ze**expo[2])
    else:
        nd = fn[host]
        for (xi, eta), w in zip(rsp, rsw):
            Xs.append(_Ntri(xi, eta) @ nd)
            dxi, det = _dNtri(xi, eta)
            Jv = np.linalg.norm(np.cross(dxi @ nd, det @ nd))
            Ws.append(w * Jv * xi**expo[0] * eta**expo[1])
    return np.array(Xs), np.array(Ws)


def _phi(cb, b, p, gl, gw, cn, fn):
    host, kind, e = cb["host"][b], cb["kind"][b], np.array(cb["expo"][3*b:3*b+3])
    if kind == 0:
        return rp._hdiv_curved_tet_potential(cn[host].ravel().tolist(), int(e[0]), int(e[1]), int(e[2]), list(p), gl, gw)
    return rp._hdiv_curved_tri_potential(fn[host].ravel().tolist(), int(e[0]), int(e[1]), list(p), gl, gw)


def _quaddot(cb, a, b, rtp, rtw, rsp, rsw, gl, gw, cn, fn):
    """(1/4pi) INT_a m_a Phi_b -- the curved outer quad of a times the curved potential of b (one direction)."""
    Xa, Wa = _outer_quad(cb, a, rtp, rtw, rsp, rsw, cn, fn)
    return INV4PI * sum(Wa[q] * _phi(cb, b, Xa[q], gl, gw, cn, fn) for q in range(len(Wa)))


def _dense_entry(cb, a, b, rtp, rtw, rsp, rsw, gl, gw, cn, fn):
    """The SYMMETRIZED entry 0.5*(QuadDot(a,b)+QuadDot(b,a)) -- matches GetInteractionMatrixElement's high-order
    symmetrization (the curved outer quadrature is not exactly symmetric, so the single direction differs by
    ~1e-4 on far pairs; the production H-matrix .entry() symmetrizes, so the dense reference must too)."""
    qab = _quaddot(cb, a, b, rtp, rtw, rsp, rsw, gl, gw, cn, fn)
    if a == b:
        return qab
    qba = _quaddot(cb, b, a, rtp, rtw, rsp, rsw, gl, gw, cn, fn)
    return 0.5 * (qab + qba)


def test_curved_gram_entries_match_dense_curved_math():
    """The production curved H-matrix Gram reproduces the independent dense curved math (== the validated
    curved potentials) to machine precision on representative self/near/far pairs across all charge-kind
    combinations."""
    geo = OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0))
    mesh = ng.Mesh(geo.GenerateMesh(maxh=1.5)); mesh.Curve(2)
    order = 1; quad = 3 * order; cg = 8     # RT1 (RT0 retired); quad = the LINEAR ChargeGram default (3*p
    #                                          PSD floor) -- the dense reference MUST use the SAME quad as the
    #                                          production Gram here (ChargeGram with no bh_table -> linear)
    gx, gw = _g01(cg); gl, gwl = gx.tolist(), gw.tolist()
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=order)
        B, G, M_mass = ChargeGram(fes, curve_order=2, curve_gauss=cg, eps=1e-12, leafsize=8)
        cb = _charge_basis_curved(fes, quad)
    n = len(cb["host"]); n_el = cb["n_el"]
    cn = np.array(cb["cell_nodes"]).reshape(-1, 10, 3)
    fn = np.array(cb["face_nodes"]).reshape(-1, 6, 3)
    rtp, rtw = _tet_ref(quad); rsp, rsw = _tri_ref(quad)
    kind = np.array(cb["kind"]); cent = np.array(cb["host"], float)  # placeholder; use centroids below
    cent = np.array([cn[cb["host"][a]].mean(0) if cb["kind"][a] == 0 else fn[cb["host"][a]].mean(0)
                     for a in range(n)])
    cells = np.where(kind == 0)[0]; faces = np.where(kind == 1)[0]

    def nearest_far(group_a, group_b):
        ca, cb_ = cent[group_a], cent[group_b]
        D = np.linalg.norm(ca[:, None, :] - cb_[None, :, :], axis=2)
        D[D < 1e-12] = np.inf
        i, j = np.unravel_index(np.argmin(D), D.shape)
        i2, j2 = np.unravel_index(np.argmax(np.where(np.isinf(D), -1, D)), D.shape)
        return (group_a[i], group_b[j]), (group_a[i2], group_b[j2])

    pairs = [(int(cells[0]), int(cells[0])), (int(faces[0]), int(faces[0]))]  # two self pairs
    for ga, gb in [(cells, cells), (cells, faces), (faces, faces)]:
        (na, nb), (fa, fb) = nearest_far(ga, gb)
        pairs += [(int(na), int(nb)), (int(fa), int(fb))]

    worst = 0.0
    for (a, b) in pairs:
        gh = G.entry(a, b)
        gd = _dense_entry(cb, a, b, rtp, rtw, rsp, rsw, gl, gwl, cn, fn)
        rel = abs(gh - gd) / (abs(gd) + 1e-300)
        worst = max(worst, rel)
        assert rel < 1e-9, f"curved Gram entry ({a},{b}) kinds({cb['kind'][a]},{cb['kind'][b]}): " \
                           f"H-matrix {gh:.6e} vs dense {gd:.6e} rel {rel:.2e}"
    assert worst < 1e-9


def test_curved_demag_solve_runs_and_converges():
    """End-to-end: Solve(curve_order=2) RUNS, converges, and returns a physically sane demag
    factor (near 1/3 for a sphere -- NOT a curving-improvement claim; the demag factor is curving-insensitive,
    this just locks that the curved solve produces a valid result)."""
    geo = OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0))
    mesh = ng.Mesh(geo.GenerateMesh(maxh=1.0))
    Hext = ng.CoefficientFunction((0, 0, 1.0))
    with ng.TaskManager():
        r = Solve(mesh, mu_r=100.0, H_ext=Hext, order=1, curve_order=2)
    assert r["curve_order"] == 2
    assert r["iters"] < 400
    assert 0.20 < r["demag"] < 0.45, r["demag"]            # sphere demag ~1/3 (loose band; curving-insensitive)
    assert abs(r["M_avg"][2]) > 1.0                         # responds to the applied field
    assert abs(r["M_avg"][0]) < 0.05 and abs(r["M_avg"][1]) < 0.05   # transverse leak small (true sphere)


def test_curved_moment_beats_flat():
    """Curving's REAL engineering win = the MAGNETIC MOMENT m = INT M dV (the dipole FAR-FIELD quantity),
    NOT the volume-averaged demag factor (curving-insensitive) NOR the per-element leak (RT1 + curved-Piola
    FE noise that averages out).  A faceted P1 sphere INSCRIBES the surface (flat faces chord across) and
    UNDER-counts the volume by ~9% at this coarse mesh -> its moment is ~9% low -> the dipole far field is
    ~9% wrong.  P2 curving pushes the mid-edge nodes onto the sphere -> volume recovered to <0.2% -> the
    moment is accurate to ~0.13%.  This locks that curved RT1 beats flat by >10x on the moment vs the
    analytic linear-sphere dipole (m_z err ~= the volume deficit, since M_avg is curving-insensitive).
    Measured margins: flat ~9.3% / curved ~0.13% (70x).  (Sugahara 2026-06-29; C:\\temp\\moment_err.py.)"""
    R, mu_r, H0 = 1.0, 100.0, 1.0e4
    M_in = 3.0 * (mu_r - 1.0) / (mu_r + 2.0) * H0           # uniform internal M of a linear sphere
    m_exact = M_in * (4.0 / 3.0 * np.pi * R ** 3)           # z-only analytic dipole moment INT M dV

    def moment_err(curve):
        """Generate a fresh mesh, (optionally) P2-curve it, solve, and return |m_comp - m_exact| / |m_exact|.
        Each solve runs in its OWN TaskManager region (no cross-solve state; Netgen is deterministic so the
        flat and curved meshes share topology -- apples to apples)."""
        mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=0.6))
        Hext = ng.CoefficientFunction((0, 0, H0))
        with ng.TaskManager():
            if curve:
                mesh.Curve(2)
            r = Solve(mesh, mu_r=mu_r, H_ext=Hext, order=1,
                                 curve_order=2 if curve else None)
            V = float(ng.Integrate(ng.CoefficientFunction(1.0), mesh))   # mesh volume (curved when curved)
        m = np.asarray(r["M_avg"], float) * V                # magnetic moment = volume-averaged M * mesh volume
        return float(np.linalg.norm([m[0], m[1], m[2] - m_exact]) / m_exact)

    e_flat = moment_err(False)
    e_curved = moment_err(True)
    assert e_curved < 5.0e-3, "curved moment err %.4f%% (want <0.5%%)" % (100 * e_curved)
    assert e_curved < e_flat / 10.0, \
        "curving win too small: flat %.4f%% curved %.4f%% (want curved < flat/10)" % (100 * e_flat, 100 * e_curved)
