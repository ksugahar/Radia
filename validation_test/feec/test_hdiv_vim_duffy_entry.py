"""Golden: the C++ Duffy singular-quadrature charge potential (RadHACApKChargeGram::PhiAtHO_Duffy) -- the
order>=3 / curved-panel path that takes over from the analytic-moment kernels when they run out (a tet
volume charge of degree>=2 needs TetMoment2; a surface charge of degree>=3 needs degree-3 moments).

PhiInner dispatches: charge degree <= 2 (tet deg<=1, face deg<=2) -> the EXACT analytic moment potential;
else -> this Duffy quadrature (6-pt Gauss-Legendre on signed radial sub-tets / sub-triangles from x0 =
closest point of the host to the field point).  Locks the C++ Duffy entries (degree-2 tet charges, the
order-3 volume-charge regime) to an INDEPENDENT brute force at ~5e-3 -- including the SELF block and the
sign on a NEGATIVELY-oriented tet (the signed sub-tets give the signed-volume integral; PhiAtHO_Duffy
multiplies by sign(host vol) to recover the physical absolute-volume charge integral).

NOTE: the Duffy is ~1e-3 accurate -- enough for curved-panel field evaluation, but NOT for the order>=3
MATERIAL solve, where the ill-conditioned high-degree monomial basis (cond(B)^2 in N=B^T G B) amplifies the
~1e-3 entry error so the demag spectrum escapes [0,1].  So hdiv_demag_solve(order>2) stays fail-loud; this
test locks the Duffy ENTRY accuracy (the curved / future path), not an order>=3 material solve.
See memory hdiv-vim-sauter-schwab-cg.
"""
import numpy as np
import pytest

rp = pytest.importorskip("radia._radia_pybind")
INV4PI = 1.0 / (4.0 * np.pi)


def _g01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w


def _tet_ref(o):
    s, ws = _g01(o)
    P, W = [], []
    for a, wa in zip(s, ws):
        for b, wb in zip(s, ws):
            for c, wc in zip(s, ws):
                P.append((a, b * (1 - a), c * (1 - a) * (1 - b)))
                W.append(wa * wb * wc * (1 - a) ** 2 * (1 - b))
    return np.array(P), np.array(W)


def _tri_ref(o):
    s, ws = _g01(o)
    P, W = [], []
    for u, wu in zip(s, ws):
        for v, wv in zip(s, ws):
            P.append((u, v * (1 - u)))
            W.append(wu * wv * (1 - u))
    return np.array(P), np.array(W)


def _affine_tet(V, y):
    A = np.array([V[1] - V[0], V[2] - V[0], V[3] - V[0]]).T
    return np.linalg.solve(A, y - V[0])


def _cloud_abs(V, n):
    g, w = _g01(n)
    D = abs(np.dot(np.cross(V[1] - V[0], V[2] - V[0]), V[3] - V[0]))
    e10, e21, e32 = V[1] - V[0], V[2] - V[1], V[3] - V[2]
    P, W, L = [], [], []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                u, v, ww = g[i], g[j], g[k]
                y = V[0] + u * (e10 + v * (e21 + ww * e32))
                P.append(y); W.append(w[i] * w[j] * w[k] * u * u * v * D); L.append(_affine_tet(V, y))
    return np.array(P), np.array(W), np.array(L)


def _brute(Va, ea, Vb, eb, na, nb):
    Pa, Wa, La = _cloud_abs(Va, na)
    Pb, Wb, Lb = _cloud_abs(Vb, nb)
    mb = Lb[:, 0] ** eb[0] * Lb[:, 1] ** eb[1] * Lb[:, 2] ** eb[2]
    s = 0.0
    for i in range(len(Pa)):
        ma = La[i, 0] ** ea[0] * La[i, 1] ** ea[1] * La[i, 2] ** ea[2]
        r = np.linalg.norm(Pa[i] - Pb, axis=1)
        s += Wa[i] * ma * np.sum(Wb * mb / r)
    return s * INV4PI


def test_cpp_duffy_degree2_entries_match_brute():
    """C++ Duffy (PhiAtHO_Duffy via PhiInner) on degree-2 tet volume charges == brute to ~5e-3, including the
    self block and the negatively-oriented-tet sign."""
    Va = np.array([[0., 0., 0.], [1., 0., 0.], [0.2, 1.0, 0.], [0.3, 0.4, 0.9]])   # +oriented
    Vb = np.array([[0., 0., 0.], [1., 0., 0.], [0.2, 1.0, 0.], [0.4, 0.3, -0.8]])  # -oriented (apex below)
    host = [0, 1]
    kind = [0, 0]
    expo = [2, 0, 0, 2, 0, 0]                                  # l0^2 on each tet (degree 2 -> Duffy path)
    exps = [(2, 0, 0), (2, 0, 0)]
    rtp, rtw = _tet_ref(8)
    rsp, rsw = _tri_ref(8)
    o = rp._ChargeGramHMatrix(
        cell_verts=np.concatenate([Va.ravel(), Vb.ravel()]).tolist(), face_verts=[], n_el=2,
        charge_host=host, charge_kind=kind, charge_expo=expo,
        ref_tet_pts=rtp.ravel().tolist(), ref_tet_w=rtw.tolist(),
        ref_tri_pts=rsp.ravel().tolist(), ref_tri_w=rsw.tolist(), build=False)
    for (a, b) in [(0, 0), (0, 1), (1, 1)]:
        cpp = o.entry(int(a), int(b))
        ref = _brute(Va if host[a] == 0 else Vb, exps[a], Va if host[b] == 0 else Vb, exps[b], 22, 24)
        rel = abs(cpp - ref) / (abs(ref) + 1e-30)
        assert rel < 5e-3, f"Duffy entry({a},{b}) C++ {cpp:.6e} vs brute {ref:.6e} rel {rel:.2e}"
        assert np.sign(cpp) == np.sign(ref), f"Duffy entry({a},{b}) sign wrong (host orientation): {cpp:.3e} vs {ref:.3e}"
