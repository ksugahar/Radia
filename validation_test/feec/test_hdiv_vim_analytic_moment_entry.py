"""Golden: the C++ high-order charge-Gram NEAR entry uses the EXACT analytic moment potential
(`RadHACApKChargeGram::PhiAtHO_Analytic`), NOT the old point-subtraction `PhiAtHO`.

The old subtraction evaluated the source monomial at the TARGET outer point (outside the source host for
ADJACENT pairs) -> ~1.5e-2 wrong on adjacent high-order entries (the M4 entry bug).  The analytic-moment
path (affine monomial -> physical polynomial A + B.y + y^T C y, contracted with PhiTet/TetMoment1 and
TriPotential/TriMoment1/TriMoment2) is EXACT through the singularity for self/adjacent/far alike.

This locks the C++ adjacent tet-tet entry to an INDEPENDENT brute force (no shared formula) at ~1e-3 -- a
band the old subtraction (~1.5e-2) would FAIL.  (NOTE: this is the ENTRY accuracy fix; the order-p MATERIAL
SOLVE remains a separate open item -- the high-order demag OPERATOR N=B^T G B is invalid, eig escapes [0,1];
see memory hdiv-highorder-material-solve-wrong.  So Solve(order>0) stays fail-loud.)
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


def _tet_cloud(V, n):
    """absolute-volume (physical) tet Gauss cloud via the collapsed-cube map; monomials in affine coords."""
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


def _brute_tet_tet(Va, ea, Vb, eb, na, nb):
    Pa, Wa, La = _tet_cloud(Va, na)
    Pb, Wb, Lb = _tet_cloud(Vb, nb)
    mb = Lb[:, 0] ** eb[0] * Lb[:, 1] ** eb[1] * Lb[:, 2] ** eb[2]
    s = 0.0
    for i in range(len(Pa)):
        ma = La[i, 0] ** ea[0] * La[i, 1] ** ea[1] * La[i, 2] ** ea[2]
        r = np.linalg.norm(Pa[i] - Pb, axis=1)
        s += Wa[i] * ma * np.sum(Wb * mb / r)
    return s * INV4PI


def test_cpp_adjacent_highorder_entry_matches_brute():
    """C++ oracle.entry for ADJACENT high-order tet-tet pairs == brute to ~1e-3 (old subtraction was ~1.5e-2)."""
    # two tets sharing face (V0,V1,V2)
    Va = np.array([[0., 0., 0.], [1., 0., 0.], [0.2, 1.0, 0.], [0.3, 0.4, 0.9]])
    Vb = np.array([[0., 0., 0.], [1., 0., 0.], [0.2, 1.0, 0.], [0.4, 0.3, -0.8]])
    host = [0, 0, 1, 1]
    kind = [0, 0, 0, 0]
    expo = [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0]            # const + linear l0 on each tet
    exps = [(0, 0, 0), (1, 0, 0), (0, 0, 0), (1, 0, 0)]
    rtp, rtw = _tet_ref(7)
    rsp, rsw = _tri_ref(7)
    oracle = rp._ChargeGramHMatrix(
        cell_verts=np.concatenate([Va.ravel(), Vb.ravel()]).tolist(),
        face_verts=[], n_el=2, charge_host=host, charge_kind=kind, charge_expo=expo,
        ref_tet_pts=rtp.ravel().tolist(), ref_tet_w=rtw.tolist(),
        ref_tri_pts=rsp.ravel().tolist(), ref_tri_w=rsw.tolist(), build=False)
    worst = 0.0
    for (a, b) in [(0, 2), (1, 2), (0, 3), (1, 3)]:
        cpp = oracle.entry(int(a), int(b))
        Vt_a = Va if host[a] == 0 else Vb
        Vt_b = Va if host[b] == 0 else Vb
        ref = _brute_tet_tet(Vt_a, exps[a], Vt_b, exps[b], 22, 24)
        rel = abs(cpp - ref) / (abs(ref) + 1e-30)
        worst = max(worst, rel)
        assert rel < 2e-3, f"entry({a},{b}) C++ {cpp:.6e} vs brute {ref:.6e} rel {rel:.2e} -- regressed to subtraction?"
    assert worst < 2e-3, f"worst adjacent high-order entry rel {worst:.2e} (old subtraction was ~1.5e-2)"
