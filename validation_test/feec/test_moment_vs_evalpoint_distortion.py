"""Reproducible backing for the MMMM paper's central claim (SA-26-010): the multipole-MOMENT closure
is far more accurate than the classical MSC EVAL-POINT (EIEM2) closure on DISTORTED (sheared) elements.

Both closures use the SAME per-edge uniform line-charge DOF on the SAME sheared mesh; they differ ONLY
in how the charges are determined:
  * MOMENT    -- radia.mmmm2d (1 monopole + 2 dipole + (nEdge-3) residual-quadrupole conditions).
  * EVAL-POINT-- sigma_f = chi * H_n(evalpt_f), evalpt_f = 0.5(edge_mid_f + centroid): the EIEM2
    collocation (nEdge equations per element), implemented here with the SAME exact 2D segment field
    the moment closure uses (so the comparison is fair -- no quadrature bias).

A unit disk (demag D = 1/2) is sheared (x,y)->(x+s y, y); the sheared disk is exactly an ELLIPSE whose
2x2 demag tensor is analytic, giving a closed-form reference on increasingly-sheared elements.

This is the 2D, analytic-reference, reproducible counterpart of the paper's 3D hex Table (whose numbers
had no committed source).  Result (locked below): the moment closure is ORDERS OF MAGNITUDE more
accurate than eval-point across shear, because the residual quadrupole carries the shear gradient the
normal-only collocation loses.  Saves test_moment_vs_evalpoint_distortion_summary.json.
"""
import json
import os

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.geom2d import SplineGeometry
from netgen.meshing import Mesh as NGMesh, MeshPoint, Element2D, Element1D, FaceDescriptor
from netgen.csg import Pnt

import radia
import radia.mmmm2d as m2

MU0 = 4e-7 * np.pi
CHI = 20.0
H0 = (1.0, 0.0)
HERE = os.path.dirname(os.path.abspath(__file__))


def _base_disk(maxh, quad=True):
    g = SplineGeometry(); g.AddCircle((0, 0), r=1.0, bc="outer")
    return ng.Mesh(g.GenerateMesh(maxh=maxh, quad_dominated=quad))


def _remap(base, new_pts):
    m = NGMesh(dim=2); fd = m.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    pid = [m.Add(MeshPoint(Pnt(float(x), float(y), 0.0))) for x, y in new_pts]
    for el in base.Elements(ng.VOL):
        m.Add(Element2D(fd, [pid[v.nr] for v in el.vertices]))
    for el in base.Elements(ng.BND):
        vs = [pid[v.nr] for v in el.vertices]; m.Add(Element1D([vs[0], vs[1]], index=1))
    m.SetBCName(0, "outer"); return ng.Mesh(m)


def _sheared(base, s):
    p = np.array([list(base[v].point)[:2] for v in base.vertices])
    o = p.copy(); o[:, 0] = p[:, 0] + s * p[:, 1]
    return _remap(base, o)


def _analytic_M(s, chi, h0):
    Q = np.array([[1.0, -s], [-s, 1.0 + s * s]]); lam, V = np.linalg.eigh(Q)
    L = 1.0 / np.sqrt(lam); D = np.array([L[1] / (L[0] + L[1]), L[0] / (L[0] + L[1])])
    return np.linalg.solve(np.eye(2) + chi * (V @ np.diag(D) @ V.T), chi * np.asarray(h0, float))


def _elem_edges(mesh):
    pts = np.array([list(mesh[v].point)[:2] for v in mesh.vertices])
    out = []
    for el in mesh.Elements(ng.VOL):
        V = pts[[v.nr for v in el.vertices]]
        x, y = V[:, 0], V[:, 1]
        if np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y) < 0:
            V = V[::-1]
        c = V.mean(0); eds = []
        for i in range(len(V)):
            P0 = V[i]; P1 = V[(i + 1) % len(V)]; t = P1 - P0; Lg = np.hypot(*t); th = t / Lg
            nout = np.array([th[1], -th[0]]); mid = 0.5 * (P0 + P1)
            if nout @ (mid - c) < 0:
                nout = -nout
            eds.append((P0, P1, mid, nout, Lg, c))
        out.append(eds)
    return out


def _seg_field(P0, P1, L, P):
    """Exact 2D field at P from a unit uniform line charge on [P0,P1] (atan, not atan2, per the
    interior-eval-point convention)."""
    t = (P1 - P0) / L; nrm = np.array([t[1], -t[0]])
    d = P - P0; u = d @ t; v = d @ nrm
    r1 = np.hypot(*(P - P0)); r2 = np.hypot(*(P - P1))
    Et = np.log(max(r1, 1e-300) / max(r2, 1e-300)) / (2 * np.pi)
    En = (np.arctan((L - u) / v) + np.arctan(u / v)) / (2 * np.pi) if abs(v) > 1e-13 else 0.0
    return Et * t + En * nrm


def _solve_evalpoint(mesh, chi, h0):
    """MSC EIEM2 eval-point closure: sigma_f = chi H_n at 0.5(mid_f + centroid).  Returns M_avg (2,)."""
    elems = _elem_edges(mesh)
    flat = [(ei, k) for ei, eds in enumerate(elems) for k in range(len(eds))]
    idx = {fk: g for g, fk in enumerate(flat)}
    E = len(flat)
    evpt = np.array([0.5 * (elems[ei][k][2] + elems[ei][k][5]) for ei, k in flat])
    nout = np.array([elems[ei][k][3] for ei, k in flat])
    A = np.eye(E); rhs = np.empty(E); h0 = np.asarray(h0, float)
    for f, (ei, k) in enumerate(flat):
        nf = nout[f]; rhs[f] = chi * (h0 @ nf)
        for g, (ej, l) in enumerate(flat):
            P0, P1, _, _, Lg, _ = elems[ej][l]
            A[f, g] -= chi * (_seg_field(P0, P1, Lg, evpt[f]) @ nf)
    sig = np.linalg.solve(A, rhs)
    num = np.zeros(2); atot = 0.0
    for ei, eds in enumerate(elems):
        Vv = np.array([e[0] for e in eds]); x, y = Vv[:, 0], Vv[:, 1]
        atot += 0.5 * abs(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
        for k, (_, _, mid, _, Lg, _) in enumerate(eds):
            num += sig[idx[(ei, k)]] * Lg * mid
    return num / atot


def test_moment_closure_beats_evalpoint_on_shear():
    shears = [0.0, 0.4, 0.8, 1.2]
    rows = []
    with ng.TaskManager():
        base = _base_disk(0.22, quad=True)
        # sanity: eval-point on the regular disk recovers the disk demag to a few %
        Ma0 = _analytic_M(0.0, CHI, H0)
        Me0 = _solve_evalpoint(base, CHI, H0)
        assert abs(Me0[0] - Ma0[0]) / abs(Ma0[0]) < 0.03, (Me0, Ma0)   # eval-point is SANE (~1%)
        for s in shears:
            mesh = _sheared(base, s) if s > 0 else base
            Ma = _analytic_M(s, CHI, H0)
            Mm = np.asarray(m2.solve_planar_demag(mesh, mu_r=1 + CHI, H_ext=H0)["M_avg"], float)
            Me = _solve_evalpoint(mesh, CHI, H0)
            em = float(np.linalg.norm(Mm - Ma) / np.linalg.norm(Ma))
            ee = float(np.linalg.norm(Me - Ma) / np.linalg.norm(Ma))
            rows.append({"shear": s, "moment_relerr": em, "evalpoint_relerr": ee, "ratio": ee / em})
    # THE CLAIM: moment closure is far more accurate than eval-point at every shear, and stays small.
    for r in rows:
        assert r["moment_relerr"] < r["evalpoint_relerr"], r
        assert r["moment_relerr"] < 2e-3, r                            # moment stays < 0.2% to s=1.2
        assert r["ratio"] > 10.0, r                                    # >= an order of magnitude better
    summary = {
        "claim": "multipole-moment closure vs MSC eval-point (EIEM2), same surface-charge DOF, "
                 "sheared quad disk (analytic ellipse reference)",
        "chi": CHI, "H0": list(H0), "radia_version": radia.__version__,
        "sanity_regular_evalpoint_relerr": float(abs(Me0[0] - Ma0[0]) / abs(Ma0[0])),
        "rows": rows,
    }
    with open(os.path.join(HERE, "test_moment_vs_evalpoint_distortion_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
