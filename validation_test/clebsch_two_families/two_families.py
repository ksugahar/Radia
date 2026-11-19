"""Two magnetic-surface families of a nonaxisymmetric analytic field.

Prototype (diagnostic only, no acceptance criterion). The analytic map
    x = u + a sin(v) sin(s),  y = v + b sin(x) sin(2 s),  z = s,  a = b = 0.5
has det(DX) = 1, so B = dX/ds = grad(u) x grad(v) with B_z = 1. The labels
(u, v) are the oracle only: the reconstruction sees the B evaluator and the
inlet seed section z = 0, where (u, v) = (x, y).

Three error sources are measured separately:
  (1) projection: B_h = curl(A_h), A_h in HCurl(order k) on a structured
      hexahedral box (simply connected, natural boundary condition), fitted by
      least squares to B with a small gauge mass term; B_h is copied into
      HDiv(order k) to measure cell divergence and one-sided face normal flux;
  (2) tracking: a fixed-step 3-stage Gauss-Legendre collocation tracker that
      stops at every element face (same physical crossing point, one-sided
      normal B_h checked), compared with an independent adaptive DOP853 trace
      of the same B_h and with transported patch flux (Stokes, loop A_h.dl);
  (3) reconstruction: surfaces of both families are rebuilt from tracked
      curves (linear across curves, cubic Hermite along them); signed and
      absolute normal flux are integrated with Gauss points between the
      tracked curves (on the curves the integrand vanishes by construction).
      Oracle-sampled surfaces with the exact B isolate pure reconstruction.

Usage: python two_families.py --out result.json [--configs 4:2,8:3,2x1x1:6] [--threads 6]
"""
import argparse
import hashlib
import json
import math
import os
import platform
import time

import numpy as np
from scipy.integrate import solve_ivp

import ngsolve as ngs
from ngsolve.meshes import MakeStructured3DMesh

A_COEF = 0.5
B_COEF = 0.5
HALF = 2.0
ZEND = math.pi
SURF_LABELS = (-0.6, 0.0, 0.6)


# ----------------------------------------------------------------- analytic
def exact_b(x, y, z):
    v = y - B_COEF * np.sin(x) * np.sin(2 * z)
    bx = A_COEF * np.sin(v) * np.cos(z)
    by = B_COEF * np.cos(x) * bx * np.sin(2 * z) + 2 * B_COEF * np.sin(x) * np.cos(2 * z)
    return np.stack([bx, by, np.ones_like(x)], -1)


def labels(x, y, z):
    v = y - B_COEF * np.sin(x) * np.sin(2 * z)
    u = x - A_COEF * np.sin(v) * np.sin(z)
    return u, v


def label_grad(x, y, z, fam, d=1e-6):
    g = []
    for dx, dy, dz in ((d, 0, 0), (0, d, 0), (0, 0, d)):
        p = labels(x + dx, y + dy, z + dz)[fam]
        m = labels(x - dx, y - dy, z - dz)[fam]
        g.append((p - m) / (2 * d))
    return np.stack(g, -1)


def forward_map(u, v, s):
    x = u + A_COEF * np.sin(v) * np.sin(s)
    y = v + B_COEF * np.sin(x) * np.sin(2 * s)
    return x, y


def exact_cf():
    x, y, z = ngs.x, ngs.y, ngs.z
    v = y - B_COEF * ngs.sin(x) * ngs.sin(2 * z)
    bx = A_COEF * ngs.sin(v) * ngs.cos(z)
    by = B_COEF * ngs.cos(x) * bx * ngs.sin(2 * z) + 2 * B_COEF * ngs.sin(x) * ngs.cos(2 * z)
    return ngs.CF((bx, by, 1.0))


# --------------------------------------------------------------- projection
def build_mesh(nx, ny, nz):
    m = MakeStructured3DMesh(hexes=True, nx=nx, ny=ny, nz=nz,
                             mapping=lambda X, Y, Z: (-HALF + 2 * HALF * X, -HALF + 2 * HALF * Y, ZEND * Z))
    return m if isinstance(m, ngs.Mesh) else ngs.Mesh(m)


def project(mesh, order, eps):
    bex = exact_cf()
    fes = ngs.HCurl(mesh, order=order, nograds=True)
    u, w = fes.TnT()
    a = ngs.BilinearForm(fes, symmetric=True)
    a += (ngs.curl(u) * ngs.curl(w) + eps * u * w) * ngs.dx
    f = ngs.LinearForm(bex * ngs.curl(w) * ngs.dx(bonus_intorder=6))
    t0 = time.perf_counter()
    a.Assemble()
    f.Assemble()
    gfa = ngs.GridFunction(fes)
    inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    gfa.vec.data = inv * f.vec
    r = f.vec.CreateVector()
    r.data = f.vec - a.mat * gfa.vec
    rel_res = ngs.Norm(r) / ngs.Norm(f.vec)
    if rel_res > 1e-6:
        raise RuntimeError("projection true residual %.3e > 1e-6" % rel_res)
    t_solve = time.perf_counter() - t0
    bh = ngs.curl(gfa)
    fesd = ngs.HDiv(mesh, order=order)
    gfb = ngs.GridFunction(fesd)
    gfb.Set(bh, bonus_intorder=4)
    n = ngs.specialcf.normal(3)
    io = 2 * order + 6
    one = ngs.GridFunction(ngs.L2(mesh, order=0))
    one.Set(1.0)
    inner = one.Other()  # 1 on interior facets, 0 on the domain boundary
    jump = (bh - bh.Other()) * inner
    diag = {
        "ndof_hcurl": fes.ndof, "ndof_hdiv": fesd.ndof, "eps_gauge": eps, "t_assemble_solve_s": t_solve,
        "true_residual_rel": rel_res,
        "hdiv_copy_l2": math.sqrt(ngs.Integrate(ngs.InnerProduct(gfb - bh, gfb - bh), mesh, order=io)),
        "cell_div_l2": math.sqrt(ngs.Integrate(ngs.div(gfb) ** 2, mesh, order=io)),
        "face_normal_jump_l2": math.sqrt(ngs.Integrate(
            (jump * n) ** 2 * ngs.dx(element_boundary=True, bonus_intorder=io), mesh)),
        "face_tangent_jump_l2": math.sqrt(ngs.Integrate(
            ngs.InnerProduct(ngs.Cross(jump, n), ngs.Cross(jump, n))
            * ngs.dx(element_boundary=True, bonus_intorder=io), mesh)),
        "boundary_flux_mismatch_rel": ngs.Integrate(ngs.sqrt(((bh - bex) * n) ** 2) * ngs.ds(bonus_intorder=io), mesh)
        / ngs.Integrate(ngs.sqrt((bex * n) ** 2) * ngs.ds(bonus_intorder=io), mesh),
        "l2_err_rel": math.sqrt(ngs.Integrate(ngs.InnerProduct(bh - bex, bh - bex), mesh, order=io)
                                / ngs.Integrate(ngs.InnerProduct(bex, bex), mesh, order=io)),
    }
    return gfa, bh, diag


def heldout(mesh, bh, rng):
    npt = 4000
    p = np.column_stack([rng.uniform(-1.8, 1.8, npt), rng.uniform(-1.8, 1.8, npt), rng.uniform(0.05, ZEND - 0.05, npt)])
    bv = np.asarray(bh(mesh(p[:, 0], p[:, 1], p[:, 2])))
    be = exact_b(p[:, 0], p[:, 1], p[:, 2])
    e = np.linalg.norm(bv - be, axis=1)
    mb, me = np.linalg.norm(bv, axis=1), np.linalg.norm(be, axis=1)
    cosang = np.clip(np.sum(bv * be, 1) / (mb * me), -1, 1)
    return {"n_points": npt, "vector_rms_rel": float(np.sqrt(np.mean(e ** 2) / np.mean(me ** 2))),
            "vector_max_rel": float(np.max(e / me)), "magnitude_rms_rel": float(np.sqrt(np.mean(((mb - me) / me) ** 2))),
            "direction_max_rad": float(np.max(np.arccos(cosang))), "direction_rms_rad": float(np.sqrt(np.mean(np.arccos(cosang) ** 2)))}


def cut_flux(mesh, bh, fx, fy, zc):
    # flux of B_h,z through [-1, 1]^2 at z = zc, Gauss per interval split at element faces
    g, w = np.polynomial.legendre.leggauss(12)

    def rule(f):
        edges = np.unique(np.concatenate([[-1.0, 1.0], f[(f > -1) & (f < 1)]]))
        return (np.concatenate([0.5 * (a + b) + 0.5 * (b - a) * g for a, b in zip(edges[:-1], edges[1:])]),
                np.concatenate([0.5 * (b - a) * w for a, b in zip(edges[:-1], edges[1:])]))
    xs, wx = rule(fx)
    ys, wy = rule(fy)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    W = np.outer(wx, wy)
    bz = np.asarray(bh(mesh(X.ravel(), Y.ravel(), np.full(X.size, zc))))[:, 2]
    flux = float(np.sum(W.ravel() * bz))
    return {"z": zc, "flux_h": flux, "flux_exact": 4.0, "rel_err": abs(flux / 4.0 - 1)}


# ----------------------------------------------------------------- tracking
SQ15 = math.sqrt(15.0)
GC = np.array([0.5 - SQ15 / 10, 0.5, 0.5 + SQ15 / 10])
GA = np.array([[5 / 36, 2 / 9 - SQ15 / 15, 5 / 36 - SQ15 / 30],
               [5 / 36 + SQ15 / 24, 2 / 9, 5 / 36 - SQ15 / 24],
               [5 / 36 + SQ15 / 30, 2 / 9 + SQ15 / 15, 5 / 36]])
GB = np.array([5 / 18, 4 / 9, 5 / 18])


def _beta():
    out = []
    for j in range(3):
        p = np.poly1d([1.0])
        for m in range(3):
            if m != j:
                p = p * np.poly1d([1.0, -GC[m]]) / (GC[j] - GC[m])
        out.append(np.polyint(p))
    return out


BETA = _beta()


class Field:
    """B_h evaluator. __call__ uses NGSolve point location (it assigns points
    within ~1e-5 of a face to one side); local() evaluates the polynomial of
    an explicitly chosen structured cell, so tracking and one-sided face
    checks are cell-local."""

    def __init__(self, mesh, cf, fx, fy, fz):
        self.mesh, self.cf, self.calls, self.points = mesh, cf, 0, 0
        self.fx, self.fy, self.fz = fx, fy, fz
        self.nx, self.ny, self.nz = len(fx) - 1, len(fy) - 1, len(fz) - 1
        t = np.asarray(mesh(np.array([0.5 * (fx[0] + fx[1])]), np.array([0.5 * (fy[0] + fy[1])]),
                            np.array([0.5 * (fz[0] + fz[1])])))
        self.dtype, self.meshptr, self.vorb = t.dtype, t["meshptr"][0], t["VorB"][0]
        rng = np.random.default_rng(7)
        p = np.column_stack([rng.uniform(fx[0], fx[-1], 500), rng.uniform(fy[0], fy[-1], 500),
                             rng.uniform(fz[0], fz[-1], 500)])
        # keep the check away from faces: the locator snaps points within ~1e-5 of a face to either side
        far = np.ones(len(p), bool)
        for i, f in enumerate((fx, fy, fz)):
            far &= np.min(np.abs(p[:, i][:, None] - f[None, :]), axis=1) > 1e-3 * (f[1] - f[0])
        p = p[far]
        c = [np.clip(np.searchsorted(f, p[:, i], side="right") - 1, 0, len(f) - 2) for i, f in enumerate((fx, fy, fz))]
        mp = np.asarray(mesh(p[:, 0], p[:, 1], p[:, 2]))
        ref = np.asarray(cf(mp))
        loc = self.local(p[:, :2], p[:, 2], *c)
        if not (np.array_equal(mp["nr"], self.number(*c)) and np.max(np.abs(loc - ref)) < 1e-12 * max(1.0, np.max(np.abs(ref)))):
            raise RuntimeError("structured cell numbering / reference map does not match the mesh")
        self.calls = self.points = 0

    def number(self, cx, cy, cz):
        return (cx * self.ny + cy) * self.nz + cz

    def __call__(self, x, y, z):
        self.calls += 1
        self.points += len(x)
        return np.asarray(self.cf(self.mesh(x, y, z)))

    def local(self, X, z, cx, cy, cz):
        self.calls += 1
        self.points += len(X)
        arr = np.zeros(len(X), dtype=self.dtype)
        arr["meshptr"] = self.meshptr
        arr["VorB"] = self.vorb
        arr["facetnr"] = -1
        arr["eb"] = 0
        arr["nr"] = self.number(cx, cy, cz)
        arr["x"] = (z - self.fz[cz]) / (self.fz[cz + 1] - self.fz[cz])
        arr["z"] = (X[:, 0] - self.fx[cx]) / (self.fx[cx + 1] - self.fx[cx])
        arr["y"] = (X[:, 1] - self.fy[cy]) / (self.fy[cy + 1] - self.fy[cy])
        return np.asarray(self.cf(arr))

    def slope(self, X, z, cells=None):
        b = self(X[:, 0], X[:, 1], z) if cells is None else self.local(X, z, *cells)
        return b[:, :2] / b[:, 2:3]


def gl_step(field, X, z, h, cells):
    """3-stage Gauss-Legendre collocation step, all stages in the given cells."""
    reps = [np.repeat(c, 3) for c in cells]
    K = np.repeat(field.slope(X, z, cells)[:, None, :], 3, axis=1)
    d = 0.0
    for _ in range(80):
        ys = X[:, None, :] + h[:, None, None] * np.einsum("ij,njd->nid", GA, K)
        zs = z[:, None] + GC[None, :] * h[:, None]
        kn = field.slope(ys.reshape(-1, 2), zs.ravel(), reps).reshape(-1, 3, 2)
        d = np.max(np.abs(kn - K)) if kn.size else 0.0
        K = kn
        if d < 1e-15:
            break
    else:
        if d > 1e-11:
            raise RuntimeError("Gauss-Legendre stages did not converge (%.2e)" % d)
    return X + h[:, None] * np.einsum("j,njd->nd", GB, K), K


def dense(X, K, h, theta, comp):
    b = np.stack([np.polyval(BETA[j], theta) for j in range(3)], -1)
    return X[:, comp] + h * np.sum(b * K[:, :, comp], -1)


def _start_cells(field, X, z):
    """Cell of each seed; a seed on a face goes to the side B_h moves it into."""
    cells = []
    for c, f in ((0, field.fx), (1, field.fy)):
        idx = np.clip(np.searchsorted(f, X[:, c], side="right") - 1, 0, len(f) - 2)
        dist = np.abs(X[:, c][:, None] - f[None, :])
        on = np.where(np.min(dist, axis=1) < 1e-12)[0]
        if on.size:
            b = field(X[on, 0], X[on, 1], z[on])
            j = np.argmin(dist[on], axis=1)
            idx[on] = np.clip(np.where(b[:, c] >= 0, j, j - 1), 0, len(f) - 2)
        cells.append(idx)
    return cells


def track(field, X0, z_st, hmax, check_faces=True):
    """Advance all seeds through z stations (which must include every z face),
    stopping at every x/y element face. Stages are evaluated in the current
    cell only; a crossing lands on the face, keeps the same physical point and
    label pair, and continues in the neighbouring cell. One-sided values at the
    landing point are evaluated exactly in both cells."""
    fx, fy, fz = field.fx, field.fy, field.fz
    if not np.all(np.min(np.abs(fz[:, None] - z_st[None, :]), axis=1) < 1e-12):
        raise ValueError("z stations must include every z element face")
    N = len(X0)
    X = X0.astype(float).copy()
    z = np.full(N, z_st[0])
    cx, cy = _start_cells(field, X, z)
    pos, der = [X.copy()], [field.slope(X, z, (cx, cy, np.zeros(N, int)))]
    cross = {"n": 0, "normal_jump_max_rel": 0.0, "tangent_jump_max_rel": 0.0, "tangent_jump_rms_rel": 0.0,
             "landing_residual_max": 0.0, "grazing_min_normal_over_b": float("inf"),
             "z_face_normal_jump_max_rel": 0.0, "z_face_tangent_jump_max_rel": 0.0, "per_seed_counts": None}
    tsum, tcnt = 0.0, 0
    counts = np.zeros(N, int)
    tol = 1e-12
    ar = np.arange
    for k in range(1, len(z_st)):
        zt, z0k = z_st[k], z_st[k - 1]
        cz_k = int(np.clip(np.searchsorted(fz, 0.5 * (z0k + zt), side="right") - 1, 0, len(fz) - 2))
        hnom = (zt - z0k) / math.ceil((zt - z0k) / hmax - 1e-9)
        while True:
            act = np.where(z < zt - 1e-13)[0]
            if act.size == 0:
                break
            h = np.minimum(hnom, zt - z[act])
            Xa, za = X[act], z[act]
            czv = np.full(act.size, cz_k)
            Xn, K = gl_step(field, Xa, za, h, (cx[act], cy[act], czv))
            lox, hix, loy, hiy = fx[cx[act]], fx[cx[act] + 1], fy[cy[act]], fy[cy[act] + 1]
            ex = (Xn[:, 0] < lox - tol) | (Xn[:, 0] > hix + tol)
            ey = (Xn[:, 1] < loy - tol) | (Xn[:, 1] > hiy + tol)
            ev = np.where(ex | ey)[0]
            if ev.size:
                th = np.full(ev.size, 2.0)
                face = np.zeros(ev.size)
                comp = np.zeros(ev.size, int)
                for c, e_, lo, hi in ((0, ex, lox, hix), (1, ey, loy, hiy)):
                    sel = np.where(e_[ev])[0]
                    if not sel.size:
                        continue
                    i = ev[sel]
                    fv = np.where(Xn[i, c] < lo[i], lo[i], hi[i])
                    a_, b_ = np.zeros(sel.size), np.ones(sel.size)
                    s0 = -np.sign(Xn[i, c] - fv)
                    for _ in range(60):
                        mid = 0.5 * (a_ + b_)
                        g = dense(Xa[i], K[i], h[i], mid, c) - fv
                        inside = np.sign(g) == s0
                        a_ = np.where(inside, mid, a_)
                        b_ = np.where(inside, b_, mid)
                    t_c = 0.5 * (a_ + b_)
                    better = t_c < th[sel]
                    th[sel] = np.where(better, t_c, th[sel])
                    face[sel] = np.where(better, fv, face[sel])
                    comp[sel] = np.where(better, c, comp[sel])
                i = ev
                cells_i = (cx[act][i], cy[act][i], czv[i])
                outward = np.sign(Xn[i][ar(i.size), comp] - face).astype(int)
                hh = h[i] * th
                for _ in range(8):
                    Xe, Ke = gl_step(field, Xa[i], za[i], hh, cells_i)
                    r = Xe[ar(i.size), comp] - face
                    if np.max(np.abs(r)) < 1e-14:
                        break
                    sl = field.slope(Xe, za[i] + hh, cells_i)[ar(i.size), comp]
                    hh = np.clip(hh - r / sl, 0.0, h[i])
                cross["landing_residual_max"] = max(cross["landing_residual_max"], float(np.max(np.abs(r))))
                Xe[ar(i.size), comp] = face
                Xn[i] = Xe
                h[i] = hh
                ze = za[i] + hh
                ncx = cells_i[0] + np.where(comp == 0, outward, 0)
                ncy = cells_i[1] + np.where(comp == 1, outward, 0)
                if np.any((ncx < 0) | (ncx >= len(fx) - 1) | (ncy < 0) | (ncy >= len(fy) - 1)):
                    raise RuntimeError("field line left the domain")
                if check_faces:
                    bm = field.local(Xe, ze, *cells_i)
                    bp = field.local(Xe, ze, ncx, ncy, cells_i[2])
                    nb = np.linalg.norm(bp, axis=1)
                    jn = np.abs(bp[ar(i.size), comp] - bm[ar(i.size), comp]) / nb
                    jt = np.sqrt(np.maximum(np.linalg.norm(bp - bm, axis=1) ** 2 - (jn * nb) ** 2, 0.0)) / nb
                    cross["normal_jump_max_rel"] = max(cross["normal_jump_max_rel"], float(np.max(jn)))
                    cross["tangent_jump_max_rel"] = max(cross["tangent_jump_max_rel"], float(np.max(jt)))
                    tsum += float(np.sum(jt ** 2))
                    tcnt += i.size
                    bn = 0.5 * np.abs(bp[ar(i.size), comp] + bm[ar(i.size), comp]) / nb
                    cross["grazing_min_normal_over_b"] = min(cross["grazing_min_normal_over_b"], float(np.min(bn)))
                g_act = act[i]
                cx[g_act], cy[g_act] = ncx, ncy
                counts[g_act] += 1
                cross["n"] += int(i.size)
            X[act] = Xn
            z[act] = za + h
        czv_all = np.full(N, cz_k)
        if check_faces and np.any(np.abs(fz - zt) < 1e-12) and k < len(z_st) - 1:
            bm = field.local(X, z, cx, cy, czv_all)
            bp = field.local(X, z, cx, cy, czv_all + 1)
            nb = np.linalg.norm(bp, axis=1)
            jn = np.abs(bp[:, 2] - bm[:, 2]) / nb
            jt = np.linalg.norm(bp[:, :2] - bm[:, :2], axis=1) / nb
            cross["z_face_normal_jump_max_rel"] = max(cross["z_face_normal_jump_max_rel"], float(np.max(jn)))
            cross["z_face_tangent_jump_max_rel"] = max(cross["z_face_tangent_jump_max_rel"], float(np.max(jt)))
        pos.append(X.copy())
        der.append(field.slope(X, z, (cx, cy, czv_all)))
    cross["tangent_jump_rms_rel"] = math.sqrt(tsum / tcnt) if tcnt else 0.0
    cross["per_seed_counts"] = {"min": int(counts.min()), "max": int(counts.max()), "mean": float(counts.mean())}
    return np.array(pos), np.array(der), cross


def track_reference(field, X0, z_st, rtol):
    def rhs(zz, s):
        p = s.reshape(-1, 2)
        zz_ = np.full(len(p), zz)
        cells = [np.clip(np.searchsorted(f, v, side="right") - 1, 0, len(f) - 2)
                 for f, v in ((field.fx, p[:, 0]), (field.fy, p[:, 1]), (field.fz, zz_))]
        return field.slope(p, zz_, cells).ravel()
    sol = solve_ivp(rhs, (z_st[0], z_st[-1]), X0.ravel(), method="DOP853", rtol=rtol, atol=rtol * 0.1, t_eval=z_st)
    if not sol.success:
        raise RuntimeError("reference trace failed: " + sol.message)
    return sol.y.T.reshape(len(z_st), -1, 2), int(sol.nfev)


# ------------------------------------------------------------ patch flux
def loop_points(x0, y0, x1, y1, m):
    t = np.linspace(0, 1, m + 1)
    e = [np.column_stack([x0 + (x1 - x0) * t, np.full(m + 1, y0)]),
         np.column_stack([np.full(m + 1, x1), y0 + (y1 - y0) * t]),
         np.column_stack([x1 + (x0 - x1) * t, np.full(m + 1, y1)]),
         np.column_stack([np.full(m + 1, x0), y1 + (y0 - y1) * t])]
    return e


def loop_flux(afield, edges, zc):
    """Signed flux through the closed polygon through the transported points:
    Stokes, closed line integral of A_h along straight chords, each chord split
    at every element face and integrated cell-locally (A_h.dl jumps where a
    chord crosses a face obliquely). Exact for the polygon up to quadrature;
    the polygon approximates the transported curve to O(spacing^2) away from
    kinks of the piecewise-smooth flow map."""
    poly = np.vstack([e[:-1] for e in edges])
    q0, q1 = poly, np.roll(poly, -1, axis=0)
    g, w = np.polynomial.legendre.leggauss(6)
    g, w = 0.5 * (g + 1), 0.5 * w
    P, D, W = [], [], []
    for a_, b_ in zip(q0, q1):
        cuts = [0.0, 1.0]
        for c, f in ((0, afield.fx), (1, afield.fy)):
            if b_[c] != a_[c]:
                tt = (f - a_[c]) / (b_[c] - a_[c])
                cuts += list(tt[(tt > 0) & (tt < 1)])
        cuts = np.unique(cuts)
        for t0, t1 in zip(cuts[:-1], cuts[1:]):
            tq = t0 + (t1 - t0) * g
            P.append(a_[None, :] + tq[:, None] * (b_ - a_)[None, :])
            D.append(np.repeat((b_ - a_)[None, :], len(g), 0))
            W.append((t1 - t0) * w)
    P, D, W = np.vstack(P), np.vstack(D), np.concatenate(W)
    mid = P  # Gauss points are interior to each piece, so floor gives its cell
    zz = np.full(len(P), zc)
    cells = [np.clip(np.searchsorted(f, v, side="right") - 1, 0, len(f) - 2)
             for f, v in ((afield.fx, mid[:, 0]), (afield.fy, mid[:, 1]), (afield.fz, zz))]
    av = afield.local(P, zz, *cells)
    return float(np.sum(W * (av[:, 0] * D[:, 0] + av[:, 1] * D[:, 1])))


# ---------------------------------------------------------- reconstruction
def hermite(p0, p1, d0, d1, dz, tau):
    h00, h10, h01, h11 = 2 * tau**3 - 3 * tau**2 + 1, tau**3 - 2 * tau**2 + tau, -2 * tau**3 + 3 * tau**2, tau**3 - tau**2
    g00, g10, g01, g11 = 6 * tau**2 - 6 * tau, 3 * tau**2 - 4 * tau + 1, -6 * tau**2 + 6 * tau, 3 * tau**2 - 2 * tau
    val = h00 * p0 + h10 * dz * d0 + h01 * p1 + h11 * dz * d1
    der = (g00 * p0 + g10 * dz * d0 + g01 * p1 + g11 * dz * d1) / dz
    return val, der


def surface_flux(pos, der, z_st, t_seed, bfun, fam, label0, q):
    """pos/der: (n_st, Ns, 2). Linear across curves, cubic Hermite along them."""
    g, w = np.polynomial.legendre.leggauss(q)
    g = 0.5 * (g + 1)
    w = 0.5 * w
    S, W, ST, SZ = [], [], [], []
    for k in range(len(z_st) - 1):
        dz = z_st[k + 1] - z_st[k]
        for tau, wt in zip(g, w):
            hv, hd = hermite(pos[k][:, None, :], pos[k + 1][:, None, :], der[k][:, None, :], der[k + 1][:, None, :], dz, tau)
            hv, hd = hv[:, 0, :], hd[:, 0, :]
            for sig, ws in zip(g, w):
                dt = np.diff(t_seed)
                s_xy = (1 - sig) * hv[:-1] + sig * hv[1:]
                st = (hv[1:] - hv[:-1]) / dt[:, None]
                sz = (1 - sig) * hd[:-1] + sig * hd[1:]
                S.append(np.column_stack([s_xy, np.full(len(dt), z_st[k] + tau * dz)]))
                ST.append(st)
                SZ.append(sz)
                W.append(wt * ws * dz * dt)
    S, ST, SZ, W = np.vstack(S), np.vstack(ST), np.vstack(SZ), np.concatenate(W)
    st3 = np.column_stack([ST, np.zeros(len(ST))])
    sz3 = np.column_stack([SZ, np.ones(len(SZ))])
    nds = np.cross(st3, sz3)
    bv = bfun(S)
    bn = np.sum(bv * nds, 1)
    norm = float(np.sum(W * np.linalg.norm(bv, axis=1) * np.linalg.norm(nds, axis=1)))
    lab = labels(S[:, 0], S[:, 1], S[:, 2])[fam]
    gr = np.linalg.norm(label_grad(S[:, 0], S[:, 1], S[:, 2], fam), axis=1)
    dist = np.abs(lab - label0) / gr
    return {"signed_over_norm": float(np.sum(W * bn)) / norm, "abs_over_norm": float(np.sum(W * np.abs(bn))) / norm,
            "abs": float(np.sum(W * np.abs(bn))), "signed": float(np.sum(W * bn)), "norm": norm,
            "dist_max": float(np.max(dist)), "dist_rms": float(np.sqrt(np.sum(W * dist ** 2 * np.linalg.norm(nds, axis=1))
                                                                        / np.sum(W * np.linalg.norm(nds, axis=1))))}


def seeds_family(fam, label0, ns):
    t = np.linspace(-1, 1, ns)
    if fam == 0:
        return t, np.column_stack([np.full(ns, label0), t])
    return t, np.column_stack([t, np.full(ns, label0)])


# -------------------------------------------------------------------- main
def run_config(nx, ny, nz, order, args, rng):
    mesh = build_mesh(nx, ny, nz)
    gfa, bh, pdiag = project(mesh, order, args.eps)
    fx = np.linspace(-HALF, HALF, nx + 1)
    fy = np.linspace(-HALF, HALF, ny + 1)
    fz = np.linspace(0, ZEND, nz + 1)
    out = {"n_xy": nx if nx == ny else None, "n_x": nx, "n_y": ny, "n_z": nz, "order": order,
           "h_x": 2 * HALF / nx, "h_y": 2 * HALF / ny, "h_z": ZEND / nz, "projection": pdiag,
           "heldout": heldout(mesh, bh, rng),
           "cuts": [cut_flux(mesh, bh, fx, fy, zc) for zc in (0.37 * ZEND, 0.61 * ZEND, 0.93 * ZEND)]}
    field = Field(mesh, bh, fx, fy, fz)
    afield = Field(mesh, gfa, fx, fy, fz)
    bfun_h = lambda S: field(S[:, 0], S[:, 1], S[:, 2])

    # (2) tracking vs independent reference on identical seeds
    t0 = time.perf_counter()
    rs = np.column_stack([rng.uniform(-1, 1, 48), rng.uniform(-1, 1, 48)])
    z_st = np.linspace(0, ZEND, nz * 2 + 1)
    trk = {}
    for hmax in args.hmax_list:
        p, _, cr = track(field, rs, z_st, hmax)
        trk[str(hmax)] = {"pos": p, "cross": cr}
    pref, nfev = track_reference(field, rs, z_st, args.ref_rtol)
    pref2, nfev2 = track_reference(field, rs, z_st, args.ref_rtol * 100)
    lab_u, lab_v = labels(pref[-1][:, 0], pref[-1][:, 1], np.full(len(rs), ZEND))
    out["tracking"] = {
        "n_seeds": len(rs), "ref_method": "DOP853", "ref_rtol": args.ref_rtol, "ref_nfev": nfev,
        "ref_self_check_vs_rtol_x100_max": float(np.max(np.abs(pref - pref2))),
        "label_drift_of_reference_at_outlet_max": float(max(np.max(np.abs(lab_u - rs[:, 0])), np.max(np.abs(lab_v - rs[:, 1])))),
        "gl": {k: {"point_err_max_vs_ref": float(np.max(np.linalg.norm(v["pos"] - pref, axis=2))), **{c: v["cross"][c] for c in v["cross"]}}
               for k, v in trk.items()},
        "t_s": time.perf_counter() - t0}

    # (2b) transported patch flux: Stokes loop integral of A_h, inlet vs outlet vs exact
    t0 = time.perf_counter()
    pf = {}
    po = args.seed_offset
    patches = [(x0 + po, y0 + po) for x0 in np.arange(-1, 1, 0.5) for y0 in np.arange(-1, 1, 0.5)]
    for m in args.loop_m:
        edges_all, counts = [], []
        for x0, y0 in patches:
            es = loop_points(x0, y0, x0 + 0.5, y0 + 0.5, m)
            edges_all.append(es)
        pts = np.vstack([e for es in edges_all for e in es])
        zp = np.unique(np.concatenate([fz, [0.5 * ZEND]]))
        P, _, cr = track(field, pts, zp, args.hmax_surface, check_faces=False)
        rows = []
        off = 0
        for es in edges_all:
            sizes = [len(e) for e in es]
            r = {}
            for si, zc in ((0, 0.0), (int(np.argmin(np.abs(zp - 0.5 * ZEND))), 0.5 * ZEND), (len(zp) - 1, ZEND)):
                ee, o = [], off
                for s_ in sizes:
                    ee.append(P[si][o:o + s_])
                    o += s_
                r[zc] = loop_flux(afield, ee, zc)
            off += sum(sizes)
            rows.append(r)
        fin = np.array([r[0.0] for r in rows])
        fmid = np.array([r[0.5 * ZEND] for r in rows])
        fout = np.array([r[ZEND] for r in rows])
        pf[str(m)] = {"inlet_vs_exact_max_rel": float(np.max(np.abs(fin / 0.25 - 1))),
                      "mid_vs_inlet_max_rel": float(np.max(np.abs(fmid / fin - 1))),
                      "outlet_vs_inlet_max_rel": float(np.max(np.abs(fout / fin - 1))),
                      "outlet_vs_exact_max_rel": float(np.max(np.abs(fout / 0.25 - 1))),
                      "sum_inlet": float(fin.sum()), "sum_outlet": float(fout.sum()), "crossings": cr["n"]}
    out["patch_flux"] = {"n_patches": len(patches), "patch_area_exact": 0.25, "by_loop_points_per_edge": pf,
                         "t_s": time.perf_counter() - t0}

    # (3) two families: oracle reconstruction (exact B) and tracked B_h reconstruction
    t0 = time.perf_counter()
    fams = []
    for fam in (0, 1):
        for label0 in [l + args.seed_offset for l in SURF_LABELS]:
            for ns in args.ns_list:
                tseed, X0 = seeds_family(fam, label0, ns)
                tseed = tseed + args.seed_offset
                X0[:, 1 - fam] += args.seed_offset
                for mz in args.mz_list:
                    z_st_s = np.linspace(0, ZEND, nz * mz + 1)
                    # oracle-sampled surface, exact B: pure reconstruction error
                    uu, vv = (np.full(ns, label0), tseed) if fam == 0 else (tseed, np.full(ns, label0))
                    po = np.array([np.column_stack(forward_map(uu, vv, zk)) for zk in z_st_s])
                    do = np.array([exact_b(p_[:, 0], p_[:, 1], np.full(ns, zk))[:, :2] for p_, zk in zip(po, z_st_s)])
                    # tracked surface, B_h
                    ph, dh, cr = track(field, X0, z_st_s, args.hmax_surface)
                    for q in args.q_list:
                        rec = {"family": "u" if fam == 0 else "v", "label": label0, "n_seed": ns, "z_per_cell": mz, "q": q,
                               "oracle_exactB": surface_flux(po, do, z_st_s, tseed, lambda S: exact_b(S[:, 0], S[:, 1], S[:, 2]), fam, label0, q),
                               "oracle_Bh": surface_flux(po, do, z_st_s, tseed, bfun_h, fam, label0, q),
                               "tracked_Bh": surface_flux(ph, dh, z_st_s, tseed, bfun_h, fam, label0, q),
                               "tracked_crossings": cr["n"], "tracked_normal_jump_max_rel": cr["normal_jump_max_rel"]}
                        fams.append(rec)
    out["surfaces"] = {"records": fams, "t_s": time.perf_counter() - t0}
    out["field_evaluations"] = {"calls": field.calls, "points": field.points}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--configs", default="4:2,4:3,8:2,8:3")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--eps", type=float, default=1e-9)
    ap.add_argument("--hmax-list", type=float, nargs="+", default=[0.2, 0.1, 0.05])
    ap.add_argument("--hmax-surface", type=float, default=0.05)
    ap.add_argument("--ref-rtol", type=float, default=1e-12)
    ap.add_argument("--loop-m", type=int, nargs="+", default=[4, 8, 16])
    ap.add_argument("--ns-list", type=int, nargs="+", default=[9, 17, 33])
    ap.add_argument("--mz-list", type=int, nargs="+", default=[2, 4])
    ap.add_argument("--q-list", type=int, nargs="+", default=[3, 6])
    ap.add_argument("--seed-offset", type=float, default=0.0,
                    help="shift patch partition and surface labels; 0 puts patch edges and u=0 seeds on element faces")
    args = ap.parse_args()
    ngs.SetNumThreads(args.threads)
    rng = np.random.default_rng(20261120)
    here = os.path.abspath(__file__)
    res = {"script": os.path.basename(here), "script_sha256": hashlib.sha256(open(here, "rb").read()).hexdigest(),
           "field": {"a": A_COEF, "b": B_COEF, "box_half": HALF, "z_end": ZEND, "surface_labels": SURF_LABELS},
           "args": vars(args), "versions": {"ngsolve": ngs.__version__, "numpy": np.__version__, "python": platform.python_version()},
           "threads": args.threads, "env_threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS")},
           "configs": []}
    t0 = time.perf_counter()
    with ngs.TaskManager():
        for c in args.configs.split(","):
            sz, k = c.split(":")
            dims = [int(d) for d in sz.split("x")]
            if len(dims) == 1:
                dims = dims * 3
            if len(dims) != 3:
                raise ValueError("config %r: use n:k or NXxNYxNZ:k" % c)
            k = int(k)
            tc = time.perf_counter()
            r = run_config(*dims, k, args, rng)
            r["t_total_s"] = time.perf_counter() - tc
            res["configs"].append(r)
            print("config %s k=%d done in %.1f s" % (sz, k, r["t_total_s"]), flush=True)
            with open(args.out, "w") as fh:
                json.dump(res, fh, indent=1, default=float)
    res["t_total_s"] = time.perf_counter() - t0
    with open(args.out, "w") as fh:
        json.dump(res, fh, indent=1, default=float)


if __name__ == "__main__":
    main()
