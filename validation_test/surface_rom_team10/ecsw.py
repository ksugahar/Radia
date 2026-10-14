"""ECSW (energy-conserving sampling and weighting) for the steel magnetic term.

Training rows: for FOM training snapshots u_s and basis columns phi_j,
  C[(s,j), e] = int_e H(B_s) . curl phi_j dx      (element e in steel)
  d[(s,j)]    = sum_e C[(s,j), e]                 (the exact reduced force)
Each snapshot block is normalised by ||d_s||.  Non-negative least squares
(Lawson-Hanson, early exit at ||C xi - d|| <= tol ||d||) picks few elements
with weights xi >= 0.
"""
import time
import numpy as np
from ngsolve import *
from team13_model import *


def steel_elements(mesh):
    return np.array([el.nr for el in mesh.Elements(VOL) if el.mat == 'steel'])


def _bfield(m):
    """B of a GridFunction: curl(A) (A* model) or curl of the A component (A-phi model)."""
    return m.bfield if hasattr(m, 'bfield') else (lambda g: curl(g))


def contributions(m, Phi, snaps):
    """C (rows = snapshots x columns, cols = steel elements) and d."""
    _, bs = energy_density(); mesh = m.mesh; steel = mesh.Materials('steel'); sel = steel_elements(mesh)
    bf = _bfield(m); cols = []
    for j in range(Phi.shape[1]):
        g = GridFunction(m.fes); g.vec.FV().NumPy()[:] = Phi[:, j]; cols.append(bf(g))
    gu = GridFunction(m.fes); blocks = []
    for u in snaps:
        gu.vec.FV().NumPy()[:] = u; B = bf(gu); nb = sqrt(1e-24 + B * B); H = bs(nb) / nb * B
        rows = []
        for cj in cols:
            vals = Integrate(H * cj, mesh, definedon=steel, element_wise=True, order=2 * m.order + 2)
            rows.append(np.asarray(vals)[sel])
        blk = np.array(rows); nrm = np.linalg.norm(blk.sum(axis=1)); blocks.append(blk / max(nrm, 1e-300))
    C = np.concatenate(blocks, axis=0); return C, C.sum(axis=1), sel


class NumpyHR:
    """Online reduced nonlinear magnetic force / tangent at sampled elements,
    pure numpy (no NGSolve assembly):  B_q = D_q q,  F = sum_q w_q D_q^T H(B_q),
    J = sum_q w_q D_q^T dH/dB(B_q) D_q,  energy = sum_q w_q W(|B_q|).
    H(B) is the piecewise-linear interpolant of the B-H table (= BSpline order 2)."""
    def __init__(self, m, Phi, elements, weights, intorder):
        mesh = m.mesh; steel = mesh.Materials('steel'); sel = steel_elements(mesh)
        pos = {e: i for i, e in enumerate(sel)}; ir = IntegrationRule(ET.TET, intorder); nq = len(ir.weights)
        pts = mesh.MapToAllElements(ir, steel)
        rows = np.concatenate([np.arange(pos[e] * nq, (pos[e] + 1) * nq) for e in elements])
        sub = pts[rows]; self.pts = sub; self.nq_el = nq
        J = specialcf.JacobianMatrix(3)(sub).reshape(-1, 3, 3); det = np.abs(np.linalg.det(J))
        wref = np.tile(np.asarray(ir.weights), len(elements)); xi = np.repeat(np.asarray(weights), nq)
        self.w = wref * det * xi; self.n_elements = len(elements); self.n_points = len(sub)
        g = GridFunction(m.fes); D = np.empty((len(sub), 3, Phi.shape[1])); bf = _bfield(m)
        for j in range(Phi.shape[1]):
            g.vec.FV().NumPy()[:] = Phi[:, j]; D[:, :, j] = np.asarray(bf(g)(sub)).reshape(-1, 3)
        self.D = D
        self.Df = np.ascontiguousarray(D.reshape(-1, D.shape[2]))   # (points*3, rank): BLAS-friendly layout
        Bt, Ht = bh_table(); self.Bt = np.r_[0.0, Bt[Bt > 0]] if Bt[0] > 0 else Bt; self.Ht = np.interp(self.Bt, Bt, Ht)
        self.slope = np.diff(self.Ht) / np.diff(self.Bt)
        self.Wt = np.r_[0.0, np.cumsum(0.5 * (self.Ht[1:] + self.Ht[:-1]) * np.diff(self.Bt))]

    def _h(self, b):
        k = np.clip(np.searchsorted(self.Bt, b) - 1, 0, len(self.slope) - 1)
        return self.Ht[k] + self.slope[k] * (b - self.Bt[k]), self.slope[k], k

    def force_tangent(self, q):
        """F = sum_q w_q D_q^T H(B_q),  J = sum_q w_q D_q^T dH/dB D_q with
        dH/dB = (h/b) I + (h' - h/b) u u^T.  All contractions are BLAS matrix products."""
        Df = self.Df; nq = len(self.w)
        B = (Df @ q).reshape(nq, 3); b = np.sqrt(np.sum(B * B, axis=1) + 1e-24); h, hp, _ = self._h(b)
        a1 = h / b; a2 = hp - a1; u = B / b[:, None]
        F = Df.T @ ((self.w * a1)[:, None] * B).ravel()
        Jm = Df.T @ (Df * np.repeat(self.w * a1, 3)[:, None])        # isotropic part
        DB = (self.D * u[:, :, None]).sum(axis=1)                      # u^T D_q per point, (points, rank)
        Jm += DB.T @ (DB * (self.w * a2)[:, None])                     # rank-one part
        return F, 0.5 * (Jm + Jm.T)

    def force_tangent_einsum(self, q):
        """Reference implementation (the previous einsum version), kept for the equivalence check."""
        B = self.D @ q; b = np.sqrt(np.sum(B * B, axis=1) + 1e-24); h, hp, _ = self._h(b)
        F = np.einsum('q,qin,qi->n', self.w, self.D, (h / b)[:, None] * B)
        u = B / b[:, None]; a1 = h / b; a2 = hp - a1
        DB = np.einsum('qin,qi->qn', self.D, u)
        Jm = np.einsum('q,qin,qim->nm', self.w * a1, self.D, self.D) + np.einsum('q,qn,qm->nm', self.w * a2, DB, DB)
        return F, 0.5 * (Jm + Jm.T)

    def energy(self, q):
        B = self.D @ q; b = np.sqrt(np.sum(B * B, axis=1)); h, hp, k = self._h(b)
        W = self.Wt[k] + 0.5 * (self.Ht[k] + h) * (b - self.Bt[k])
        return float(self.w @ W)


def contributions_fast(m, Phi, snaps, intorder=None, chunk=None):
    """Same C and d as `contributions`, built from point data instead of element-wise
    NGSolve integrals: the basis curls are evaluated once at every steel quadrature
    point (NumpyHR with all elements, unit weights), each snapshot's B once, and the
    element sums are numpy reductions.  C[(s,j), e] = sum_{p in e} w_p D_pj . H(B_s,p)."""
    intorder = 2 * m.order + 2 if intorder is None else intorder
    sel = steel_elements(m.mesh)
    hr = NumpyHR(m, Phi, sel, np.ones(len(sel)), intorder)
    bf = _bfield(m); g = GridFunction(m.fes); ne, nq, nr = len(sel), hr.nq_el, Phi.shape[1]
    blocks = []
    for u in snaps:
        g.vec.FV().NumPy()[:] = u
        Bs = np.asarray(bf(g)(hr.pts)).reshape(-1, 3)
        b = np.sqrt(np.sum(Bs * Bs, axis=1) + 1e-24); h, _, _ = hr._h(b)
        X = (hr.w * h / b)[:, None] * Bs                                  # w_p H(B_p)
        P = np.einsum('pin,pi->pn', hr.D, X, optimize=True)              # per-point contribution to each column
        blk = P.reshape(ne, nq, nr).sum(axis=1).T                        # (columns, elements)
        nrm = np.linalg.norm(blk.sum(axis=1)); blocks.append(blk / max(nrm, 1e-300))
    C = np.concatenate(blocks, axis=0)
    return C, C.sum(axis=1), sel


def nnls_early(C, d, tol, max_elems=5000, log=print):
    n = C.shape[1]; P = np.zeros(n, bool); xi = np.zeros(n); r = d.copy(); nd = np.linalg.norm(d); t0 = time.time()
    while np.linalg.norm(r) > tol * nd and P.sum() < max_elems:
        g = C.T @ r; g[P] = -np.inf; j = int(np.argmax(g))
        if g[j] <= 0: break
        P[j] = True
        while True:
            idx = np.where(P)[0]; z = np.linalg.lstsq(C[:, idx], d, rcond=None)[0]
            if np.all(z > 0): xi[:] = 0; xi[idx] = z; break
            neg = z <= 0; a = np.min(xi[idx][neg] / np.maximum(xi[idx][neg] - z[neg], 1e-300))
            xi[idx] = xi[idx] + a * (z - xi[idx]); P[idx[xi[idx] <= 1e-14]] = False; xi[~P] = 0
        r = d - C @ xi
        if P.sum() % 100 == 0: log(f'  nnls {P.sum()} elems, rel res {np.linalg.norm(r)/nd:.2e}, {time.time()-t0:.0f}s')
    idx = np.where(P)[0]
    return idx, xi[idx], float(np.linalg.norm(d - C @ xi) / nd)
