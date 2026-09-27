"""Surface-trace POD ROM on the A-phi model (compound vector [A | phi]).

Basis columns are full compound vectors:
  surface modes : POD of the gamma trace (steel-surface A DOFs), extended
                  - into the air by the exact linear air operator (k_rest), and
                  - into the steel interior (interior A and all phi DOFs) by the
                    reference operator nu_ref curl-curl + sigma/dt |u + grad p|^2;
  interior modes: POD of the zero-trace steel remainder (interior A and phi).
u = Phi q + i(t) e_src (e_src: unit-current air solution, zero in the steel) is exact in the air.
Exactly precomputed quadratic pieces (no approximation):
  air/Kelvin      Kr = Phi^T K_rest Phi,  c = Phi^T K_rest e_src - Phi^T f
  sigma term      1/(2dt)|Phi q - P_A Phi q_n|^2_Msig  -> Mr q - Mc q_n  (P_A keeps the A part of the previous step)
  steel eps/reg   S = Phi^T K_steel_lin Phi
Online nonlinear part: int_steel w(|curl A|), full steel assembly (hr=None) or ECSW (ecsw.NumpyHR).
"""
import time
import numpy as np
from ngsolve import *
from team13_model import *

NU_REF = NU_0 / 1000.0


def bits(n, idx):
    b = BitArray(n); b.Clear()
    for i in idx: b.Set(int(i))
    return b


class ExtenderAphi:
    def __init__(self, m):
        self.m = m; n = m.fes.ndof
        self.air_inv = m.k_rest.mat.Inverse(bits(n, m.idx_air), inverse='pardiso')
        (u, p), (v, q) = m.fes.TnT(); steel = m.mesh.Materials('steel'); s = SIGMA_PLACEHOLDER
        # (u, p) = (grad chi, -chi) has zero sigma and curl energy, so the lifting operator needs its own
        # regularisation.  Any SPD lifting is admissible: it only chooses how a trace is continued into the
        # steel; the air part of each column stays the exact harmonic extension.
        self.k_ref = BilinearForm(NU_REF * curl(u) * curl(v) * dx(steel) + s / m.dt * (u + grad(p)) * (v + grad(q)) * dx(steel)
                                  + 1e-6 * s / m.dt * (u * v + p * q) * dx(steel), symmetric=True).Assemble()
        self.int_inv = self.k_ref.mat.Inverse(bits(n, m.idx_int), inverse='pardiso')
        self.tmp = m.gfu.vec.CreateVector(); self.rhs = m.gfu.vec.CreateVector(); self.sol = m.gfu.vec.CreateVector()

    def _solve(self, mat, inv, g_full):
        self.tmp.FV().NumPy()[:] = g_full; mat.Mult(self.tmp, self.rhs); self.rhs.data *= -1
        self.sol.data = inv * self.rhs; return self.sol.FV().NumPy().copy()

    def extend(self, g):
        m = self.m; e = np.zeros(m.fes.ndof); e[m.idx_gamma] = g
        ea = self._solve(m.k_rest.mat, self.air_inv, e); ei = self._solve(self.k_ref.mat, self.int_inv, e)
        e[m.idx_air] = ea[m.idx_air]; e[m.idx_int] = ei[m.idx_int]; return e

    def interior_lift(self, g):
        m = self.m; e = np.zeros(m.fes.ndof); e[m.idx_gamma] = g
        return self._solve(self.k_ref.mat, self.int_inv, e)[m.idx_int]

    def source(self):
        """Unit-current air solution of the first coil (1D), or of every coil (columns) if the model has several."""
        E = self.sources()
        return E[:, 0] if E.shape[1] == 1 else E

    def sources(self):
        m = self.m; cols = []
        for fs in m.fs:
            self.sol.data = self.air_inv * fs.vec; e = np.zeros(m.fes.ndof); e[m.idx_air] = self.sol.FV().NumPy()[m.idx_air]; cols.append(e)
        return np.array(cols).T


def _pod(blocks, rank):
    X = np.concatenate([b / np.linalg.norm(b) for b in blocks], axis=1)
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    return U[:, :rank], float(np.sum(s[:rank] ** 2) / np.sum(s ** 2))


def build_basis_zero_lift(m, ext, train, k, r):
    """Surface modes carry values only on the surface: trace POD, exact air extension,
    ZERO continuation into the steel (no assumed reference permeability).  The steel
    interior (interior A and all phi) is carried entirely by the interior modes, a POD of
    the full interior part of the snapshots."""
    Vg, rg = _pod([S[1:, m.idx_gamma].T for S in train], k)
    Wi, ri = _pod([S[1:, m.idx_int].T for S in train], r) if r else (np.zeros((len(m.idx_int), 0)), 1.0)
    cols = []
    for j in range(Vg.shape[1]):
        e = np.zeros(m.fes.ndof); e[m.idx_gamma] = Vg[:, j]
        ea = ext._solve(m.k_rest.mat, ext.air_inv, e); e[m.idx_air] = ea[m.idx_air]; cols.append(e)
    for j in range(Wi.shape[1]):
        e = np.zeros(m.fes.ndof); e[m.idx_int] = Wi[:, j]; cols.append(e)
    return np.array(cols).T, dict(gamma_retention=rg, interior_retention=ri)


def build_basis_global(m, train, rank):
    """Ablation: plain POD of the full snapshot vectors (air, surface and steel together);
    used with u = Phi q (no e_src), so the air part is a combination of snapshot air fields."""
    P, ret = _pod([S[1:].T for S in train], rank)
    return P, dict(retention=ret)


def _mult_cols(m, mat, X):
    Y = np.empty_like(X); x = m.gfu.vec.CreateVector(); y = x.CreateVector()
    for j in range(X.shape[1]):
        x.FV().NumPy()[:] = X[:, j]; mat.Mult(x, y); Y[:, j] = y.FV().NumPy()
    return Y


def _sym(A): return 0.5 * (A + A.T)


def split_orthonormalize(m, Phi, rel_tol=1e-10):
    """Split every column into its A part and its phi part and orthonormalise each block
    separately (in the energy metric).  Only the A part of the previous step enters the
    sigma term, so the reduced previous-step operator Phi^T M_sigma (P_A Phi) must act
    inside the basis; with mixed [A | phi] columns P_A Phi leaves the span and backward
    Euler becomes unstable (observed growth of about x2.5 per step).  With a block basis
    P_A Phi = the A block, and the reduced scheme inherits the full-space stability.
    The air part lives in the A block, so the exact-air structure is kept."""
    PA, ia = energy_orthonormalize(m, Phi * m.isA[:, None], rel_tol)
    PP, ip = energy_orthonormalize(m, Phi * (~m.isA)[:, None], rel_tol)
    return np.column_stack((PA, PP)), dict(kept_A=ia['kept'], kept_phi=ip['kept'], kept=ia['kept'] + ip['kept'],
                                           dropped=ia['dropped'] + ip['dropped'],
                                           energy_condition=max(ia['energy_condition'], ip['energy_condition']))


def energy_orthonormalize(m, Phi, rel_tol=1e-10):
    """Remove (near) zero-energy combinations of the basis columns.
    The A-phi pair (grad chi, -chi) carries no sigma or curl energy; snapshots contain it
    with an arbitrary, eps-controlled amplitude, so a plain POD basis can span directions
    with ~zero energy and the reduced Newton matrix becomes singular.  The columns are
    re-combined to be orthonormal in K = M_sigma/dt + nu_ref C_steel + K_rest and
    directions with eigenvalue below rel_tol * max are dropped.  Linear combinations keep
    the exact-air structure (each column's air part is the harmonic extension of its trace)."""
    KP = _mult_cols(m, m.m_sig.mat, Phi) / m.dt + NU_REF * _mult_cols(m, m.c_steel.mat, Phi) + _mult_cols(m, m.k_rest.mat, Phi)
    G = _sym(Phi.T @ KP); lam, V = np.linalg.eigh(G)
    keep = lam > rel_tol * lam.max()
    return Phi @ (V[:, keep] / np.sqrt(lam[keep])[None, :]), dict(kept=int(keep.sum()), dropped=int((~keep).sum()),
                                                                  energy_condition=float(lam.max() / lam[keep].min()))


class FastROMAphi:
    def __init__(self, m, Phi, e_src, hr=None):
        # e_src: (n,) for one coil or (n, ns) with one unit-current air solution per coil; i(t) is then a vector
        self.m = m; self.Phi = Phi; E = np.asarray(e_src); self.E = E[:, None] if E.ndim == 1 else E
        t0 = time.time(); w, _ = energy_density()
        self.nhr = hr if (hr is not None and hasattr(hr, 'force_tangent')) else None
        if self.nhr is None:
            u = m.fes.TrialFunction()[0]
            self.a = BilinearForm(m.fes, symmetric=True)
            self.a += Variation(w(sqrt(1e-12 + curl(u) * curl(u))) * dx(m.mesh.Materials('steel'), bonus_intorder=2))
            self.n_elements = int(sum(1 for el in m.mesh.Elements(VOL) if el.mat == 'steel'))
        else:
            self.n_elements = hr.n_elements
        self.sidx = m.idx_steelA; self.Ps = Phi[self.sidx]
        PA = Phi * m.isA[:, None]                                    # A part only (previous step)
        F = np.array([fs.vec.FV().NumPy() for fs in m.fs]).T             # (n, ns) unit-current load vectors
        if F.shape[1] != self.E.shape[1]: raise ValueError(('sources', F.shape[1], 'air solutions', self.E.shape[1]))
        KP = _mult_cols(m, m.k_rest.mat, Phi); KE = _mult_cols(m, m.k_rest.mat, self.E)
        # linear part 1/2 u^T K u - (F i)^T u with u = Phi q + E i:
        #   1/2 q^T Kr q + q^T C i + 1/2 i^T EE i - i^T FE i
        self.Kr = _sym(Phi.T @ KP); self.C = Phi.T @ KE - Phi.T @ F
        self.EE = _sym(self.E.T @ KE); self.FE = F.T @ self.E
        MP = _mult_cols(m, m.m_sig.mat, Phi); MPA = _mult_cols(m, m.m_sig.mat, PA)
        self.Mr = _sym(Phi.T @ MP) / m.dt; self.Mc = (Phi.T @ MPA) / m.dt; self.Mcc = _sym(PA.T @ MPA) / m.dt
        self.S = _sym(Phi.T @ _mult_cols(m, m.k_steel_lin.mat, Phi))
        self.setup_s = time.time() - t0
        self.g = GridFunction(m.fes); self.r = self.g.vec.CreateVector(); self.col = self.g.vec.CreateVector(); self.out = self.g.vec.CreateVector()

    def _set(self, q):
        v = self.g.vec.FV().NumPy(); v[:] = 0; v[self.sidx] = self.Ps @ q

    def _cur(self, i):
        return np.atleast_1d(np.asarray(i, dtype=float))

    def _lin(self, q, qn, i):
        return self.Mr @ q - self.Mc @ qn + self.S @ q + self.Kr @ q + self.C @ self._cur(i)

    def energy(self, q, qn, i):
        if self.nhr is not None: en = self.nhr.energy(q)
        else: self._set(q); en = self.a.Energy(self.g.vec)
        c = self._cur(i)
        return (en + 0.5 * q @ self.Mr @ q - q @ self.Mc @ qn + 0.5 * qn @ self.Mcc @ qn + 0.5 * q @ (self.S + self.Kr) @ q
                + q @ self.C @ c + 0.5 * c @ self.EE @ c - c @ self.FE @ c)

    def _newton_system(self, q, qn, i, T):
        H0 = self.Mr + self.S + self.Kr
        if self.nhr is not None:
            t = time.time(); F, Jn = self.nhr.force_tangent(q); T['assemble'] += time.time() - t
            rr = F + self._lin(q, qn, i); Jr = Jn + H0
        else:
            t = time.time(); self._set(q); self.a.AssembleLinearization(self.g.vec); self.a.Apply(self.g.vec, self.r); T['assemble'] += time.time() - t
            t = time.time(); nr = len(q)
            rr = self.Ps.T @ self.r.FV().NumPy()[self.sidx] + self._lin(q, qn, i)
            JP = np.empty((len(self.sidx), nr)); colv = self.col.FV().NumPy()
            for j in range(nr):
                colv[:] = 0; colv[self.sidx] = self.Ps[:, j]; self.a.mat.Mult(self.col, self.out); JP[:, j] = self.out.FV().NumPy()[self.sidx]
            Jr = _sym(self.Ps.T @ JP) + H0; T['project'] += time.time() - t
        t = time.time(); dq = np.linalg.solve(Jr, rr); T['solve'] += time.time() - t
        return dq, abs(dq @ rr)

    def run(self, fn, steps, tol=1e-10, linesearch=False):
        dt = self.m.dt; nr = self.Phi.shape[1]; q = np.zeros(nr); qs = [q.copy()]; its = []; rejects = 0
        T = dict(assemble=0., project=0., solve=0., linesearch=0., other=0.); t_start = time.time()
        for k in range(1, steps + 1):
            i = fn(k * dt); qn = qs[-1]; prev = None
            for it in range(40):
                dq, dec = self._newton_system(q, qn, i, T)
                if it == 0: d0 = max(dec, 1e-300)
                if dec <= tol * d0 or dec < 1e-24: break
                ls = linesearch
                if prev is not None and dec > prev[1] and not linesearch:
                    rejects += 1; q, dq = prev[0], prev[2]; ls = True
                t = time.time(); tau = 1.
                if ls:
                    E0 = self.energy(q, qn, i)
                    while self.energy(q - tau * dq, qn, i) > E0 and tau >= 1e-4: tau *= .5
                T['linesearch'] += time.time() - t
                prev = (q.copy(), dec, dq.copy()); q = q - tau * dq
            else: raise RuntimeError(('A-phi ROM Newton', k))
            its.append(it); qs.append(q.copy())
        total = time.time() - t_start; T['other'] = total - sum(T.values())
        U = np.array([self.Phi @ qq for qq in qs])
        U[1:] += np.array([self._cur(fn(k * dt)) for k in range(1, steps + 1)]) @ self.E.T
        return U, dict(seconds=total, setup_s=self.setup_s, breakdown=T, newton_mean=float(np.mean(its)),
                       newton_max=int(max(its)), rejected_full_steps=rejects, hr_elements=self.n_elements)
