"""Surface-trace + interior POD Galerkin ROM for the TEAM-13 transient FOM.

Global reduced basis (full DOF vectors):
  surface modes : trace POD V (on steel-surface DOFs Gamma), each extended
                  - into the air/coil/Kelvin region by the exact linear air
                    operator (harmonic extension, K_rest), and
                  - into the steel interior by the reference steel operator
                    nu_ref curl-curl + sigma/dt mass;
  interior modes: POD of the zero-trace steel-interior remainder.
Source: i(t) * e_src, with e_src the air solution for unit current and zero
trace (exact, since the air region is linear).  Hence u = Phi q + i(t) e_src
is exact in the air for any trace; approximation enters only through the
Gamma-trace and interior truncation.

Newton on Phi^T R(Phi q + i e_src) = 0 with the full nonlinear assembly (no
hyper-reduction).  Usage:
  python rom.py <maxh> <order> <steps> <train waves,comma> <test waves,comma> <k_list> <r_list>
"""
import sys, json, time, platform, hashlib
from pathlib import Path
import numpy as np
import ngsolve
from ngsolve import *
from team13_model import *
from transient import Model, waveform, T_END

ROOT = Path(__file__).resolve().parent
NU_REF = NU_0 / 1000.0


def bits(n, idx):
    b = BitArray(n); b.Clear()
    for i in idx: b.Set(int(i))
    return b


class Extender:
    def __init__(self, m):
        self.m = m; n = m.fes.ndof
        self.air_inv = m.k_rest.mat.Inverse(bits(n, m.idx_air), inverse='pardiso')
        u, v = m.fes.TnT(); steel = m.mesh.Materials('steel')
        self.k_ref = BilinearForm(NU_REF * curl(u) * curl(v) * dx(steel) + SIGMA_PLACEHOLDER / m.dt * u * v * dx(steel),
                                  symmetric=True).Assemble()
        self.int_inv = self.k_ref.mat.Inverse(bits(n, m.idx_int), inverse='pardiso')
        self.tmp = m.gfu.vec.CreateVector(); self.rhs = m.gfu.vec.CreateVector(); self.sol = m.gfu.vec.CreateVector()

    def _solve(self, mat, inv, g_full):
        self.tmp.FV().NumPy()[:] = g_full; mat.Mult(self.tmp, self.rhs); self.rhs.data *= -1
        self.sol.data = inv * self.rhs; return self.sol.FV().NumPy().copy()

    def extend(self, g):
        """Full vector from a trace g (len n_gamma): air + steel-interior harmonic extensions."""
        m = self.m; e = np.zeros(m.fes.ndof); e[m.idx_gamma] = g
        ea = self._solve(m.k_rest.mat, self.air_inv, e); ei = self._solve(self.k_ref.mat, self.int_inv, e)
        e[m.idx_air] = ea[m.idx_air]; e[m.idx_int] = ei[m.idx_int]; return e

    def interior_lift(self, g):
        m = self.m; e = np.zeros(m.fes.ndof); e[m.idx_gamma] = g
        return self._solve(self.k_ref.mat, self.int_inv, e)[m.idx_int]

    def source(self):
        """Air solution for unit current with zero trace and zero steel values."""
        m = self.m; self.sol.data = self.air_inv * m.f.vec; e = np.zeros(m.fes.ndof); e[m.idx_air] = self.sol.FV().NumPy()[m.idx_air]
        return e


def pod(blocks, rank):
    X = np.concatenate([b / np.linalg.norm(b) for b in blocks], axis=1)
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    return U[:, :rank], float(np.sum(s[:rank] ** 2) / np.sum(s ** 2))


def build_basis(m, ext, train, k, r, cache):
    """train: list of snapshot arrays (steps+1, ndof)."""
    Vg, rg = pod([S[1:, m.idx_gamma].T for S in train], k)
    if 'lift' not in cache:
        cache['lift'] = [np.array([ext.interior_lift(s[m.idx_gamma]) for s in S[1:]]).T for S in train]
    R = [S[1:, m.idx_int].T - L for S, L in zip(train, cache['lift'])]
    Wi, ri = pod(R, r) if r else (np.zeros((len(m.idx_int), 0)), 1.0)
    cols = [ext.extend(Vg[:, j]) for j in range(k)]
    for j in range(Wi.shape[1]):
        e = np.zeros(m.fes.ndof); e[m.idx_int] = Wi[:, j]; cols.append(e)
    return np.array(cols).T, dict(gamma_retention=rg, interior_retention=ri)


def run_rom(m, Phi, e_src, wave, steps, tol=1e-10):
    fn = waveform(wave); dt = T_END / steps; n, nr = Phi.shape
    g = m.gfu; q = np.zeros(nr); r = g.vec.CreateVector(); col = g.vec.CreateVector(); out = g.vec.CreateVector()
    snaps = [np.zeros(n)]; m.An.vec[:] = 0; its = []; t0 = time.time()
    def setvec(vec, qq, i): vec.FV().NumPy()[:] = Phi @ qq + i * e_src
    for k in range(1, steps + 1):
        i = fn(k * dt); m.An.vec.FV().NumPy()[:] = snaps[-1]
        for it in range(40):
            setvec(g.vec, q, i); m.a.AssembleLinearization(g.vec); m.residual(g.vec, i, r)
            rr = Phi.T @ r.FV().NumPy(); JP = np.empty((n, nr))
            for j in range(nr):
                col.FV().NumPy()[:] = Phi[:, j]; m.a.mat.Mult(col, out); JP[:, j] = out.FV().NumPy()
            Jr = Phi.T @ JP; dq = np.linalg.solve(0.5 * (Jr + Jr.T), rr); dec = abs(dq @ rr)
            if it == 0: dec0 = max(dec, 1e-300)
            if dec <= tol * dec0 or dec < 1e-24: break
            setvec(g.vec, q, i); E0 = m.energy(g.vec, i); tau = 1.
            while True:
                setvec(out, q - tau * dq, i)
                if m.energy(out, i) <= E0 or tau < 1e-4: break
                tau *= .5
            q = q - tau * dq
        else: raise RuntimeError(('ROM Newton', k))
        its.append(it); snaps.append(Phi @ q + i * e_src)
    return np.array(snaps), dict(newton_mean=float(np.mean(its)), newton_max=int(max(its)), seconds=time.time() - t0)


def steel_error(m, A, B):
    num = den = 0.
    for a, b in zip(A[1:], B[1:]):
        num += m.steel_B_norm2(a - b); den += m.steel_B_norm2(b)
    return float(np.sqrt(num / den))


if __name__ == '__main__':
    maxh, order, steps = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    train_w, test_w = sys.argv[4].split(','), sys.argv[5].split(',')
    ks = [int(v) for v in sys.argv[6].split(',')]; rs = [int(v) for v in sys.argv[7].split(',')]
    half = len(sys.argv) > 8 and sys.argv[8] == 'half'
    fdir = ROOT / 'results' / (f'fom_h{maxh*1000:g}_p{order}_n{steps}' + ('_half' if half else ''))
    with TaskManager():
        m = Model(maxh, order, T_END / steps, half=half); ext = Extender(m); e_src = ext.source()
        train = [np.load(fdir / f'{w}.npy') for w in train_w]; test = {w: np.load(fdir / f'{w}.npy') for w in test_w}
        cache = {}; rows = []
        for k in ks:
            for r in rs:
                Phi, info = build_basis(m, ext, train, k, r, cache)
                for w, ref in test.items():
                    u, st = run_rom(m, Phi, e_src, w, steps)
                    err = steel_error(m, u, ref)
                    rows.append(dict(k=k, r=r, wave=w, steel_B_error=err, **info, **st))
                    print(f'k{k} r{r} {w}: steel B error {100*err:.3f}%  newton {st["newton_mean"]:.1f}/{st["newton_max"]}  {st["seconds"]:.0f}s', flush=True)
    res = dict(case='team13_surface_rom', half=half, maxh_steel=maxh, order=order, steps=steps, train=train_w, test=test_w,
               ndof=m.fes.ndof, n_gamma=len(m.idx_gamma), n_interior=len(m.idx_int), nu_ref=NU_REF, rows=rows,
               host=platform.node(), ngsolve=ngsolve.__version__,
               sources={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in ['team13_model.py', 'transient.py', 'rom.py']})
    (fdir / f'rom_{"-".join(train_w)}.json').write_text(json.dumps(res, indent=2))
