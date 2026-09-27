"""Automatic selection of the training excitation patterns (TEAM_CASE=team10, TEAM_QUARTER=1, TEAM_COIL2=1).
Start from the main-coil rise only, then add one training trajectory per round from a candidate pool.
Selection rules compared at equal budget:
  residual : RB-greedy.  Run the current ROM on every candidate and pick the largest FOM residual,
             eta = sqrt(sum_k R_k^T K0^-1 R_k / sum_k F_k^T K0^-1 F_k),  K0 = FOM tangent at u = 0
             (one factorisation; no candidate FOM is solved before it is selected)
  gp       : BayPOD-like.  GP (RBF) on the pattern features (a1, a2, shape); pick the largest posterior
             variance.  Physics-blind: it depends only on where the selected features lie.
  hand     : S2 of run_pattern_study.py (c1_rise, c2_pulse), for reference.
Candidate snapshots are cheap FOMs (30 steps, Newton 1) computed beforehand; a rule only reads a
candidate's snapshots after selecting it.  c1_rise_neg is the exact mirror -c1_rise (odd B-H law), so a
sound indicator must never pick it.
After every round: steel B error on the exact tests (c2_rise, mix_add, mix_oppose) for surface_lift and
global bases, and eta on the tests (checks that eta tracks the true error).
Usage: python run_pattern_select.py [maxh] [order] [rounds]
"""
import os, sys, json, time, platform, hashlib
from pathlib import Path
import numpy as np
import ngsolve
from ngsolve import *
from team13_model import *
from transient import waveform, T_END
from aphi_model import ModelAphi
from rom import steel_error, build_basis
from rom_aphi import ExtenderAphi, FastROMAphi, split_orthonormalize, build_basis_global, bits
from ecsw import contributions_fast, nnls_early, NumpyHR

assert CASE == 'team10' and QUARTER and COIL2, 'set TEAM_CASE=team10 TEAM_QUARTER=1 TEAM_COIL2=1'
ROOT = Path(__file__).resolve().parent; R = ROOT / 'results'
FEAT = {'c1_rise': (1, 0, 0), 'c1_pulse': (1, 0, 1), 'c2_pulse': (0, 1, 1), 'mixp_add': (1, 1, 1),
        'mixp_oppose': (1, -1, 1), 'c2_pulse_x2': (0, 2, 1), 'c1_rise_half': (.5, 0, 0), 'c1_rise_neg': (-1, 0, 0)}
START = 'c1_rise'; HAND = ['c1_rise', 'c2_pulse']; TESTS = ['c2_rise', 'mix_add', 'mix_oppose']
K, RR = 16, 8; TOL = 3e-3; STEPS = 60; GP_LEN = 1.0   # evaluation basis (tests)
# selection basis: only has to rank the candidates, so it is kept small (10 % target passes at k4 r2)
_sb = os.environ.get('SELECT_BASIS', '4,2,1e-2').split(','); K_SEL, R_SEL, TOL_SEL = int(_sb[0]), int(_sb[1]), float(_sb[2])
# accuracy target 10 %.  Stop threshold on the steel indicator: not calibrated yet (the old 0.05 was for the dual norm)
ETA_STOP = float(os.environ.get('SELECT_ETA_STOP', '0.05')); STRIDE = 5
CHECK = os.environ.get('SELECT_CHECK', '1') == '1'   # round 0: compare steel/dual indicators and strides, time each


class Indicator:
    """Dual-norm FOM residual of a ROM trajectory, normalised by the source term.
    kind='dual' : full residual, norm K0^-1 (K0 = FOM tangent at u = 0; one full factorisation, full Apply per step)
    kind='steel': for u = Phi q + E i with the surface_lift basis the air rows of the residual vanish by construction
                  (harmonic extension, E solves the air with the source), so only the steel rows (gamma + interior)
                  are assembled: steel Variations + K_rest u on the gamma rows.  Norm: the steel lifting operator
                  k_ref (57k DOF factorisation).  Normalisation: the steel load of the bare source solution E c,
                  (K_rest E c)|steel, which is linear in c -> 2x2 Gram matrix.
    """
    def __init__(self, m, ext, E, kind='steel'):
        self.m = m; self.kind = kind; n = m.fes.ndof
        self.g = m.gfu.vec.CreateVector(); self.r = m.gfu.vec.CreateVector(); self.z = m.gfu.vec.CreateVector()
        if kind == 'dual':
            m.gfu.vec[:] = 0; m.gn.vec[:] = 0; m.a.AssembleLinearization(m.gfu.vec)
            self.Kinv = m.a.mat.Inverse(m.fes.FreeDofs(), inverse='pardiso'); loads = [fs.vec.FV().NumPy().copy() for fs in m.fs]
        else:
            (u, p), _ = m.fes.TnT(); w, _ = energy_density(); steel = m.mesh.Materials('steel'); s = SIGMA_PLACEHOLDER
            An = m.gn.components[0]; reg = 1e-12 * s / m.dt
            self.a_s = BilinearForm(m.fes, symmetric=True)
            self.a_s += Variation(0.5 * s / m.dt * (u - An + grad(p)) * (u - An + grad(p)) * dx(steel))
            self.a_s += Variation(w(sqrt(1e-12 + curl(u) * curl(u))) * dx(steel, bonus_intorder=2))
            self.a_s += Variation(0.5 * m.gauge * NU_0 * u * u * dx(steel))
            self.a_s += Variation(0.5 * reg * p * p * dx(steel))
            self.idx = np.sort(np.concatenate([m.idx_gamma, m.idx_int]))
            self.Kinv = ext.k_ref.mat.Inverse(bits(n, self.idx), inverse='pardiso')
            self.mask = np.zeros(n); self.mask[self.idx] = 1.
            loads = []
            for j in range(E.shape[1]):
                self.g.FV().NumPy()[:] = E[:, j]; m.k_rest.mat.Mult(self.g, self.r); loads.append(self.r.FV().NumPy() * self.mask)
        Y = []
        for l in loads: self.r.FV().NumPy()[:] = l; self.z.data = self.Kinv * self.r; Y.append(self.z.FV().NumPy().copy())
        self.G = np.array([[float(np.dot(y, l)) for l in loads] for y in Y])

    def residual(self, u_prev, u, c):
        m = self.m; m.gn.vec.FV().NumPy()[:] = u_prev; self.g.FV().NumPy()[:] = u
        if self.kind == 'dual': m.residual(self.g, c, self.r)
        else:
            self.a_s.Apply(self.g, self.r); m.k_rest.mat.Mult(self.g, self.z)
            self.r.FV().NumPy()[:] = (self.r.FV().NumPy() + self.z.FV().NumPy()) * self.mask
        return self.r

    def __call__(self, U, fn, stride=1):
        """stride > 1 samples every stride-th step (the residual still uses the true previous step)."""
        m = self.m; num = den = 0.
        for k in range(stride, U.shape[0], stride):
            c = m.currents(fn(k * m.dt)); r = self.residual(U[k - 1], U[k], c)
            self.z.data = self.Kinv * r; num += InnerProduct(self.z, r); den += float(c @ self.G @ c)
        return float(np.sqrt(max(num, 0.) / den))

    def air_fraction(self, U, fn, k):
        """|full residual on air rows| / |full residual| at step k (should be ~0 for the surface_lift basis)."""
        m = self.m; c = m.currents(fn(k * m.dt)); m.gn.vec.FV().NumPy()[:] = U[k - 1]; self.g.FV().NumPy()[:] = U[k]
        m.residual(self.g, c, self.r); r = self.r.FV().NumPy(); fr = np.array([m.fes.FreeDofs()[i] for i in range(len(r))])
        return float(np.linalg.norm(r[m.idx_air]) / np.linalg.norm(r[fr]))


def mem():
    try:
        import psutil; p = psutil.Process(); vm = psutil.virtual_memory()
        return f'rss {p.memory_info().rss/2**30:.1f} GB, private {p.memory_info().private/2**30:.1f} GB, avail {vm.available/2**30:.1f} GB'
    except Exception as ex: return repr(ex)


def gp_pick(selected, pool):
    X = np.array([FEAT[w] for w in selected], float); kern = lambda A, B: np.exp(-((A[:, None] - B[None]) ** 2).sum(-1) / (2 * GP_LEN ** 2))
    Kinv = np.linalg.inv(kern(X, X) + 1e-8 * np.eye(len(X))); var = {}
    for w in pool:
        x = np.array([FEAT[w]], float); ks = kern(X, x)[:, 0]; var[w] = float(1 - ks @ Kinv @ ks)
    return max(var, key=var.get), var


if __name__ == '__main__':
    maxh = float(sys.argv[1]) if len(sys.argv) > 1 else .006; order = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    if os.environ.get('NGS_THREADS'): SetNumThreads(int(os.environ['NGS_THREADS']))
    base = f'team10_aphiF_h{maxh*1000:g}_p{order}'; log = []; t0 = time.time()
    with TaskManager():
        m = ModelAphi(maxh, order, T_END / STEPS); ext = ExtenderAphi(m); E = ext.sources(); Z = np.zeros_like(E); ind = Indicator(m, ext, E, 'steel')
        ind_dual = Indicator(m, ext, E, 'dual') if CHECK else None
        refs = {w: np.load(R / f'{base}_n60_quarter_coil2' / f'{w}.npy') for w in TESTS}
        tdir = R / f'{base}_n30_quarter_coil2_nw1'; snaps_cache = {}; lift_cache = {}

        def snaps(w):
            if w not in snaps_cache:
                snaps_cache[w] = -snaps(START) if w == 'c1_rise_neg' else np.load(tdir / f'{w}.npy')
            return snaps_cache[w]

        def lifts(w):
            if w not in lift_cache:
                lift_cache[w] = np.array([ext.interior_lift(s[m.idx_gamma]) for s in snaps(w)[1:]]).T
            return lift_cache[w]

        def make_rom(sel, basis, k=K, r=RR, tol=TOL):
            train = [snaps(w) for w in sel]
            if basis == 'surface_lift': P, info = build_basis(m, ext, train, k, r, {'lift': [lifts(w) for w in sel]}); es = E
            else: P, info = build_basis_global(m, train, k + r); es = Z
            Phi, oi = split_orthonormalize(m, P)
            C, d, sel_el = contributions_fast(m, Phi, [S[j] for S in train for j in range(1, S.shape[0])])
            els, wts, _ = nnls_early(C, d, tol)
            return FastROMAphi(m, Phi, es, hr=NumpyHR(m, Phi, sel_el[els], wts, 2 * m.order + 2)), dict(info, **oi, ecsw_elements=len(els))

        def evaluate(rule, rnd, sel, eta_pool=None, extra=None):
            for basis in ['surface_lift', 'global']:
                rom, info = make_rom(sel, basis); row = dict(rule=rule, round=rnd, selected=list(sel), basis=basis, info=info, tests={})
                for w in TESTS:
                    try:
                        u, st = rom.run(waveform(w), STEPS)
                        row['tests'][w] = dict(steel_B_error=steel_error(m, u, refs[w]), eta=ind(u, waveform(w), STRIDE) if basis == 'surface_lift' else None,  # steel indicator assumes zero air residual
                                              seconds=st['seconds'])
                    except Exception as ex: row['tests'][w] = dict(failure=repr(ex))
                if basis == 'surface_lift' and eta_pool is not None: row['eta_pool'] = eta_pool
                if extra: row.update(extra)
                log.append(row)
                print(f'{rule:8s} r{rnd} {basis:12s} {"+".join(sel):40s} ' + ' '.join(
                    f'{w}={100*v["steel_B_error"]:.2f}%' + (f'/eta{v["eta"]:.3g}' if v['eta'] is not None else '') if 'steel_B_error' in v else f'{w}=FAIL' for w, v in row['tests'].items()), flush=True)

        # residual-greedy; stops when every candidate has eta < ETA_STOP (or after `rounds` additions)
        sel = [START]
        for rnd in range(rounds + 1):
            pool = [w for w in FEAT if w not in sel]; eta = {}; chk = {}; tm = dict(build=0., rom=0., eta=0.)
            if rnd < rounds:
                t = time.time(); rom, _ = make_rom(sel, 'surface_lift', K_SEL, R_SEL, TOL_SEL); tm['build'] = time.time() - t
                for w in pool:
                    try:
                        t = time.time(); u, _ = rom.run(waveform(w), STEPS); tm['rom'] += time.time() - t
                        t = time.time(); eta[w] = ind(u, waveform(w), STRIDE); tm['eta'] += time.time() - t
                        if rnd == 0 and CHECK:
                            c = {}
                            for name, fn_eta in [('steel_s1', lambda: ind(u, waveform(w))), ('dual_s5', lambda: ind_dual(u, waveform(w), STRIDE)),
                                                 ('dual_s1', lambda: ind_dual(u, waveform(w)))]:
                                t = time.time(); c[name] = fn_eta(); c[name + '_s'] = time.time() - t
                            c['air_fraction_k30'] = ind.air_fraction(u, waveform(w), 30); chk[w] = c
                    except Exception as ex:
                        import traceback; eta[w] = float('inf'); print('   candidate', w, 'FAILED', repr(ex), mem(), flush=True); traceback.print_exc()
                print('  mem', mem(), flush=True)
                print('  eta pool', {w: round(v, 4) for w, v in eta.items()}, {k: round(v, 1) for k, v in tm.items()}, flush=True)
                for w, c in chk.items(): print('   check', w, {k: (round(v, 4) if not k.endswith('_s') else round(v, 1)) for k, v in c.items()}, flush=True)
            evaluate('residual', rnd, sel, eta_pool=eta, extra=dict(select_seconds=tm, indicator_check=chk))
            if rnd == rounds or max(eta.values()) < ETA_STOP: break
            sel = sel + [max(eta, key=eta.get)]
        # BayPOD-like GP variance on the pattern features
        if os.environ.get('SELECT_ONLY_RESIDUAL') == '1': rounds = -1
        sel = [START]
        for rnd in range(rounds + 1):
            pool = [w for w in FEAT if w not in sel]; pick, var = gp_pick(sel, pool)
            if rnd > 0: evaluate('gp', rnd, sel, extra=dict(gp_var=var))
            if rnd < rounds: sel = sel + [pick]
        if rounds >= 0: evaluate('hand', 1, HAND)
    res = dict(case='team10_two_coil_pattern_selection', start=START, features=FEAT, tests=TESTS, k=K, r=RR, ecsw_tol=TOL,
               gp_length=GP_LEN, eta_stop=ETA_STOP, stride=STRIDE, select_basis=[K_SEL, R_SEL, TOL_SEL], target=0.10, rows=log, ndof=m.fes.ndof, total_s=time.time() - t0, host=platform.node(), ngsolve=ngsolve.__version__,
               sources={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in
                        ['team13_model.py', 'transient.py', 'aphi_model.py', 'rom.py', 'rom_aphi.py', 'ecsw.py', 'run_pattern_select.py']})
    (R / f'team10_pattern_select_p{order}{os.environ.get("SELECT_TAG", "")}.json').write_text(json.dumps(res, indent=2))
