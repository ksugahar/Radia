"""Smallest basis meeting the 10 % target on the two-coil tests (TEAM_CASE=team10, TEAM_QUARTER=1, TEAM_COIL2=1).
Training = the automatically selected set (c1_rise + c2_pulse_x2, one greedy round of run_pattern_select.py).
Sweep (k surface, r interior) and the ECSW tolerance; report steel B error and online seconds per test.
Usage: python run_basis_10pct.py [maxh] [order]
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
from rom_aphi import ExtenderAphi, FastROMAphi, split_orthonormalize, build_basis_global
from ecsw import contributions_fast, nnls_early, NumpyHR

assert CASE == 'team10' and QUARTER and COIL2, 'set TEAM_CASE=team10 TEAM_QUARTER=1 TEAM_COIL2=1'
ROOT = Path(__file__).resolve().parent; R = ROOT / 'results'
TRAIN = ['c1_rise', 'c2_pulse_x2']; TESTS = ['c2_rise', 'mix_add', 'mix_oppose']; STEPS = 60
KR = [(4, 2), (6, 3), (8, 4), (10, 5), (16, 8)]; TOLS = [3e-3, 1e-2, 3e-2]

if __name__ == '__main__':
    maxh = float(sys.argv[1]) if len(sys.argv) > 1 else .006; order = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    if os.environ.get('NGS_THREADS'): SetNumThreads(int(os.environ['NGS_THREADS']))
    base = f'team10_aphiF_h{maxh*1000:g}_p{order}'; rows = []; t0 = time.time()
    with TaskManager():
        m = ModelAphi(maxh, order, T_END / STEPS); ext = ExtenderAphi(m); E = ext.sources(); Z = np.zeros_like(E)
        refs = {w: np.load(R / f'{base}_n60_quarter_coil2' / f'{w}.npy') for w in TESTS}
        train = [np.load(R / f'{base}_n30_quarter_coil2_nw1' / f'{w}.npy') for w in TRAIN]
        snaps = [S[j] for S in train for j in range(1, S.shape[0])]; cache = {}
        for k, r in KR:
            for basis in ['surface_lift', 'global']:
                if basis == 'surface_lift': P, info = build_basis(m, ext, train, k, r, cache); es = E
                else: P, info = build_basis_global(m, train, k + r); es = Z
                Phi, oi = split_orthonormalize(m, P); info = dict(info, **oi)
                t = time.time(); C, d, sel = contributions_fast(m, Phi, snaps); t_c = time.time() - t
                for tol in TOLS:
                    t = time.time(); els, wts, _ = nnls_early(C, d, tol); t_n = time.time() - t
                    hr = NumpyHR(m, Phi, sel[els], wts, 2 * m.order + 2)
                    for w in TESTS:
                        try:
                            u, st = FastROMAphi(m, Phi, es, hr=hr).run(waveform(w), STEPS)
                            err, sec, fail = steel_error(m, u, refs[w]), st['seconds'], None
                        except Exception as ex: err, sec, fail = None, None, repr(ex)
                        rows.append(dict(basis=basis, k=k, r=r, columns=Phi.shape[1], ecsw_tol=tol, ecsw_elements=len(els),
                                         contributions_s=t_c, nnls_s=t_n, test=w, steel_B_error=err, online_s=sec, failure=fail, **info))
                        print(f'k{k} r{r} {basis:12s} tol{tol:g} el{len(els):5d} {w:10s}: ' +
                              (f'{100*err:.2f}% {sec:.1f}s' if err is not None else f'FAILED {fail}'), flush=True)
    res = dict(case='team10_two_coil_basis_10pct', target=0.10, train=TRAIN, tests=TESTS, rows=rows, ndof=m.fes.ndof,
               total_s=time.time() - t0, host=platform.node(), ngsolve=ngsolve.__version__,
               sources={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in
                        ['team13_model.py', 'transient.py', 'aphi_model.py', 'rom.py', 'rom_aphi.py', 'ecsw.py', 'run_basis_10pct.py']})
    (R / f'team10_basis_10pct_p{order}.json').write_text(json.dumps(res, indent=2))
