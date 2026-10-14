"""Calibrate the stop rule of the residual-greedy pattern selection (TEAM_CASE=team10, TEAM_QUARTER=1, TEAM_COIL2=1).
At the 10 % target the production basis is the selection basis, so both are the same small surface_lift ROM.
Training sets: c1_rise plus every subset of 0-2 pool members (29 sets).  For each set and basis:
  tests      : true steel B error of c2_rise, mix_add, mix_oppose (exact 60-step FOMs) and their eta
  pool       : eta of every pool member not in the set (the quantity a stop rule can see)
  floor      : eta of the set's own training patterns (truncation floor of the ROM)
Rules then checked offline (analyse()):  A  max(pool eta) < tau;   B  max(pool eta) / max(floor) < rho.
A rule is safe at a threshold if it never stops on a set whose worst test error is >= 10 %.
Usage: python run_stop_calibration.py [maxh] [order]
"""
import os, sys, json, time, platform, hashlib, itertools
from pathlib import Path
import numpy as np
import ngsolve
from ngsolve import *
from team13_model import *
from transient import waveform, T_END
from aphi_model import ModelAphi
from rom import steel_error, build_basis
from rom_aphi import ExtenderAphi, FastROMAphi, split_orthonormalize
from ecsw import contributions_fast, nnls_early, NumpyHR
from run_pattern_select import Indicator, FEAT, START, TESTS, STEPS, STRIDE

ROOT = Path(__file__).resolve().parent; R = ROOT / 'results'
POOL = [w for w in FEAT if w != START]; BASES = [(4, 2, 1e-2), (8, 4, 1e-2)]; TARGET = 0.10


def analyse(rows):
    out = {}
    for k, r, _ in BASES:
        rs = [x for x in rows if x['k'] == k and x['r'] == r and x['pool']]
        fail = np.array([x['worst_test'] >= TARGET for x in rs]); A = np.array([max(x['pool'].values()) for x in rs])
        B = A / np.array([max(x['floor'].values()) for x in rs])
        res = {}
        for name, v in [('A_max_pool_eta', A), ('B_ratio_to_floor', B)]:
            safe = float(v[fail].min()) if fail.any() else float('inf')     # stop only below every failing set
            res[name] = dict(safe_threshold=safe, stops_at_safe=int((v < safe).sum()), passing_sets=int((~fail).sum()), sets=len(rs),
                             values_fail=sorted(map(float, v[fail])), values_pass=sorted(map(float, v[~fail])))
        out[f'k{k}r{r}'] = res
    return out


if __name__ == '__main__':
    maxh = float(sys.argv[1]) if len(sys.argv) > 1 else .006; order = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    if os.environ.get('NGS_THREADS'): SetNumThreads(int(os.environ['NGS_THREADS']))
    base = f'team10_aphiF_h{maxh*1000:g}_p{order}'; rows = []; t0 = time.time()
    with TaskManager():
        m = ModelAphi(maxh, order, T_END / STEPS); ext = ExtenderAphi(m); E = ext.sources(); ind = Indicator(m, ext, E, 'steel')
        refs = {w: np.load(R / f'{base}_n60_quarter_coil2' / f'{w}.npy') for w in TESTS}
        tdir = R / f'{base}_n30_quarter_coil2_nw1'; snaps, lifts = {}, {}
        for w in FEAT:
            snaps[w] = -np.load(tdir / f'{START}.npy') if w == 'c1_rise_neg' else np.load(tdir / f'{w}.npy')
            lifts[w] = np.array([ext.interior_lift(s[m.idx_gamma]) for s in snaps[w][1:]]).T
        sets = [[START] + list(c) for n in range(3) for c in itertools.combinations(POOL, n)]
        for sel in sets:
            train = [snaps[w] for w in sel]
            for k, r, tol in BASES:
                t = time.time()
                P, info = build_basis(m, ext, train, k, r, {'lift': [lifts[w] for w in sel]}); Phi, oi = split_orthonormalize(m, P)
                C, d, sel_el = contributions_fast(m, Phi, [S[j] for S in train for j in range(1, S.shape[0])])
                els, wts, _ = nnls_early(C, d, tol, log=lambda *a, **kw: None)
                rom = FastROMAphi(m, Phi, E, hr=NumpyHR(m, Phi, sel_el[els], wts, 2 * m.order + 2)); t_build = time.time() - t
                row = dict(selected=sel, k=k, r=r, ecsw_tol=tol, ecsw_elements=len(els), build_s=t_build, tests={}, pool={}, floor={})
                try:
                    for w in TESTS:
                        u, _ = rom.run(waveform(w), STEPS); row['tests'][w] = dict(steel_B_error=steel_error(m, u, refs[w]), eta=ind(u, waveform(w), STRIDE))
                    for w in FEAT:
                        u, _ = rom.run(waveform(w), STEPS); v = ind(u, waveform(w), STRIDE)
                        (row['floor'] if w in sel else row['pool'])[w] = v
                    row['worst_test'] = max(x['steel_B_error'] for x in row['tests'].values())
                except Exception as ex: row['failure'] = repr(ex)
                row['seconds'] = time.time() - t; rows.append(row)
                if 'failure' in row: print(f'k{k}r{r} {"+".join(sel):45s} FAILED {row["failure"]}', flush=True); continue
                print(f'k{k}r{r} {"+".join(sel):45s} worst {100*row["worst_test"]:5.2f}%  max pool eta {max(row["pool"].values()) if row["pool"] else float("nan"):.4f}'
                      f'  floor {max(row["floor"].values()):.4f}  {row["seconds"]:.0f}s', flush=True)
    ana = analyse([x for x in rows if 'failure' not in x])
    print(json.dumps(ana, indent=1), flush=True)
    res = dict(case='team10_two_coil_stop_calibration', target=TARGET, bases=BASES, stride=STRIDE, pool=POOL, tests=TESTS, rows=rows,
               analysis=ana, ndof=m.fes.ndof, total_s=time.time() - t0, platform_class=platform.system(), ngsolve=ngsolve.__version__,
               sources={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in
                        ['team13_model.py', 'transient.py', 'aphi_model.py', 'rom.py', 'rom_aphi.py', 'ecsw.py', 'run_pattern_select.py', 'run_stop_calibration.py']})
    (R / f'team10_stop_calibration_p{order}.json').write_text(json.dumps(res, indent=2))
