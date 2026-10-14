"""Changing excitation PATTERN: two independently driven coils (TEAM_CASE=team10, TEAM_QUARTER=1, TEAM_COIL2=1).
Training sets (cheap snapshots, 30 steps, Newton 1):
  S1 = main coil only (c1_rise, c1_pulse)             -> the web coil pattern is never seen
  S2 = main coil only (c1_rise) + web coil only (c2_pulse) -> both patterns seen, never together
Tests (exact 60-step FOMs): c2_rise, mix_add, mix_oppose (nonlinear: no superposition).
Bases at equal size (k surface + r interior -> A/phi block split):
  surface_lift  trace POD + exact air + nu_ref continuation, u = Phi q + E i   (E: air solution per coil)
  surface_zero  trace POD + exact air + zero continuation,  u = Phi q + E i
  global        plain POD of full snapshots,                u = Phi q
  global_src    plain POD + the exact source columns,       u = Phi q + E i   (separates source handling from trace structure)
Usage: python run_pattern_study.py [maxh] [order]
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
from rom_aphi import ExtenderAphi, FastROMAphi, split_orthonormalize, build_basis_zero_lift, build_basis_global
from ecsw import contributions_fast, nnls_early, NumpyHR

assert CASE == 'team10' and QUARTER and COIL2, 'set TEAM_CASE=team10 TEAM_QUARTER=1 TEAM_COIL2=1'
ROOT = Path(__file__).resolve().parent; R = ROOT / 'results'
TRAIN = {'S1': ['c1_rise', 'c1_pulse'], 'S2': ['c1_rise', 'c2_pulse']}
TESTS = ['c2_rise', 'mix_add', 'mix_oppose']; KR = [(10, 5), (16, 8)]; TOL = 3e-3; STEPS = 60


def sections_history(m, U):
    g = GridFunction(m.fes); out = []
    for u in U[1:]:
        g.vec.FV().NumPy()[:] = u; out.append(team10_sections(m.bfield(g), m.mesh))
    return out


if __name__ == '__main__':
    maxh = float(sys.argv[1]) if len(sys.argv) > 1 else .006; order = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    if os.environ.get('NGS_THREADS'): SetNumThreads(int(os.environ['NGS_THREADS']))
    base = f'team10_aphiF_h{maxh*1000:g}_p{order}'; rows = []
    with TaskManager():
        m = ModelAphi(maxh, order, T_END / STEPS); ext = ExtenderAphi(m); E = ext.sources(); Z = np.zeros_like(E)
        assert E.shape[1] == 2
        refs = {w: np.load(R / f'{base}_n60_quarter_coil2' / f'{w}.npy') for w in TESTS}
        tdir = R / f'{base}_n30_quarter_coil2_nw1'
        for s, waves in TRAIN.items():
            train = [np.load(tdir / f'{w}.npy') for w in waves]
            snaps = [S[j] for S in train for j in range(1, S.shape[0])]
            for k, r in KR:
                bases = {}
                P, info = build_basis(m, ext, train, k, r, {}); bases['surface_lift'] = (P, info, E)
                P, info = build_basis_zero_lift(m, ext, train, k, r); bases['surface_zero'] = (P, info, E)
                P, info = build_basis_global(m, train, k + r); bases['global'] = (P, info, Z); bases['global_src'] = (P, info, E)
                for name, (P0, info, es) in bases.items():
                    Phi, oi = split_orthonormalize(m, P0); info = dict(info, **oi)
                    C, d, sel = contributions_fast(m, Phi, snaps)
                    els, wts, res = nnls_early(C, d, TOL); hr = NumpyHR(m, Phi, sel[els], wts, 2 * m.order + 2)
                    for w in TESTS:
                        for variant, h in [('all_elements', None), (f'ecsw_{TOL:g}', hr)]:
                            try:
                                u, st = FastROMAphi(m, Phi, es, hr=h).run(waveform(w), STEPS)
                                err = steel_error(m, u, refs[w]); sec = sections_history(m, u); fail = None
                            except Exception as ex:
                                err, sec, st, fail = None, None, dict(seconds=None), repr(ex)
                            rows.append(dict(set=s, train=waves, basis=name, k=k, r=r, columns=Phi.shape[1], test=w, variant=variant,
                                             steel_B_error=err, seconds=st['seconds'], ecsw_elements=len(els) if h is not None else None,
                                             failure=fail, sections=sec, **info))
                            print(f'{s} k{k} r{r} {name:13s} {w:10s} {variant:13s}: ' + (f'{100*err:.3f}%' if err is not None else f'FAILED {fail}'), flush=True)
    res = dict(case='team10_two_coil_pattern_study', train=TRAIN, tests=TESTS, rows=rows, ndof=m.fes.ndof,
               fom_sections={w: [h['B_pos'] for h in json.loads((R / f'{base}_n60_quarter_coil2' / f'{w}.json').read_text())['history']] for w in TESTS},
               platform_class=platform.system(), ngsolve=ngsolve.__version__,
               sources={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in
                        ['team13_model.py', 'transient.py', 'aphi_model.py', 'rom.py', 'rom_aphi.py', 'ecsw.py', 'run_pattern_study.py']})
    (R / f'team10_pattern_study_p{order}.json').write_text(json.dumps(res, indent=2))
