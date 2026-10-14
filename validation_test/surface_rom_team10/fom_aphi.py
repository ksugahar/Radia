"""A-phi FOM with snapshot output (TEAM_CASE=team10, TEAM_QUARTER=1 expected).
Usage: python fom_aphi.py <maxh> <order> <wave> <steps> [max_newton]
Output: results/team10_aphiF_h{maxh}_p{order}_n{steps}[_quarter][_nw{k}][_T{ms}ms]/{wave}.npy|json
"""
import sys, os, json, platform, hashlib
from pathlib import Path
import numpy as np
import ngsolve
from ngsolve import *
from team13_model import *
from transient import waveform, T_END
from aphi_model import ModelAphi, run_fom

ROOT = Path(__file__).resolve().parent

if __name__ == '__main__':
    maxh, order, wave, steps = float(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
    max_newton = int(sys.argv[5]) if len(sys.argv) > 5 else 40
    if os.environ.get('NGS_THREADS'): SetNumThreads(int(os.environ['NGS_THREADS']))
    with TaskManager():
        m = ModelAphi(maxh, order, T_END / steps)
        snaps, hist, info = run_fom(m, waveform(wave), steps, max_newton=max_newton)
    tag = (f'team10_aphiF_h{maxh*1000:g}_p{order}_n{steps}' + ('_quarter' if QUARTER else '_half') + ('_coil2' if COIL2 else '')
           + (f'_nw{max_newton}' if max_newton < 40 else '') + (f'_T{T_END*1000:g}ms' if abs(T_END - 0.12) > 1e-12 else ''))
    out = ROOT / 'results' / tag; out.mkdir(parents=True, exist_ok=True)
    np.save(out / f'{wave}.npy', snaps)
    res = dict(case=f'team10_aphi_fom_{wave}', formulation='A-phi, HCurl nograds x H1(steel)', maxh_steel=maxh, order=order,
               steps=steps, dt=T_END / steps, quarter=QUARTER, gauge=m.gauge, max_newton=max_newton, sigma=SIGMA_PLACEHOLDER,
               ndof=m.fes.ndof, ndof_A=m.Vh.ndof, ndof_phi=m.Qh.ndof, n_gamma=len(m.idx_gamma), n_interior=len(m.idx_int),
               n_air=len(m.idx_air), solve_s=info['seconds'], linear_solves=info['linear_solves'], history=hist,
               platform_class=platform.system(), ngsolve=ngsolve.__version__,
               linear_residual_limit=info['linear_residual_limit'],
               refinement_solves=info['refinement_solves'],
               max_relative_linear_residual=info['max_relative_linear_residual'],
               sources={q: hashlib.sha256((ROOT / q).read_bytes()).hexdigest() for q in ['team13_model.py', 'aphi_model.py', 'fom_aphi.py', 'transient.py']})
    (out / f'{wave}.json').write_text(json.dumps(res, indent=2))
    print('done', wave, 'ndof', m.fes.ndof, 'gamma', len(m.idx_gamma), 'int', len(m.idx_int), 'air', len(m.idx_air), 'sec', round(info['seconds'], 1))
