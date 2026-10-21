"""Public synthetic ESIM cell refinement evidence, not full-coil acceptance."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy.special import ive
from radia.esim_cell_problem import ESIMFiniteSlabSolver, ESIMCellProblemSolver


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    mu0 = 4e-7*np.pi
    linear = [[h, 100*mu0*h] for h in (0., 10., 100., 1e3, 1e4, 1e5, 5e5)]
    fields = np.r_[0., np.geomspace(1., 1e6, 35)]
    nonlinear = [[h, mu0*h+1.2*np.tanh(h/1500.)] for h in fields]
    rows = []
    # Acceptance is declared before executing: 0.3% vs analytic, 0.3% seed
    # variation, with an empirical 0.1% nested mesh indicator (not a proof).
    gates = dict(analytic_relative_limit=.003, seed_relative_limit=.003,
                 mesh_indicator_limit=.001, true_linear_residual_limit=1e-6)
    for frequency in (150e3, 222e3, 300e3):
        for geometry in ('slab', 'cylinder', 'halfspace'):
            for nodes in ((100, 800, 1000) if geometry == 'halfspace' else (24, 60, 200)):
                if geometry == 'halfspace':
                    cell = ESIMCellProblemSolver(bh_curve=linear, sigma=5e6, frequency=frequency, n_nodes=nodes)
                    expected = (1+1j)*np.sqrt(2*np.pi*frequency*100*mu0/(2*5e6))
                else:
                    cell = ESIMFiniteSlabSolver(.005, mu_r=100, sigma=5e6, frequency=frequency, n_nodes=nodes, geometry=geometry)
                    gamma = (1+1j)/cell.delta
                    expected = cell.rho*gamma*(ive(1, gamma*.005)/ive(0, gamma*.005) if geometry == 'cylinder' else np.tanh(gamma*.005))
                result = cell.solve(1000.)
                error = abs(result['Z']/expected-1)
                rows.append(dict(geometry=geometry, frequency_Hz=frequency, seed_nodes=nodes,
                    Z_real=float(result['Z'].real), Z_imag=float(result['Z'].imag),
                    P_W_per_m2=result['P_prime'], analytic_relative_error=float(error),
                    accepted=bool(error <= gates['analytic_relative_limit']),
                    mesh_history=result['mesh_history'], depth_history=result.get('depth_history', []),
                    linear_residual=result['linear_residual']))
    nonlinear_rows = []
    for field in (200., 1000., 4000., 8000., 20000., 50000., 100000.):
        values = []
        for nodes in (24, 60, 200):
            result = ESIMFiniteSlabSolver(.005, bh_curve=nonlinear, sigma=5e6,
                frequency=222e3, n_nodes=nodes, geometry='cylinder').solve(field)
            values.append(dict(seed_nodes=nodes, Z_real=float(result['Z'].real),
                Z_imag=float(result['Z'].imag), P_W_per_m2=result['P_prime'],
                history=result['mesh_history'], iterations=result['iterations']))
        z = [complex(v['Z_real'], v['Z_imag']) for v in values]
        error = max(abs(v/z[-1]-1) for v in z)
        power_error = max(abs(v['P_W_per_m2']/values[-1]['P_W_per_m2']-1) for v in values)
        nonlinear_rows.append(dict(H_peak_A_per_m=field, values=values,
            impedance_seed_relative_error=error, power_seed_relative_error=power_error,
            accepted=bool(max(error, power_error) <= gates['seed_relative_limit'])))
    iteration_limits = []
    for field in (50000., 100000.):
        cell = ESIMFiniteSlabSolver(.005, bh_curve=nonlinear, sigma=5e6,
            frequency=222e3, n_nodes=60, geometry='cylinder')
        try:
            limited = cell.solve(field, max_iter=50)
            iteration_limits.append(dict(field_A_per_m=field, cap=50,
                status='converged', iterations=limited['iterations']))
        except RuntimeError as error:
            iteration_limits.append(dict(field_A_per_m=field, cap=50,
                status='failed_loudly', error=str(error)))
    source_root = Path(__file__).resolve().parents[2]
    sources = {p: hashlib.sha256((source_root/p).read_bytes()).hexdigest() for p in
        ('src/radia/esim_cell_problem.py', 'src/radia/esim_depth.py')}
    report = dict(gates=gates, linear=rows, synthetic_nonlinear=nonlinear_rows,
        iteration_limits=iteration_limits,
        accepted=all(r['accepted'] for r in rows+nonlinear_rows),
        elapsed_seconds=time.perf_counter()-started, source_hashes=sources,
        versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
        convention='SI; H is peak phasor; P=0.5 Re(Z) |H|^2',
        limitations=['Synthetic materials only; private reported case not rerun.',
            'Independent volumetric FEM/full-coil acceptance NOT_PERFORMED.',
            'Phasor ESIM does not certify time-resolved nonlinear BH cycles.',
            'Nested mesh change is empirical and not a rigorous error bound.'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    if not report['accepted']:
        raise SystemExit('Adaptive ESIM acceptance failed; inspect result JSON')


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        import sys
        if '--output' in sys.argv:
            output = Path(sys.argv[sys.argv.index('--output')+1])
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(json.dumps(dict(accepted=False, error=repr(error),
                status='FAILED; incomplete matrix is not PASS'), indent=2)+'\n', encoding='utf-8')
        raise
