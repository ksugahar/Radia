"""Self-authored weak-coupling case; run on an admitted idle compute host.

Cold timings include assembly, table construction, final direct certification
and field export. No speed ratio is a correctness gate. Requires a fixed Radia
native runtime matching the transferred Python source; no editable install.
"""
import argparse
import json
from pathlib import Path
import sys
import time
from unittest.mock import patch

import numpy as np
import ngsolve as ng
from netgen.occ import Cylinder, Pnt, Vec, OCCGeometry, Glue

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src/radia/panels'))
from radia.panels import calc_inductance as calc


def run_case(directory, *, maxh=.008, current=20., repeats=2, threads=1, before_run=None):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    ng.SetNumThreads(threads)
    with ng.TaskManager():
        body = Cylinder(Pnt(0, 0, -.0125), Vec(0, 0, 1), r=.025, h=.025)
        for face in body.faces:
            face.name = 'sibc'
        body.mat('workpiece')
        fixture = ng.Mesh(OCCGeometry(Glue(list(body.faces))).GenerateMesh(maxh=maxh))
        # Independent analytic saturating law, not a measured/commercial dataset.
        h = np.r_[0., np.geomspace(.1, 1e6, 160)]
        bh = directory/'analytic-bh.txt'
        np.savetxt(bh, np.c_[h, 4e-7*np.pi*h + .02*np.tanh(h/200.)])
        angle = np.linspace(0, 2*np.pi, 181)
        path = np.stack([.035*np.cos(angle), .035*np.sin(angle), np.zeros_like(angle)], axis=1)
        source = dict(source_type='filament', paths=[np.stack([path[:-1], path[1:]], axis=1)],
                      I_fil=np.array([current+0j]))
        original_mesh = ng.Mesh
        output = []
        references = []
        for repeat in range(repeats):
            for mode in ('direct', 'table'):
                if before_run is not None:
                    before_run(repeat, mode)
                case = directory/f'{repeat}-{mode}'
                case.mkdir(exist_ok=True)
                label = str(case/'case.vol')
                fixture.ngmesh.Save(label)
                args = calc.build_argparser().parse_args([
                    '--coil-solver', 'peec', '--coil-step', 'unused.step', '--vol', label,
                    '--wp-label', 'sibc', '--frequency', '1000', '--current', str(current),
                    '--sigma', '5.8e7', '--mu-r', '1', '--h1-order', '1',
                    '--wp-bem-backend', 'intree-dense', '--impedance-model', 'esim',
                    '--bh-file', str(bh), '--half-thickness', '.025', '--esim-per-panel',
                    '--esim-max-iter', '50', '--esim-tol', '.001',
                    '--esim-panel-evaluator', mode])

                def load(value, *a, **kw):
                    return fixture if isinstance(value, str) and value == label else original_mesh(value, *a, **kw)

                started = time.perf_counter()
                with patch.object(ng, 'Mesh', load):
                    result = calc._solve_workpiece_weak_coupled(args, source)
                elapsed = time.perf_counter()-started
                zs = np.array(result['esim_per_panel_Z_s_real']) + 1j*np.array(result['esim_per_panel_Z_s_imag'])
                fields = np.array(result['esim_per_panel_H_t'])
                assert result['esim_converged']
                assert result['esim_fixed_point_relative_error'] <= args.esim_tol
                assert result['esim_panel_evaluation']['certification_cell_calls'] == len(zs)
                assert np.ptp(np.abs(zs))/np.mean(np.abs(zs)) > 1e-3
                if mode == 'direct':
                    references.append((result, zs, fields))
                reference, refzs, refh = references[-1]
                assert result['esim_panel_evaluation']['configuration_sha256'] == reference['esim_panel_evaluation']['configuration_sha256']
                error_z = float(np.max(abs(zs-refzs)/np.maximum(abs(refzs), 1e-30)))
                error_h = float(np.max(abs(fields-refh))/max(float(np.max(abs(refh))), 1e-30))
                error_p = abs(result['P_wp']/reference['P_wp']-1)
                assert max(error_z, error_h, error_p) < 2*args.esim_tol
                assert result['esim_iterations'] == reference['esim_iterations']
                assert result['esim_anderson_restarts'] == reference['esim_anderson_restarts']
                assert result['esim_anderson_clips'] == reference['esim_anderson_clips']
                history = [item['dZ_max'] for item in result['esim_history']]
                refhistory = [item['dZ_max'] for item in reference['esim_history']]
                assert np.max(np.abs(np.array(history)-refhistory)) < 2*args.esim_tol
                output.append(dict(repeat=repeat, mode=mode, cold_wall_seconds=elapsed,
                    panels=len(zs), iterations=result['esim_iterations'], history_dZ_max=history,
                    anderson_m=result['esim_anderson_m'], anderson_restarts=result['esim_anderson_restarts'],
                    anderson_clips=result['esim_anderson_clips'], power_W=result['P_wp'],
                    reaction_balance=result['wp_power_balance_relative_error'],
                    fixed_point_error=result['esim_fixed_point_relative_error'],
                    relative_Z_error=error_z, relative_H_error=error_h, relative_power_error=error_p,
                    evaluation=result['esim_panel_evaluation']))
                # Recover completed runs even if a later mode fails a gate.
                (directory/'completed-runs.json').write_text(json.dumps(output, indent=2,
                    allow_nan=False), encoding='utf-8')
        return dict(case='analytic-saturating-cylinder-filament', frequency_Hz=1000,
                    current_A=current, maxh_m=maxh, threads=threads, repeats=repeats,
                    ngsolve_version=ng.__version__, independent_fem_validation='not-performed', runs=output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--maxh', type=float, default=.008)
    opts = parser.parse_args()
    result = run_case(opts.job_dir, repeats=opts.repeats, threads=opts.threads, maxh=opts.maxh)
    opts.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
