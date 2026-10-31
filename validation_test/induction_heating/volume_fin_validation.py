"""Public rectilinear fin diagnostics: finite-volume R/L and regional pairings.

This is NOT an adopted fin validation: matched converged volumetric FEM and
the user-owned four 1% gates remain to be established. No hidden SIBC/return
geometry is supplied. Run under an explicitly bounded caller TaskManager.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import ngsolve as ng
import radia
from radia import peec_matrices
from radia.fin_observables import FinMeasurementRegions
from radia.volume_peec import build_volume_peec
from radia.volume_fin import solve_volume_fin, measure_volume_fin


def public_fin(level):
    """Independently authored T conductor; exact union of supplied SI boxes."""
    if isinstance(level, bool) or not isinstance(level, int) or level < 1:
        raise ValueError('Positive integer refinement level required')
    def subdivide(knots):
        return np.r_[*(np.linspace(a, b, level+1)[:-1] for a, b in zip(knots[:-1], knots[1:])), knots[-1]]*1e-3
    x = subdivide([-3., -1., 0., 1., 3.])
    y = subdivide([-1., 0., .5, 1., 1.5])
    z = subdivide([-.25, .25])
    xc, yc, zc = np.meshgrid((x[:-1]+x[1:])/2, (y[:-1]+y[1:])/2,
                              (z[:-1]+z[1:])/2, indexing='ij')
    metal = (yc < 0) | (abs(xc) < .001)
    source = [(0, j, k, 0, -1) for j in range(level) for k in range(level)]
    sink = [(len(x)-2, j, k, 0, 1) for j in range(level) for k in range(level)]
    return build_volume_peec((x, y, z), metal, 5.8e7,
                             source_faces=source, sink_faces=sink)


def campaign(output, levels, frequencies, quadrature_order):
    root = Path(__file__).resolve().parents[2]
    head = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
    regions = FinMeasurementRegions((0., 0., 0.), (0., 1., 0.), (1., 0., 0.),
                                     .0002, .0012, (-.0008, .0008))
    convention = ('Free-space nonmagnetic volume Neumann partial potential; '
                  'no return conductor in this diagnostic. Whole conductor incl. '
                  'outside-window participates; half each reciprocal mutual pair.')
    artifacts = {}
    for module in (radia, peec_matrices):
        path = Path(module.__file__).resolve()
        if not path.is_relative_to(root):
            raise RuntimeError(f'Wrong worktree import: {path}')
        artifacts[module.__name__] = dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    report = dict(status='INCOMPLETE', source=head, artifacts=artifacts,
                  reference_status='matched-converged-FEM-pending',
                  adoption_status='not-adopted', physics_target_relative=.01,
                  physics_verdict=None, quadrature_order=quadrature_order,
                  geometry='independently-authored-public-rectilinear-T', cases=[])
    def save():
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    save()
    with ng.TaskManager():
        for level in levels:
            topology = public_fin(level)
            for frequency in frequencies:
                started = time.perf_counter()
                solved = solve_volume_fin(topology, frequency_hz=frequency,
                                           drive_current_peak_A=1., quadrature_order=quadrature_order)
                solve_s = time.perf_counter()-started
                measured = measure_volume_fin(topology, solved, regions,
                                               quadrature_order=quadrature_order, convention=convention)
                report['cases'].append(dict(level=level, cells=len(topology.cell_indices),
                    branches=len(topology.branches), frequency_hz=frequency,
                    R_ohm=solved.resistance_ohm, L_partial_H=solved.partial_inductance_H,
                    Z_ohm=[solved.impedance_ohm.real, solved.impedance_ohm.imag],
                    KCL=solved.kcl_relative_residual, KVL=solved.kvl_relative_residual,
                    solve_s=solve_s, total_s=time.perf_counter()-started, regions=measured))
                save()  # preserve raw values before any reviewer judgment
    report['status'] = 'DIAGNOSTICS_COMPLETE_NOT_ADOPTED'
    save()
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--levels', nargs='+', type=int, required=True)
    parser.add_argument('--frequencies', nargs='+', type=float, required=True)
    parser.add_argument('--quadrature-order', type=int, required=True)
    args = parser.parse_args()
    campaign(args.output, args.levels, args.frequencies, args.quadrature_order)
