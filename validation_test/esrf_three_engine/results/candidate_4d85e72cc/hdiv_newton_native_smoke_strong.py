"""Check the installed candidate's real native HDiv nonlinear route."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time

import ngsolve as ng
import numpy as np
import radia
from ngsolve.meshes import MakeStructured3DMesh
from radia.vim import Solve

parser = argparse.ArgumentParser()
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--drive', type=float, default=1000.)
args = parser.parse_args()
ng.SetNumThreads(4)
mu0 = 4e-7*np.pi
linear_bh = np.array([[0., 0.], [1e3, mu0*100.*1e3], [1e5, mu0*100.*1e5]])
h = np.array([0., 10., 100., 1000., 10000., 100000.])
nonlinear_bh = np.column_stack((h, mu0*(h+1e6*h/(h+300.))))
report = {'radia': radia.__file__, 'version': importlib.metadata.version('radia'),
          'ngsolve': ng.__version__, 'drive': args.drive, 'rows': [], 'completed': False}
for kind in ('tet', 'hex', 'wedge'):
    for order in (1, 2):
        mesh = MakeStructured3DMesh(nx=1, ny=1, nz=1, hexes=kind=='hex', prism=kind=='wedge',
                                   mapping=lambda x,y,z: (.02*x, .02*y, .02*z))
        start = time.perf_counter()
        with ng.TaskManager():
            linear = Solve(mesh, mu_r=100., H_ext=ng.CF((0., 0., 1000.)), order=order,
                           gram_eps=1e-9, tol=1e-9)
            parity = Solve(mesh, bh_table=linear_bh, H_ext=ng.CF((0., 0., 1000.)), order=order,
                           gram_eps=1e-9, tol=1e-9, nl_tol=2e-5, nl_maxit=80)
            nonlinear = Solve(mesh, bh_table=nonlinear_bh, H_ext=ng.CF((0., 0., args.drive)),
                              order=order, gram_eps=1e-9, tol=1e-9, nl_tol=2e-5, nl_maxit=80)
        a,b = np.asarray(linear['M_avg']), np.asarray(parity['M_avg'])
        error = float(np.linalg.norm(a-b)/np.linalg.norm(a))
        stats = nonlinear['nonlinear_solve_stats']
        residual = float(stats['nonlinear_final_relative_residual'])
        assert error < 1e-6, (kind, order, error)
        assert np.isfinite(residual) and residual <= 2e-5, stats
        assert stats['nonlinear_converged_final_stage'], stats
        if args.drive >= 1e5:
            assert stats['nonlinear_tangent_assemblies'] > 0, stats
        row = dict(kind=kind, order=order, linear_parity=error, residual=residual,
                   seconds=time.perf_counter()-start, nonlinear_stats=stats)
        report['rows'].append(row)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
        print(json.dumps(row), flush=True)
report['completed'] = True
args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
