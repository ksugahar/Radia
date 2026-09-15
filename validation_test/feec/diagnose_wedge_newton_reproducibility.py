"""Read-only solver diagnostic using one installed wheel and a fixed wedge."""
import argparse
import hashlib
import importlib.metadata as md
import json
import os
import platform
import subprocess
from pathlib import Path
import time
import traceback

import ngsolve as ng
import numpy as np
import radia
from ngsolve.meshes import MakeStructured3DMesh
from radia.vim import Solve

p = argparse.ArgumentParser()
p.add_argument('--output', type=Path, required=True)
p.add_argument('--threads', type=int, default=4)
p.add_argument('--preconditioner', default='auto')
p.add_argument('--spectrum', action='store_true')
a = p.parse_args()
ng.SetNumThreads(a.threads)
native = Path(radia.__file__).parent / '_radia_pybind.pyd'
report = dict(version=md.version('radia'), import_path=radia.__file__,
              native_sha256=hashlib.sha256(native.read_bytes()).hexdigest(),
              dependencies={n: md.version(n) for n in ['numpy','scipy','ngsolve','netgen-mesher','mkl']},
              threads=a.threads, preconditioner=a.preconditioner,
              environment={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_THREADING_LAYER','MKL_CBWR']},
              passed=False)
mesh = MakeStructured3DMesh(nx=1,ny=1,nz=1,hexes=False,prism=True,
                            mapping=lambda x,y,z:(.02*x,.02*y,.02*z))
mesh_data = {'points':[list(v.point) for v in mesh.vertices],
             'elements':[[v.nr for v in el.vertices] for el in mesh.Elements(ng.VOL)]}
report['mesh_geometry'] = mesh_data
report['mesh_vertex_connectivity_sha256'] = hashlib.sha256(json.dumps(mesh_data,sort_keys=True).encode()).hexdigest()
report['mesh_identity_scope'] = 'Vertex/connectivity digest only; not a full .vol identity.'
report['mesh_recipe'] = dict(factory='ngsolve.meshes.MakeStructured3DMesh', nx=1, ny=1, nz=1,
                             hexes=False, prism=True, mapping_scale_m=0.02,
                             curvature_operation=None, identification_operation=None)
report['mesh_materials'] = list(mesh.GetMaterials())
report['mesh_boundaries'] = list(mesh.GetBoundaries())
root = Path(radia.__file__).parent
report['python_sources'] = {str(f.relative_to(root)).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest()
                            for f in sorted(root.rglob('*.py'))}
h = np.array([0.,10.,100.,1000.,10000.,100000.])
bh = np.column_stack((h,4e-7*np.pi*(h+1e6*h/(h+300.))))
report['bh_table'] = bh.tolist()
report['controls'] = dict(drive=300000., order=2, gram_eps=1e-9, tol=1e-9, nl_tol=2e-5, nl_maxit=80)
start = time.perf_counter()
print('solve start: ' + report['version'], flush=True)
try:
    with ng.TaskManager():
        r = Solve(mesh,bh_table=bh,H_ext=ng.CF((0.,0.,300000.)),order=2,
                  gram_eps=1e-9,tol=1e-9,nl_tol=2e-5,nl_maxit=80,
                  preconditioner=a.preconditioner)
    report['stats'] = r['nonlinear_solve_stats']
    report['M_avg'] = np.asarray(r['M_avg']).tolist()
    residual = float(report['stats']['nonlinear_final_relative_residual'])
    report['passed'] = bool(report['stats']['nonlinear_converged_final_stage'] and np.isfinite(residual) and residual <= 2e-5)
except Exception as exc:
    print('solve exception: ' + str(exc), flush=True)
    report['error'] = str(exc)
    report['traceback'] = traceback.format_exc()
    frames = []
    tb = exc.__traceback__
    while tb:
        frames.append(tb.tb_frame)
        tb = tb.tb_next
    for frame in frames:
        loc = frame.f_locals
        if frame.f_code.co_name == '_solve_nonlinear_energy_cpp':
            report['outer_relative_residual'] = float(loc['relative_residual'])
            report['outer_stats'] = dict(loc['stats'])
        if frame.f_code.co_name == '_solve_W':
            report['inner_timings'] = dict(loc['res'].get('timings', {}))
            report['inner_iterations'] = int(loc['it'])
            report['inner_requested_tolerance'] = float(loc['solve_tol'])
    # The failed solve returns its iterate. Evaluate its actual W+N residual
    # without changing the solve or accepting that iterate as a solution.
    outer = next((f.f_locals for f in frames if f.f_code.co_name == '_solve_nonlinear_energy_cpp'), None)
    inner = next((f.f_locals for f in frames if f.f_code.co_name == '_solve_W'), None)
    if outer is not None and inner is not None:
        print('evaluating returned inner iterate', flush=True)
        from scipy.sparse import coo_matrix
        rows,cols,values = inner['W_matrix'].COO()
        x = np.asarray(inner['res']['m'], dtype=float)
        w = coo_matrix((values,(rows,cols)),shape=(len(x),len(x))).tocsr()
        with ng.TaskManager():
            action = w @ x + outer['_N_apply'](x)
        rhs = np.asarray(inner['rhs'])
        report['inner_true_relative_residual'] = float(np.linalg.norm(rhs-action)/np.linalg.norm(rhs))
        if a.spectrum:
            if len(x) > 500:
                raise ValueError('Dense diagnostic limited to 500 unknowns')
            from scipy.linalg import eigvalsh
            with ng.TaskManager():
                eye = np.eye(len(x))
                demag = np.column_stack([outer['_N_apply'](eye[:,i]) for i in range(len(x))])
            system = w.toarray() + demag
            symmetric = (system + system.T)*0.5
            diagonal = np.diag(symmetric)
            report['frozen_system'] = dict(
                relative_asymmetry=float(np.linalg.norm(system-system.T)/np.linalg.norm(system)),
                material_eigenvalues=eigvalsh((w.toarray()+w.toarray().T)*0.5).tolist(),
                system_eigenvalues=eigvalsh(symmetric).tolist(),
                diagonal_min=float(diagonal.min()), diagonal_max=float(diagonal.max()))
            if np.all(diagonal > 0):
                scale = np.sqrt(diagonal)
                report['frozen_system']['jacobi_eigenvalues'] = eigvalsh(symmetric/scale[:,None]/scale[None,:]).tolist()
report['seconds'] = time.perf_counter()-start
report['platform'] = dict(platform=platform.platform(), processor=platform.processor(), python=platform.python_version())
from threadpoolctl import threadpool_info
report['threadpools'] = threadpool_info()
for pool in report['threadpools']:
    dll = Path(pool['filepath'])
    pool['sha256'] = hashlib.sha256(dll.read_bytes()).hexdigest()
cpu = subprocess.check_output(['pwsh','-NoProfile','-Command',
    'Get-CimInstance Win32_Processor | Select-Object Name,NumberOfCores,NumberOfLogicalProcessors | ConvertTo-Json -Compress'], text=True)
report['cpu'] = json.loads(cpu)
a.output.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
print(json.dumps(report,indent=2,allow_nan=False))
raise SystemExit(0 if report['passed'] else 1)
