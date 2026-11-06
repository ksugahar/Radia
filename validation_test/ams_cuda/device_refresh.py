"""Device AMS refresh against the CPU HypreBasedAMS reuse_hierarchy Update.

A public in-memory Netgen unit cube; NGSolve assembles the HCurl operators. The
first operator builds the hierarchy (CPU, reuse_hierarchy). Each update keeps the
sparsity pattern and changes the coefficients spatially, so every Galerkin level,
l1 norm and coarsest inverse changes. The device refresh is compared level by
level with the CPU Update of the same values; PCG on the refreshed operator checks
the coarsest inverse. Records only; acceptance criteria stay with the user.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time


def main(args):
    import numpy as np
    import ngsolve as ng
    from netgen.csg import unit_cube
    source_root = Path(__file__).resolve().parents[2]
    sys.path[:0] = [str(source_root), str(source_root / 'packages/radia-simulink/src')]
    from radia_simulink.ams_exchange import write_exchange
    import matlab.engine
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=False)
    build = args.build.resolve(); matlab_dir = args.matlab.resolve()
    receipt = json.loads((build / 'build.json').read_text())
    binary = build / 'radia_cuda_ams.mexw64'
    digest = hashlib.sha256(binary.read_bytes()).hexdigest()
    if digest != receipt['sha256']['build-ams-cuda/radia_cuda_ams.mexw64']:
        raise RuntimeError('AMS binary receipt mismatch')
    here = Path(__file__)
    report = dict(schema='radia.ams-device-refresh-validation.v1', completed=False,
                  scope='dimensionless linear 3D HCurl unit cube; device refresh vs CPU reuse Update',
                  native_source_commit=receipt['source_commit'], binary_sha256=digest,
                  helper_sha256=hashlib.sha256(here.read_bytes()).hexdigest(),
                  matlab_helper_sha256=hashlib.sha256(here.with_name('deviceRefresh.m').read_bytes()).hexdigest(),
                  ngsolve=ng.__version__, maxh=args.maxh, cycle=args.cycle)
    eng = None
    try:
        ng.SetNumThreads(1)
        mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=args.maxh))
        fes = ng.HCurl(mesh, order=1, nograds=True, dirichlet='.*'); u, v = fes.TnT()
        free = np.asarray(list(fes.FreeDofs()), dtype=bool)
        if not (64 <= free.sum() and fes.ndof < 50000):
            raise RuntimeError('Case must be nontrivial and below the frozen DoF cap')
        x, y, z = ng.x, ng.y, ng.z
        coefficients = [(1.0, 1.0),
                        (1.0 + 0.8 * ng.sin(3 * x) * ng.cos(2 * y), 0.3 + 4 * z * z),
                        (2.5 + ng.cos(5 * z), 1e-3 + x * y)]

        def assemble(nu, sigma):
            form = ng.BilinearForm(fes); form += (nu * ng.curl(u) * ng.curl(v) + sigma * u * v) * ng.dx
            form.Assemble(); return form.mat

        def csr(matrix):
            values, columns, rowptr = matrix.CSR()
            return dict(rows=float(matrix.height), cols=float(matrix.width),
                        rowptr=np.array(rowptr, dtype=np.int32, copy=True),
                        columns=np.array(columns, dtype=np.int32, copy=True), values=np.array(values, copy=True))

        matrices = [csr(assemble(*c)) for c in coefficients]
        for m in matrices[1:]:
            if not (np.array_equal(m['rowptr'], matrices[0]['rowptr']) and np.array_equal(m['columns'], matrices[0]['columns'])):
                raise RuntimeError('Update operators must keep the first sparsity pattern')
        gradient, _ = fes.CreateGradient()
        coords = np.array([[mesh.ngmesh.Points()[i + 1][j] for j in range(3)] for i in range(mesh.nv)])
        model = dict(schema='radia.ams-native-input.v1', A=matrices[0], G=csr(gradient), coordinates=coords,
                     free=free.astype(np.int32), cycle=float(args.cycle))
        # The last update returns to the first operator: a second refresh on the same handle.
        updates = [m['values'] for m in matrices[1:]] + [matrices[0]['values']]
        write_exchange(out / 'input.h5', dict(model=model, updates=updates))
        report.update(ndof=fes.ndof, free_dofs=int(free.sum()), nonzeros=int(matrices[0]['values'].size),
                      vertices=mesh.nv, elements=mesh.ne,
                      coefficients=['nu=1, sigma=1', 'nu=1+0.8 sin(3x)cos(2y), sigma=0.3+4z^2',
                                    'nu=2.5+cos(5z), sigma=1e-3+xy', 'nu=1, sigma=1 (second refresh)'])
        start = time.perf_counter(); eng = matlab.engine.start_matlab('-nodesktop -nosplash')
        report['matlab_startup_seconds'] = time.perf_counter() - start
        eng.addpath(str(here.parent), nargout=0)
        eng.deviceRefresh(str(out), str(build), str(matlab_dir), nargout=0)
        report['gpu'] = json.loads((out / 'gpu-result.json').read_text())
        report['completed'] = True
    except BaseException as e:
        report['error'] = repr(e); raise
    finally:
        if eng is not None: eng.quit()
        (out / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps(report, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build', type=Path, required=True); p.add_argument('--matlab', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--maxh', type=float, default=.12); p.add_argument('--cycle', type=int, default=1, choices=(1, 7))
    main(p.parse_args())
