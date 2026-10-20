"""Independent CPU energy-Newton fixture for the bounded public BH experiment."""
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--field', type=float, required=True, help='Uniform z-directed field in A/m')
    parser.add_argument('--ho-far-factor', type=float, default=2.0, help='CPU near/far threshold; inf explicitly selects all-high-order reference')
    parser.add_argument('--nonlinear-tolerance', type=float, default=1e-10)
    parser.add_argument('--linear-tolerance', type=float, default=1e-10)
    parser.add_argument('--max-iterations', type=int, default=50)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    if not (args.ho_far_factor>0 and 0<args.nonlinear_tolerance<=1e-6 and 0<args.linear_tolerance<=1e-6 and 1<=args.max_iterations<=100 and 1<=args.threads<=16):
        parser.error('invalid integration/convergence/thread bounds')
    if not np.isfinite(args.field):
        parser.error('field must be finite')
    args.out.mkdir(parents=True, exist_ok=False)
    import ngsolve as ng
    from radia.vim import Solve
    from hdiv_hdf5 import write_hdiv_h5
    ng.SetNumThreads(args.threads)
    h = np.array([0., 1., 10., 100., 1000., 1e4, 1e5])
    # Independently authored synthetic monotone saturation law; not measured steel.
    table = np.column_stack((h, 4e-7*np.pi*(h+1e6*1000*h/(1e6+1000*h))))
    result = dict(accepted=False, scope='Bounded synthetic BH independent CPU fixture', units='SI')
    try:
        with ng.TaskManager():
            mesh = ng.Mesh(str(args.mesh))
            start = time.perf_counter()
            ans = Solve(mesh, bh_table=table, H_ext=ng.CoefficientFunction((0,0,args.field)),
                        order=2, gram_eps=1e-10, leaf=64, eta=2, tol=args.linear_tolerance, ho_far_factor=args.ho_far_factor,
                        nl_tol=args.nonlinear_tolerance, nl_maxit=args.max_iterations, newton_inner_tol='fixed',
                        newton_warmstart='linear', preconditioner='mass-riesz', curve_gauss=8)
            seconds = time.perf_counter()-start
            write_hdiv_h5(args.out/'cpu.h5', dict(table=table,
                coefficients=ans['_m_coefficients'], mean=np.asarray(ans['M_avg']),
                applied_field=np.array([0.,0.,args.field])),
                dict(units='SI', labels=list(mesh.GetBoundaries())))
            result.update(accepted=True, whole_seconds_diagnostic=seconds,
                          stats=ans['nonlinear_solve_stats'], threads=args.threads,
                          integration_profile=('all-high-order' if np.isinf(args.ho_far_factor) else 'near-far-split'),
                          ho_far_factor=('inf' if np.isinf(args.ho_far_factor) else args.ho_far_factor),
                          mesh_sha256=hashlib.sha256(args.mesh.read_bytes()).hexdigest())
    except BaseException as ex:
        result['error'] = repr(ex)
        raise
    finally:
        (args.out/'result.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
