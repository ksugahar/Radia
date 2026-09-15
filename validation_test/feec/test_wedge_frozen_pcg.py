"""Native true-residual regression for the saved difficult WEDGE tangent."""
import json
from pathlib import Path

import ngsolve as ng
import numpy as np
from ngsolve.meshes import MakeStructured3DMesh
from scipy.sparse import coo_matrix

from radia.vim import Solve
from radia.vim._solve import _h_solve_auto_prec


def test_frozen_wedge_jacobi_reaches_true_residual():
    evidence = Path(__file__).parents[1] / "esrf_three_engine/results/candidate_d0d0bc4b5"
    reference = json.loads((evidence / "lab-v5-frozen.json").read_text(encoding="utf-8"))
    tolerance = reference["inner_requested_tolerance"]
    with np.load(evidence / "lab-v5-frozen.npz", allow_pickle=False) as data:
        weight, demag, rhs, initial = (data[name].copy() for name in ("W", "N", "rhs", "initial_x"))
    ng.SetNumThreads(4)
    mesh = MakeStructured3DMesh(nx=1, ny=1, nz=1, hexes=False, prism=True,
                               mapping=lambda x, y, z: (.02*x, .02*y, .02*z))
    with ng.TaskManager():
        baseline = Solve(mesh, mu_r=100., H_ext=ng.CF((0., 0., 1000.)),
                         order=2, gram_eps=1e-9, tol=1e-9)
        operator = baseline["_charge_gram"]
        identity = np.eye(len(rhs))
        actual_n = np.column_stack([
            operator.apply_configured_demag(np.ascontiguousarray(identity[:, i]), True)
            for i in range(len(rhs))])
        assert np.linalg.norm(actual_n-demag)/np.linalg.norm(demag) < 1e-10
        result = _h_solve_auto_prec(operator, coo_matrix(weight), len(rhs), 1.,
                                    rhs, tolerance, 4000, x0=initial)
        solution = np.asarray(result["m"])
        action = weight @ solution + operator.apply_configured_demag(solution, True)
    residual = np.linalg.norm(rhs-action)/np.linalg.norm(rhs)
    report = dict(iterations=int(result["iters"]), tolerance=tolerance,
                  true_relative_residual=float(residual), timings=result["timings"])
    print(json.dumps(report, indent=2))
    assert result["timings"]["last_solve_converged"]
    assert np.isfinite(residual) and residual <= tolerance
    with ng.TaskManager():
        exhausted = _h_solve_auto_prec(operator, coo_matrix(weight), len(rhs), 1.,
                                       rhs, tolerance, 1, x0=initial)
        stopped = np.asarray(exhausted["m"])
        stopped_action = weight @ stopped + operator.apply_configured_demag(stopped, True)
    assert not exhausted["timings"]["last_solve_converged"]
    assert exhausted["iters"] == 1
    assert np.linalg.norm(rhs-stopped_action)/np.linalg.norm(rhs) > tolerance


if __name__ == "__main__":
    test_frozen_wedge_jacobi_reaches_true_residual()
