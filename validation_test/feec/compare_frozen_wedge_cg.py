"""Compare uninterrupted and periodically restarted CG on identical saved data."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve
from scipy.sparse.linalg import LinearOperator, cg


def compare(system_path, reference_path, maxiter=4000, restart_period=1000):
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    tolerance = float(reference["inner_requested_tolerance"])
    with np.load(system_path, allow_pickle=False) as saved:
        matrix = saved["W"] + saved["N"]
        rhs = saved["rhs"].copy()
        initial = saved["initial_x"].copy()
        returned = saved["returned_x"].copy()
    diagonal = np.diag(matrix)
    if not np.all(diagonal > 0) or not np.isfinite(matrix).all():
        raise ValueError("Finite matrix and positive Jacobi diagonal required")
    scale = np.linalg.norm(rhs)
    if not scale > 0:
        raise ValueError("Nonzero source required")
    preconditioner = LinearOperator(matrix.shape, matvec=lambda x: x / diagonal)
    counter = [0]

    def count(_):
        counter[0] += 1

    solution, info = cg(matrix, rhs, x0=initial, rtol=tolerance, atol=0,
                        maxiter=maxiter, M=preconditioner, callback=count)
    uninterrupted = dict(info=int(info), iterations=counter[0],
        true_relative_residual=float(np.linalg.norm(rhs-matrix@solution)/scale))
    solution = initial.copy()
    history = []
    completed = 0
    while completed < maxiter:
        budget = min(restart_period, maxiter-completed)
        counter[0] = 0
        solution, info = cg(matrix, rhs, x0=solution, rtol=tolerance, atol=0,
                            maxiter=budget, M=preconditioner, callback=count)
        completed += counter[0]
        residual = float(np.linalg.norm(rhs-matrix@solution)/scale)
        history.append(dict(iterations=completed, info=int(info),
                            true_relative_residual=residual))
        if info <= 0 or counter[0] == 0:
            break
    direct = solve((matrix+matrix.T)/2, rhs, assume_a="pos")
    return dict(schema="radia.frozen-cg-diagnostic.v1",
        system_sha256=hashlib.sha256(system_path.read_bytes()).hexdigest(),
        tolerance=tolerance, maxiter=maxiter, restart_period=restart_period,
        uninterrupted=uninterrupted, restarted=history,
        native_iterate_dense_residual=float(np.linalg.norm(rhs-matrix@returned)/scale),
        direct_true_residual=float(np.linalg.norm(rhs-matrix@direct)/scale),
        scope="Saved dense system diagnostic; not native solver acceptance.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.system, args.reference)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps(result, indent=2))
