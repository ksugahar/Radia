"""Audit the closed-torus notebook's matrix and basis-dependent excitation.

This is a diagnostic, not an inductance acceptance test. No .vol is written.
"""
from __future__ import annotations
import argparse
import ast
import hashlib
import json
import platform
import time
from pathlib import Path
import numpy as np


def run(notebook, output):
    import ngsolve as ng
    from ngsolve.bem import LaplaceSL
    ng.SetNumThreads(2)
    raw = Path(notebook).read_bytes()
    nb = json.loads(raw)
    source = "\n".join("".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code")
    tree = ast.parse(source)
    definitions = [n for n in tree.body if isinstance(n, ast.FunctionDef)
                   and n.name in {"create_torus_mesh", "extract_dense"}]
    namespace = {"np": np}
    exec(compile(ast.Module(body=definitions, type_ignores=[]), str(notebook), "exec"), namespace)
    started = time.perf_counter()
    with ng.TaskManager():
        mesh = namespace["create_torus_mesh"](0.05, 0.005, 0.5)
        space = ng.HDivSurface(mesh, order=0)
        if space.ndof > 1500:
            raise RuntimeError(f"Diagnostic exceeds bounded dense budget: {space.ndof}")
        label = sorted(set(mesh.GetBoundaries()))[0]
        trial, test = space.TnT()
        operator = LaplaceSL(trial.Trace()*ng.ds(label))*test.Trace()*ng.ds(label)
        matrix = namespace["extract_dense"](operator.mat, space.ndof)*4e-7*np.pi
    assert np.isfinite(matrix).all()
    singular = np.linalg.svd(matrix, compute_uv=False)
    tol = singular[0]*max(matrix.shape)*np.finfo(float).eps
    rank = int(np.count_nonzero(singular > tol))
    e = np.ones(space.ndof)/space.ndof
    signs = np.where(np.arange(space.ndof)%2, -1., 1.)
    transformed = signs[:,None]*matrix*signs[None,:]
    def extract(a, load):
        current = np.linalg.lstsq(a, load, rcond=None)[0]
        return float(1/(load@current)), float(np.linalg.norm(a@current-load)/np.linalg.norm(load))
    legacy, residual = extract(matrix, e)
    covariant, covariant_residual = extract(transformed, signs*e)
    naive, naive_residual = extract(transformed, e)
    np.testing.assert_allclose(covariant, legacy, rtol=1e-8)
    report = {
        "schema": "radia.closed-torus-diagnostic.v1", "host": platform.node(),
        "python": platform.python_version(), "ngsolve": ng.__version__, "numpy": np.__version__,
        "notebook_sha256": hashlib.sha256(raw).hexdigest(),
        "diagnostic_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "radius_m": 0.05, "wire_radius_m": 0.005, "curvaturesafety": 0.5,
        "threads": 2, "ndof": space.ndof, "free_dofs": sum(space.FreeDofs()),
        "boundaries": list(mesh.GetBoundaries()), "selected_boundary": label,
        "rank": rank, "rank_tolerance": float(tol),
        "singular_min": float(singular[-1]), "singular_max": float(singular[0]),
        "relative_symmetry_error": float(np.linalg.norm(matrix-matrix.T)/np.linalg.norm(matrix)),
        "zero_columns": int(np.count_nonzero(np.linalg.norm(matrix, axis=0)==0)),
        "legacy_L_H": legacy, "legacy_relative_residual": residual,
        "basis_transformed_covariant_L_H": covariant,
        "basis_transformed_covariant_relative_residual": covariant_residual,
        "basis_transformed_ones_L_H": naive,
        "basis_transformed_ones_relative_residual": naive_residual,
        "elapsed_s": time.perf_counter()-started,
        "scope": "One coarse matrix diagnostic. e=ones/n is a coefficient-space constraint, not a specified physical loop current. No general closed-surface nullspace conclusion or physical inductance acceptance.",
    }
    Path(output).write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--notebook", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.notebook, args.output)
