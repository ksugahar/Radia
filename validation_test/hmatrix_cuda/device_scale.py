"""Device compressed H-matrix setup and solve at scale on public spheres; JSON record, no speed or scale threshold.

Each mesh is a public sphere, checked by the authoritative check-vol through
``tools/prepare_hdiv_mesh_receipt.py`` run in ``--receipt-python``. The only pass
condition is the native true residual of the solve. Timings are diagnostic.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import platform
import subprocess
import time


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def public_sphere(path, maxh, curve_order):
    import ngsolve as ng
    from netgen.csg import CSGeometry, Pnt, Sphere
    geometry = CSGeometry()
    geometry.Add(Sphere(Pnt(0, 0, 0), 1).mat("magnetic_body").bc("exterior"))
    mesh = ng.Mesh(geometry.GenerateMesh(maxh=maxh))
    if curve_order == 2:
        mesh.Curve(2)
    mesh.ngmesh.Save(str(path))


def checked_inputs(path, contract, receipt_python):
    repo = Path(__file__).resolve().parents[2]
    receipt = Path(str(path) + ".receipt.json")
    subprocess.run([str(receipt_python), str(repo / "tools/prepare_hdiv_mesh_receipt.py"), "--mesh", str(path),
                    "--label-contract", str(contract), "--output", str(receipt)], check=True)
    provenance = json.loads(Path(str(receipt) + ".provenance.json").read_text(encoding="utf-8"))
    return dict(label_contract=str(contract), mesh_check_receipt=str(receipt),
                mesh_check_receipt_sha256=sha256(receipt), mesh_check_provenance=provenance)


def run_case(path, args, checked):
    from hdiv_process_memory import ProcessMemory
    from radia.hdiv_gpu import PreparedHDivLinear
    with ProcessMemory(os.getpid()) as memory:
        start = time.perf_counter()
        with PreparedHDivLinear(path, backend="gpu", order=args.order, gram_compression="aca",
                                curve_order=args.curve_order, **checked) as prepared:
            setup_seconds = time.perf_counter() - start
            start = time.perf_counter()
            solution = prepared.solve((0., 0., args.field), args.inv_chi, tolerance=args.tolerance,
                                      max_iterations=args.max_iterations, preconditioner=args.preconditioner)
            solve_seconds = time.perf_counter() - start
            stats = prepared.stats()
    residual = float(solution["relative_residual"])
    return dict(passed=residual <= args.tolerance, relative_residual=residual, iterations=int(solution["iterations"]),
                preconditioner=args.preconditioner, dofs=stats["dofs"], family=stats["family"],
                gram_aca=stats["gram_aca"], hmatrix=stats["assembly"]["hmatrix"],
                setup_peak_device_bytes=stats["assembly"]["setup_peak_device_bytes"],
                operator_owned_device_bytes=stats["operator"]["owned_device_bytes"],
                solver_owned_device_bytes=solution["solver_owned_device_bytes"],
                mean_magnetization=[float(value) for value in solution["mean_magnetization"]],
                setup_seconds_diagnostic=setup_seconds, solve_seconds_diagnostic=solve_seconds,
                peak_process_private_bytes=max(sample["private_bytes"] for sample in memory.samples),
                mesh_sha256=sha256(path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--maxh", type=float, nargs="+", required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--receipt-python", type=Path, required=True)
    parser.add_argument("--order", type=int, choices=[0, 1, 2], default=1)
    parser.add_argument("--curve-order", type=int, choices=[1, 2], default=1)
    parser.add_argument("--field", type=float, default=1e5)
    parser.add_argument("--inv-chi", type=float, default=.1)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--max-iterations", type=int, default=2000)
    parser.add_argument("--preconditioner", default="mass-diagonal")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    import ngsolve as ng
    result = dict(schema="radia.hdiv-device-scale.v1", accepted=False, units="SI", backend="gpu",
                  scope="Device-built compressed Gram, B, B^T and mass; no matrix crosses; timings are diagnostic",
                  order=args.order, curve_order=args.curve_order, ngsolve_version=ng.__version__,
                  machine=platform.machine(), source_sha=os.environ.get("RADIA_SOURCE_SHA"),
                  cpu_affinity=sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
                  cases=[])
    try:
        for maxh in args.maxh:
            path = args.out / f"sphere-{maxh:g}.vol".replace(".", "p", 1)
            with ng.TaskManager():
                public_sphere(path, maxh, args.curve_order)
            result["cases"].append(dict(run_case(path, args, checked_inputs(path, args.contract, args.receipt_python)),
                                        maxh=maxh))
        result["accepted"] = all(case["passed"] for case in result["cases"])
    except BaseException as exc:
        result["error"] = repr(exc)
        raise
    finally:
        (args.out / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")


if __name__ == "__main__":
    main()
