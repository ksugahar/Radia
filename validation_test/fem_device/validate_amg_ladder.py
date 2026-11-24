"""Device-built H1 AMG with CG against sparse Cholesky and the CPU compact AMG on identical inputs.

Each case assembles the same form on the device and on the CPU. The field reference is the sparse
Cholesky solution of the NGSolve matrix. Iterations of the device AMG sit beside those of the CPU
compact AMG, which uses the same strength threshold, level limit and coarsest size. Device
coarsening is PMIS and the CPU route's is its own, so the two hierarchies differ by design.

No criterion is adopted here. SOLUTION_REL is a PROPOSED detector. Writes one JSON with the source
commit, conditions, values, the transfer ledger and diagnostic timings. Usage:
    python validate_amg_ladder.py --out result.json [--cases cube/8 cube/16 jump/0.2 ...]
"""
import argparse
import json
from pathlib import Path
import platform
import subprocess
import time

import numpy as np
import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh

import radia.fem_gpu as route

ROOT = Path(__file__).resolve().parents[2]
SOLUTION_REL = 1e-6
TOLERANCE = 1e-10


def cube(n):
    return MakeStructured3DMesh(hexes=False, nx=n, ny=n, nz=n, mapping=lambda x, y, z: (x + 0.1 * y, y + 0.1 * z, z)), "left|right", 1.0


def jump(maxh):
    from netgen.occ import Box, Glue, OCCGeometry, Pnt, X

    first = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    first.mat("a")
    first.faces.Min(X).name = "left"
    second = Box(Pnt(1, 0, 0), Pnt(2, 1, 1))
    second.mat("b")
    second.faces.Max(X).name = "right"
    mesh = ng.Mesh(OCCGeometry(Glue([first, second])).GenerateMesh(maxh=maxh))
    return mesh, "left|right", {"a": 1.0, "b": 1000.0}


def run_case(text, min_coarse):
    family, value = text.split("/")
    mesh, dirichlet, stiff = cube(int(value)) if family == "cube" else jump(float(value))
    fes = ng.H1(mesh, order=1, dirichlet=dirichlet)
    coefficients = route.Coefficients(stiff=stiff, stiff_aux=stiff)
    v = fes.TestFunction()
    form = ng.LinearForm(fes)
    form += (1 + ng.x * ng.y) * v * ng.dx
    form.Assemble()
    rhs = np.array(form.vec.FV().NumPy(), copy=True)
    lift_function = ng.GridFunction(fes)
    lift_function.Set(1 + ng.x, definedon=mesh.Boundaries(dirichlet))
    lift = np.array(lift_function.vec.FV().NumPy(), copy=True)
    started = time.perf_counter()
    system = route.DeviceSystem(fes, coefficients)
    system.setup("amg", min_coarse=min_coarse)
    device = system.solve(rhs, lift=lift, tolerance=TOLERANCE, max_iterations=1000)
    device_seconds = time.perf_counter() - started
    with ng.TaskManager():
        direct = route.cpu_solve(fes, coefficients, rhs, lift=lift, preconditioner="direct")
    cpu = route.cpu_solve(fes, coefficients, rhs, lift=lift, preconditioner="amg", tolerance=TOLERANCE, max_iterations=1000, min_coarse=min_coarse)
    difference = float(np.linalg.norm(device["solution"] - direct["solution"]) / np.linalg.norm(direct["solution"]))
    return {"case": text, "ndof": int(fes.ndof), "elements": int(mesh.ne), "solution_difference_vs_sparsecholesky": difference,
            "device_iterations": device["iterations"], "cpu_iterations": cpu["iterations"], "device_true_residual": device["true_residual"],
            "cpu_true_residual": cpu["true_residual"], "direct_true_residual": direct["true_residual"], "hierarchy": device["hierarchy"],
            "nonzeros": device["nonzeros"], "device": device["device"], "device_peak_bytes": device["device_peak_bytes"], "ledger": device["ledger"],
            "timing_diagnostic": {**device["timing_diagnostic"], "route_total_seconds": device_seconds, "cpu_direct_seconds": direct["timing_diagnostic"]["total_seconds"],
                                  "cpu_amg_seconds": cpu["timing_diagnostic"]["total_seconds"]},
            "detector_status": "PASS" if difference <= SOLUTION_REL else "FAIL"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cases", nargs="*", default=["cube/8", "cube/16", "cube/32", "cube/48", "jump/0.2", "jump/0.1", "jump/0.05"])
    parser.add_argument("--min-coarse", type=int, default=route.AMG_MIN_COARSE)
    args = parser.parse_args()
    cases = []
    for text in args.cases:
        case = run_case(text, args.min_coarse)
        cases.append(case)
        print(f"{text}: ndof {case['ndof']} device {case['device_iterations']} CPU {case['cpu_iterations']} "
              f"difference {case['solution_difference_vs_sparsecholesky']:.3e} true residual {case['device_true_residual']:.3e}", flush=True)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    result = {"schema": "radia.fem-device-amg-ladder.v1", "source_commit": head, "ngsolve": ng.__version__, "python": platform.python_version(),
              "conditions": {"tolerance": TOLERANCE, "theta": route.AMG_THETA, "min_coarse": args.min_coarse, "solution_rel_proposed": SOLUTION_REL},
              "cases": cases, "decision": "PASS" if all(c["detector_status"] == "PASS" for c in cases) else "FAIL"}
    args.out.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"decision": result["decision"]}))
    raise SystemExit(0 if result["decision"] == "PASS" else 1)


if __name__ == "__main__":
    main()
