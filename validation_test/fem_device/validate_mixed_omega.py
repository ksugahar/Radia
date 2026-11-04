"""Device-route matching-trace mixed Omega against the CPU route on identical inputs.

The CPU route ``solve_magnetostatic_matching_trace_total_reduced_omega`` assembles the source load,
the interface lift and the boundary values; its sparse Cholesky solution is the field reference and
its BDDC with AMG coarse grid gives the CPU iteration count. The same load and boundary values then
cross to the device, where the form and the preconditioner are built from mesh arrays: AMG at p=1
and BDDC at p=2,3.

Parity with MIXED_OMEGA_BACKENDS.md stays OPEN for the formulation as a whole: this case covers the
condensed finite-domain form only, not the Kelvin saddle-point driver or the nonlinear tangent.
SOLUTION_REL and FIELD_REL are PROPOSED detectors. Writes one JSON. Timings are diagnostic. Usage:
    python validate_mixed_omega.py --out result.json [--cases 0.3/1 0.3/2 0.3/3 ...]
"""
import argparse
import json
from pathlib import Path
import platform
import subprocess
import time

import numpy as np
import ngsolve as ng

import radia.fem_gpu as route
from radia.kelvin_material import MU_0
from radia.kelvin_solver import solve_magnetostatic_matching_trace_total_reduced_omega as solve_cpu

ROOT = Path(__file__).resolve().parents[2]
SOLUTION_REL = 1e-6
FIELD_REL = 1e-6
MU_R = {"reduced": 1.0, "total": 100.0}
BONUS = 4


def model(maxh):
    from netgen.occ import Box, Glue, OCCGeometry, Pnt, X

    reduced = Box(Pnt(-1, -1, -1), Pnt(0, 1, 1))
    reduced.mat("reduced")
    reduced.faces.name = "natural"
    reduced.faces.Min(X).name = "norm_boundary"
    reduced.faces.Max(X).name = "source_total_interface"
    total = Box(Pnt(0, -1, -1), Pnt(1, 1, 1))
    total.mat("total")
    total.faces.name = "natural"
    total.faces.Max(X).name = "norm_boundary"
    total.faces.Min(X).name = "source_total_interface"
    mesh = ng.Mesh(OCCGeometry(Glue([reduced, total])).GenerateMesh(maxh=maxh))
    source_h = ng.CoefficientFunction((ng.x * ng.x + 0.1 * ng.y, 0.2 * ng.y, -0.15 * ng.z))
    return mesh, source_h


def run_case(text, tolerance):
    maxh, order = float(text.split("/")[0]), int(text.split("/")[1])
    mesh, source_h = model(maxh)
    trace_space = ng.H1(mesh, order=order, definedon=mesh.Boundaries("source_total_interface"))
    trace = ng.GridFunction(trace_space)
    trace.Set(-ng.z, definedon=mesh.Boundaries("source_total_interface"))
    common = dict(mu_r_by_material=MU_R, reduced_materials=("reduced",), total_materials=("total",), interface_boundary="source_total_interface", order=order,
                  dirichlet_boundary="norm_boundary", bonus_intorder=BONUS)
    with ng.TaskManager():
        direct = solve_cpu(mesh, source_h, trace, solver="direct", return_system=True, **common)
    started = time.perf_counter()
    with ng.TaskManager():
        iterative = solve_cpu(mesh, source_h, trace, solver="cg", cg_preconditioner="bddc_h1amg", cg_tolerance=tolerance, **common)
    cpu_seconds = time.perf_counter() - started
    fes = iterative["fes"]
    reference = np.array(direct["solution"].vec.FV().NumPy(), copy=True)
    if fes.ndof != direct["fes"].ndof:
        raise RuntimeError(f"the direct and iterative spaces differ: {direct['fes'].ndof} and {fes.ndof} DOFs")
    rhs = np.array(direct["system"]["linear_form"].vec.FV().NumPy(), copy=True)
    free = np.array(list(fes.FreeDofs()), dtype=bool)
    boundary = reference.copy()
    boundary[free] = 0.0
    mu = {name: float(MU_0 * value) for name, value in MU_R.items()}
    coefficients = route.Coefficients(stiff=mu, stiff_aux=mu)
    started = time.perf_counter()
    device = route.solve_h1(fes, coefficients, rhs, Backend="gpu", lift=boundary, bonus=(BONUS, 0, 0), tolerance=tolerance, max_iterations=2000)
    device_seconds = time.perf_counter() - started
    entry = solve_cpu(mesh, source_h, trace, solver="cg", cg_preconditioner="bddc_h1amg", cg_tolerance=tolerance, Backend="gpu", **common)
    entry_difference = float(np.linalg.norm(entry["solution"].vec.FV().NumPy() - reference) / np.linalg.norm(reference))
    difference = float(np.linalg.norm(device["solution"] - reference) / np.linalg.norm(reference))
    g = ng.GridFunction(fes)
    g.vec.FV().NumPy()[:] = device["solution"] - reference
    h = ng.GridFunction(fes)
    h.vec.FV().NumPy()[:] = reference
    field = float(np.sqrt(ng.Integrate(ng.InnerProduct(ng.grad(g), ng.grad(g)), mesh, order=2 * order + 2) / ng.Integrate(ng.InnerProduct(ng.grad(h), ng.grad(h)), mesh, order=2 * order + 2)))
    return {"case": text, "order": order, "maxh": maxh, "ndof": int(fes.ndof), "elements": int(mesh.ne), "preconditioner": device["preconditioner"], "device_iterations": device["iterations"],
            "cpu_bddc_h1amg_iterations": iterative["iterations"], "device_true_residual": device["true_residual"], "direct_relative_residual": direct["linear_residual_relative"],
            "production_entry_difference_vs_sparsecholesky": entry_difference, "production_entry_iterations": entry["iterations"],
            "solution_difference_vs_sparsecholesky": difference, "grad_phi_difference_vs_sparsecholesky": field, "hierarchy": device["hierarchy"], "nonzeros": device["nonzeros"],
            "device": device["device"], "device_peak_bytes": device["device_peak_bytes"], "ledger": device["ledger"],
            "timing_diagnostic": {**device["timing_diagnostic"], "route_total_seconds": device_seconds, "cpu_cg_route_seconds": cpu_seconds},
            "detector_status": "PASS" if max(difference, entry_difference) <= SOLUTION_REL and field <= FIELD_REL else "FAIL"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cases", nargs="*", default=["0.5/1", "0.5/2", "0.5/3", "0.2/1", "0.2/2", "0.2/3", "0.1/2"])
    parser.add_argument("--tolerance", type=float, default=1e-10)
    args = parser.parse_args()
    cases = []
    for text in args.cases:
        case = run_case(text, args.tolerance)
        cases.append(case)
        print(f"{text}: ndof {case['ndof']} device {case['device_iterations']} CPU {case['cpu_bddc_h1amg_iterations']} difference {case['solution_difference_vs_sparsecholesky']:.3e} "
              f"grad difference {case['grad_phi_difference_vs_sparsecholesky']:.3e} true residual {case['device_true_residual']:.3e}", flush=True)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    result = {"schema": "radia.fem-device-mixed-omega.v1", "source_commit": head, "ngsolve": ng.__version__, "python": platform.python_version(),
              "conditions": {"mu_r": MU_R, "bonus_intorder": BONUS, "tolerance": args.tolerance, "solution_rel_proposed": SOLUTION_REL, "field_rel_proposed": FIELD_REL},
              "parity_status": "OPEN for the formulation; this covers the condensed finite-domain form only", "cases": cases,
              "decision": "PASS" if all(c["detector_status"] == "PASS" for c in cases) else "FAIL"}
    args.out.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"decision": result["decision"]}))
    raise SystemExit(0 if result["decision"] == "PASS" else 1)


if __name__ == "__main__":
    main()
