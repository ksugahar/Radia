"""Device-route replay of the linear volumetric F1 slab against its approved criteria.

The model, the quantity definitions and the gate rules are those of adoption record
validation_records/induction_heating/20261024_ams_linear_slab; every threshold is read from its
criterion.json, unchanged. The six cases run with device assembly and a device preconditioner:
complex AMS at p=1, BDDC over AMS at p=2,3. The reference is the independent analytic diffusion
solution, and the same-matrix sparse Cholesky field of the CPU form for parity.

Writes one JSON with the source commit, conditions, values, gates and the transfer ledger.
Timings are diagnostic and never gate. Usage:
    python validate_f1_slab.py --out result.json [--cases 2/1 2/2 ...]
"""
import argparse
import cmath
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh

import radia.fem_gpu as route

ROOT = Path(__file__).resolve().parents[2]
RECORD = ROOT / "validation_records/induction_heating/20261024_ams_linear_slab"
MU, SIGMA, DELTA, H0 = 4 * math.pi * 1e-7, 5e6, 1e-3, 1000.0
OMEGA = 2 / (MU * SIGMA * DELTA ** 2)
HALF, WIDTH, K = 2 * DELTA, DELTA / 4, (1 + 1j) / DELTA
Z_REFERENCE = (K / SIGMA) * cmath.tanh(K * HALF)
POWER_REFERENCE = WIDTH * WIDTH * H0 ** 2 * Z_REFERENCE.real


def analytic():
    Ay = MU * H0 * (ng.exp(K * ng.x) - ng.exp(-K * ng.x)) / (2 * K * cmath.cosh(K * HALF))
    Hz = H0 * (ng.exp(K * ng.x) + ng.exp(-K * ng.x)) / (2 * cmath.cosh(K * HALF))
    return ng.CF((0, Ay, 0)), ng.CF((0, 0, Hz))


def l2(field, mesh, order):
    return math.sqrt(ng.Integrate(ng.InnerProduct(field, field).real, mesh, order=order))


def face_impedance(g, mesh, name, sign, order):
    region = mesh.Boundaries(name)
    electric = ng.BoundaryFromVolumeCF(-1j * OMEGA * g[1])
    magnetic = ng.BoundaryFromVolumeCF(ng.curl(g)[2] / MU)
    return ng.Integrate(-sign * electric, mesh, ng.BND, definedon=region, order=order) / ng.Integrate(magnetic, mesh, ng.BND, definedon=region, order=order)


def quantities(g, mesh, p):
    A_exact, H_exact = analytic()
    order = 2 * p + 12
    J_exact = -1j * OMEGA * SIGMA * A_exact
    H, J = ng.curl(g) / MU, -1j * OMEGA * SIGMA * g
    errors = {"H": l2(H - H_exact, mesh, order) / l2(H_exact, mesh, order), "J": l2(J - J_exact, mesh, order) / l2(J_exact, mesh, order)}
    errors["B"] = errors["H"]
    z_back, z_front = face_impedance(g, mesh, "back", -1, order), face_impedance(g, mesh, "front", 1, order)
    errors["Z_back"] = abs(z_back - Z_REFERENCE) / abs(Z_REFERENCE)
    errors["Z_front"] = abs(z_front - Z_REFERENCE) / abs(Z_REFERENCE)
    power = float(ng.Integrate(ng.InnerProduct(J, J).real, mesh, order=order)) / (2 * SIGMA)
    errors["power"] = abs(power - POWER_REFERENCE) / POWER_REFERENCE
    return errors, power, {"back": [z_back.real, z_back.imag], "front": [z_front.real, z_front.imag]}


def hmax_over_half(mesh):
    longest = 0.0
    for element in mesh.Elements(ng.VOL):
        points = np.array([mesh[v].point for v in element.vertices])
        for i in range(len(points)):
            for j in range(i + 1, len(points)):
                longest = max(longest, float(np.linalg.norm(points[i] - points[j])))
    return longest / HALF


def run_case(m, p, cycles, tolerance):
    mesh = MakeStructured3DMesh(hexes=False, nx=16 * m, ny=m, nz=m, secondorder=False,
                                mapping=lambda x, y, z: (4 * DELTA * x - HALF, WIDTH * (y - .5), WIDTH * (z - .5)))
    space = (ng.HCurl(mesh, order=1, nograds=True, complex=True, dirichlet=".*") if p == 1
             else ng.HCurl(mesh, order=p, complex=True, dirichlet=".*"))
    coefficients = route.Coefficients(stiff=1 / MU, mass_im=OMEGA * SIGMA, stiff_aux=1 / MU, mass_aux=OMEGA * SIGMA)
    bonus = (2 * p + 2, 2 * p + 2, 0)
    A_exact, _ = analytic()
    g = ng.GridFunction(space)
    g.Set(A_exact, definedon=mesh.Boundaries(".*"))
    lift = np.array(g.vec.FV().NumPy(), copy=True)
    zeros = np.zeros(space.ndof, dtype=complex)
    started = time.perf_counter()
    device = route.solve_hcurl(space, coefficients, zeros, Backend="gpu", lift=lift, bonus=bonus, tolerance=tolerance, max_iterations=2000, cycles=cycles)
    device_seconds = time.perf_counter() - started
    free = np.array(list(space.FreeDofs()), dtype=bool)
    retained = float(np.abs(device["solution"][~free] - lift[~free]).max())
    g.vec.FV().NumPy()[:] = device["solution"]
    errors, power, faces = quantities(g, mesh, p)
    with ng.TaskManager():
        direct = route.cpu_solve(space, coefficients, zeros, lift=lift, preconditioner="direct", bonus=bonus)
    h = ng.GridFunction(space)
    h.vec.FV().NumPy()[:] = direct["solution"]
    order = 2 * p + 12
    parity = {"A": l2(g - h, mesh, order) / l2(h, mesh, order),
              "H": l2(ng.curl(g) - ng.curl(h), mesh, order) / l2(ng.curl(h), mesh, order),
              "J": l2(g - h, mesh, order) / l2(h, mesh, order)}
    return {"m": m, "p": p, "ndof": int(space.ndof), "tetrahedra": int(mesh.ne), "hmax_over_a": hmax_over_half(mesh), "errors": errors,
            "parity": parity, "true_residual": device["true_residual"], "direct_true_residual": direct["true_residual"], "Z_faces": faces,
            "power_W": power, "iterations": device["iterations"], "preconditioner": device["preconditioner"], "krylov": "COCR", "tolerance": tolerance,
            "coarse_cycles": cycles if p > 1 else None, "dirichlet_max_abs_difference": retained, "hierarchy": device["hierarchy"],
            "nonzeros": device["nonzeros"], "ledger": device["ledger"], "device": device["device"], "device_peak_bytes": device["device_peak_bytes"],
            "timing_diagnostic": {**device["timing_diagnostic"], "route_total_seconds": device_seconds, "cpu_direct_seconds": direct["timing_diagnostic"]["total_seconds"]}}


def gates(cases, criterion):
    finest = criterion["convergence"]["finest_relative_error"]
    ratio_gate = criterion["convergence"]["finest_adjacent_p_ratio"]
    solver = criterion["solver_gates"]
    results = []

    def gate(name, actual, bound, passed):
        results.append({"name": name, "actual": actual, "bound": bound, "status": "PASS" if passed else "FAIL"})

    def plateau(a, b):
        return min(a, b) <= 0 or max(a, b) <= 64 * sys.float_info.epsilon or abs(a - b) <= 16 * max(math.ulp(a), math.ulp(b))

    table = {(c["m"], c["p"]): c for c in cases}
    for c in cases:
        gate(f"{c['m']}/{c['p']}/true_residual", c["true_residual"], solver["iterative_and_direct_true_residual_upper"],
             c["true_residual"] <= solver["iterative_and_direct_true_residual_upper"])
        gate(f"{c['m']}/{c['p']}/direct_true_residual", c["direct_true_residual"], solver["iterative_and_direct_true_residual_upper"],
             c["direct_true_residual"] <= solver["iterative_and_direct_true_residual_upper"])
        for key, value in c["parity"].items():
            gate(f"parity/{c['m']}/{c['p']}/{key}", value, solver["same_matrix_physical_L2_parity_A_H_J_upper"],
                 value <= solver["same_matrix_physical_L2_parity_A_H_J_upper"])
    if not all((m, p) in table for m in (2, 4) for p in (1, 2, 3)):
        results.append({"name": "coverage", "status": "INCOMPLETE", "cases": sorted(table)})
        return results
    for key in table[(2, 1)]["errors"]:
        for p in (1, 2, 3):
            e = table[(4, p)]["errors"][key]
            if p in finest["record_only_p"]:
                results.append({"name": f"finest/{p}/{key}", "actual": e, "status": "RECORDED_ONLY"})
            else:
                gate(f"finest/{p}/{key}", e, finest["upper_bound"], e <= finest["upper_bound"])
            a, b = table[(2, p)], table[(4, p)]
            ei, ej = a["errors"][key], b["errors"][key]
            ratio_h = a["hmax_over_a"] / b["hmax_over_a"]
            if plateau(ei, ej) or ratio_h <= 1:
                results.append({"name": f"h/2/4/{p}/{key}", "actual_errors": [ei, ej], "status": "UNRESOLVED"})
                continue
            rate = math.log(ei / ej) / math.log(ratio_h)
            gate(f"h/2/4/{p}/{key}", rate, p - .5, ej < ei and rate >= p - .5)
        for m in (2, 4):
            for p in (1, 2):
                ei, ej = table[(m, p)]["errors"][key], table[(m, p + 1)]["errors"][key]
                if plateau(ei, ej):
                    results.append({"name": f"p_monotone/{m}/{p}/{key}", "actual_errors": [ei, ej], "status": "UNRESOLVED"})
                    continue
                gate(f"p_monotone/{m}/{p}/{key}", [ei, ej], "strict decreasing", 0 < ej < ei)
                if m == 4:
                    ratio = ej / ei
                    if key in ratio_gate["gate_quantities"]:
                        gate(f"p_ratio/{p}/{key}", ratio, ratio_gate["upper_bound"], ratio <= ratio_gate["upper_bound"])
                    else:
                        results.append({"name": f"p_ratio/{p}/{key}", "actual": ratio, "status": "RECORDED_ONLY"})
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cases", nargs="*", default=["2/1", "2/2", "2/3", "4/1", "4/2", "4/3"])
    parser.add_argument("--cycles", type=int, default=4)
    parser.add_argument("--tolerance", type=float, default=1e-11)
    args = parser.parse_args()
    criterion = json.loads((RECORD / "criterion.json").read_text(encoding="utf-8"))
    recorded = {(c["m"], c["p"]): c for c in json.loads((RECORD / "result.json").read_text(encoding="utf-8"))["cases"]}
    cases = []
    for text in args.cases:
        m, p = (int(v) for v in text.split("/"))
        case = run_case(m, p, args.cycles, args.tolerance)
        reference = recorded[(m, p)]
        case["recorded_cpu"] = {"iterations": reference["iterations"], "power_W": reference["power_W"], "errors": reference["errors"]}
        cases.append(case)
        print(f"m{m} p{p}: ndof {case['ndof']} iterations {case['iterations']} (CPU recorded {reference['iterations']}) "
              f"true residual {case['true_residual']:.3e} parity A {case['parity']['A']:.3e} H error {case['errors']['H']:.6e}", flush=True)
    evaluated = gates(cases, criterion)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    summary = {status: sum(1 for g in evaluated if g["status"] == status) for status in {g["status"] for g in evaluated}}
    result = {"schema": "radia.fem-device-f1-slab.v1", "source_commit": head, "record": "induction_heating/20261024_ams_linear_slab",
              "criterion_approvals": [item["id"] for item in criterion["approval_chain"]], "ngsolve": ng.__version__, "python": platform.python_version(),
              "conditions": {"mu": MU, "sigma": SIGMA, "delta": DELTA, "H0": H0, "omega": OMEGA, "half_depth": HALF, "width": WIDTH,
                             "Z_reference": [Z_REFERENCE.real, Z_REFERENCE.imag], "P_reference_W": POWER_REFERENCE},
              "cases": cases, "gates": evaluated, "summary": summary,
              "decision": "PASS" if all(g["status"] in ("PASS", "RECORDED_ONLY") for g in evaluated) else "FAIL_OR_UNRESOLVED"}
    args.out.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"decision": result["decision"], "summary": summary}))
    raise SystemExit(0 if result["decision"] == "PASS" else 1)


if __name__ == "__main__":
    main()
