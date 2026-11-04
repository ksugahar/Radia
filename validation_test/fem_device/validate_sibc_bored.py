"""Device-route complex curl-curl with surface impedance on the bored workpiece air mesh.

The mesh, the coil and the CPU solve are those of validation_test/induction_heating/bored_workpiece.
The production CPU route ``solve_fem`` with sparse Cholesky supplies the reference field and the
line-integral right-hand side. The same right-hand side then crosses to the device, where the form
and the preconditioner are built from mesh arrays: AMS at p=1 and BDDC over AMS at p=2.

The bored_workpiece README gates FEM against its BEM power at 2 %, and its FEM power balance at
1e-5. Those values are read from validate.LIMITS, unchanged; no approval ID for them was found, so
this script reports them as the case's gates, not as an adopted criterion. SOLUTION_REL is a
PROPOSED detector. Writes one JSON. Timings are diagnostic. Usage:
    python validate_sibc_bored.py --out result.json --work scratch_dir [--cases 0.006/1 0.006/2 ...]
"""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import ngsolve as ng

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src/radia/panels"))
sys.path.insert(0, str(ROOT / "src"))
import radia.fem_gpu as route  # noqa: E402

SOLUTION_REL = 1e-6
BEM_FINEST_W = 73.634332  # bored_workpiece README, finest BEM surface
NU_0 = 1 / (4e-7 * math.pi)


def load_case_module():
    spec = importlib.util.spec_from_file_location("bored_workpiece_validate", ROOT / "validation_test/induction_heating/bored_workpiece/validate.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_mesh(case, h, work, air_radius=.18):
    from netgen.occ import Cylinder, Glue, OCCGeometry, Pnt, Sphere, Z

    with ng.TaskManager():
        part = case.part()
        part.faces.name = "sibc"
        part.faces.maxh = h
        ring = Cylinder(Pnt(0, 0, -.06), Z, .065, .12)
        near = ring - part
        near.mat("air")
        near.maxh = .012
        sphere = Sphere(Pnt(0, 0, 0), air_radius)
        sphere.faces.name = "outer"
        far = sphere - ring
        far.mat("air")
        mesh = ng.Mesh(OCCGeometry(Glue([near, far])).GenerateMesh(maxh=.045, curvaturesafety=1.5))
        mesh.Curve(2)
        vol = work / f"air_{h:g}_{air_radius:g}.vol"
        mesh.ngmesh.Save(str(vol))
    from cubit_mesh_export.check import check_consistency
    contract = {"schema": "radia.vol-label-contract.v1", "application": "bored_workpiece", "required": {"materials": ["air"], "boundaries": ["sibc", "outer"]}}
    checked = check_consistency(str(vol), contract=contract)
    if not checked["passed"]:
        raise RuntimeError("air mesh failed check-vol: " + str(checked["warnings"]))
    return vol


def cpu_reference(case, kelvin, vol, order):
    """The production CPU solve with its assembled right-hand side and field captured at return."""
    from radia.em_material import EMMaterial

    captured = {}

    def trace(frame, event, arg):
        if frame.f_code is kelvin.solve_fem.__code__:
            frame.f_trace_lines = False

            def inner(inner_frame, inner_event, inner_arg):
                if inner_event == "return" and isinstance(inner_arg, dict) and "P_total" in inner_arg:
                    v = inner_frame.f_locals
                    captured.update(fes=v["fes"], mesh=v["mesh"], rhs=np.array(v["rhs_vec"].FV().NumPy(), copy=True), field=np.array(v["gfu"].vec.FV().NumPy(), copy=True),
                                    omega=float(v["omega"]), Z_s=complex(v["Z_s"]), reg=float(v["reg"]), shift_eps=float(v["shift_eps"]), sibc=v["sibc_bnd"])
                return inner

            return inner
        return None

    paths, currents = case.coils()
    previous = sys.gettrace()
    sys.settrace(trace)
    started = time.perf_counter()
    try:
        with ng.TaskManager():
            result = kelvin.solve_fem(vol_file=str(vol), fes_order=order, frequency=case.FREQ, mat=EMMaterial(name="linear alloy", sigma=case.SIGMA, mu_r=case.MUR),
                                      half_thickness=.005, solver="sparsecholesky", nthreads=0, filaments=(paths, currents))
    finally:
        sys.settrace(previous)
    if "error" in result:
        raise RuntimeError(result["error"])
    if not captured:
        raise RuntimeError("the CPU route returned without its right-hand side being captured")
    return result, captured, time.perf_counter() - started


def surface_power(g, mesh, sibc, Z_s):
    """P = Re(Z_s)/2 * int |n x H|^2 dS with H = curl A / mu_0, as the CPU route evaluates it."""
    normal = ng.specialcf.normal(3)
    H = NU_0 * ng.curl(g.Trace())
    tangential = ng.Cross(normal, H)
    value = ng.Integrate(ng.InnerProduct(tangential, tangential).real, mesh, ng.BND, definedon=mesh.Boundaries(sibc), order=2 * g.space.globalorder + 4)
    return 0.5 * Z_s.real * float(value)


def run_case(case, kelvin, text, work, cycles, tolerance):
    h, order = float(text.split("/")[0]), int(text.split("/")[1])
    vol = build_mesh(case, h, work)
    cpu, data, cpu_seconds = cpu_reference(case, kelvin, vol, order)
    fes, mesh = data["fes"], data["mesh"]
    omega, Z_s, sibc = data["omega"], data["Z_s"], data["sibc"]
    robin = 1j * omega / Z_s
    # The CPU route's system and auxiliary forms, as in solve_fem: curl-curl with bonus order 4, gauge
    # mass, Robin trace; the auxiliary mass adds |omega| NU_0.
    coefficients = route.Coefficients(stiff=NU_0, mass=data["reg"] * NU_0, trace={sibc: robin.real}, trace_im={sibc: robin.imag}, stiff_aux=NU_0,
                                      mass_aux=(data["shift_eps"] + abs(omega)) * NU_0, trace_aux={sibc: abs(robin)})
    bonus = (4, 0, 0)
    started = time.perf_counter()
    system = route.DeviceSystem(fes, coefficients, bonus=bonus)
    system.setup("ams" if order == 1 else "bddc", cycles=cycles)
    device = system.solve(data["rhs"], tolerance=tolerance, max_iterations=2000)
    device_seconds = time.perf_counter() - started
    free = np.array(list(fes.FreeDofs()), dtype=bool)
    reference = data["field"].copy()
    difference = float(np.linalg.norm(device["solution"][free] - reference[free]) / np.linalg.norm(reference[free]))
    g = ng.GridFunction(fes)
    g.vec.FV().NumPy()[:] = device["solution"]
    power = surface_power(g, mesh, sibc, Z_s)
    g_cpu = ng.GridFunction(fes)
    g_cpu.vec.FV().NumPy()[:] = reference
    power_cpu_same_formula = surface_power(g_cpu, mesh, sibc, Z_s)
    source_work = float(.5 * omega * np.vdot(device["solution"], data["rhs"]).imag)
    limits = case.LIMITS
    gates = {"true_residual": {"actual": device["true_residual"], "bound": limits["true_residual"], "pass": device["true_residual"] <= limits["true_residual"]},
             "fem_power_balance": {"actual": abs(source_work / power - 1), "bound": limits["fem_power_balance"], "pass": abs(source_work / power - 1) <= limits["fem_power_balance"]},
             "bem_fem_power": {"actual": abs(power - BEM_FINEST_W) / BEM_FINEST_W, "bound": limits["bem_fem_power"], "pass": abs(power - BEM_FINEST_W) / BEM_FINEST_W <= limits["bem_fem_power"]},
             "solution_vs_sparsecholesky_proposed": {"actual": difference, "bound": SOLUTION_REL, "pass": difference <= SOLUTION_REL}}
    return {"case": text, "maxh_sibc_m": h, "order": order, "ndof": int(fes.ndof), "elements": int(mesh.ne), "device_iterations": device["iterations"], "true_residual": device["true_residual"],
            "preconditioner": device["preconditioner"], "krylov": "COCR", "coarse_cycles": cycles if order > 1 else None, "solution_difference_vs_sparsecholesky": difference, "power_W": power,
            "power_cpu_same_formula_W": power_cpu_same_formula, "power_cpu_route_W": float(cpu["P_total"]), "source_work_W": source_work, "bem_finest_W": BEM_FINEST_W,
            "gates": gates, "hierarchy": device["hierarchy"], "nonzeros": device["nonzeros"], "device": device["device"], "device_peak_bytes": device["device_peak_bytes"], "ledger": device["ledger"],
            "timing_diagnostic": {**device["timing_diagnostic"], "route_total_seconds": device_seconds, "cpu_route_seconds": cpu_seconds}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--cases", nargs="*", default=["0.006/1", "0.006/2"])
    parser.add_argument("--cycles", type=int, default=4)
    parser.add_argument("--tolerance", type=float, default=1e-8)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    case = load_case_module()
    import calc_fem_kelvin as kelvin
    cases = []
    for text in args.cases:
        row = run_case(case, kelvin, text, args.work, args.cycles, args.tolerance)
        cases.append(row)
        print(f"{text}: ndof {row['ndof']} iterations {row['device_iterations']} true residual {row['true_residual']:.3e} difference {row['solution_difference_vs_sparsecholesky']:.3e} "
              f"power {row['power_W']:.6f} W (CPU {row['power_cpu_route_W']:.6f}, BEM {BEM_FINEST_W})", flush=True)
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    result = {"schema": "radia.fem-device-sibc-bored.v1", "source_commit": head, "ngsolve": ng.__version__, "python": platform.python_version(), "approval_id_for_gates": None,
              "note": "bored_workpiece limits are read unchanged; no approval ID was found for them", "cases": cases,
              "decision": "PASS" if all(all(g["pass"] for g in c["gates"].values()) for c in cases) else "FAIL"}
    args.out.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"decision": result["decision"]}))
    raise SystemExit(0 if result["decision"] == "PASS" else 1)


if __name__ == "__main__":
    main()
