"""Device-route complex curl-curl with surface impedance on the bored workpiece air mesh.

The mesh, the coil and the CPU solve are those of validation_test/induction_heating/bored_workpiece.
The production CPU route ``solve_fem`` with sparse Cholesky supplies the reference field and the
line-integral right-hand side. ``solve_fem(..., Backend="gpu")`` then builds the form and the
preconditioner on the device from mesh arrays: AMS at p=1 and BDDC over AMS at p=2.

The bored_workpiece README gates FEM against its BEM power at 2 %, and its FEM power balance at
1e-5. Those values are read from validate.LIMITS, unchanged; no approval ID for them was found, so
this script reports them as the case's gates, not as an adopted criterion. That case evaluates its
balance at p=2, so it is gated at p=2 here; at p=1 sparse Cholesky itself gives 4e-5, so the balance is
recorded and the source work is compared with sparse Cholesky. FIELD_REL and POWER_REL
are PROPOSED detectors. The air carries only the small gauge mass, so the gradient part of A is weakly
determined and the raw vector difference is recorded, not gated: the comparison uses the curl in the
air and the tangential trace on the surface. Writes one JSON. Timings are diagnostic. Usage:
    python validate_sibc_bored.py --reference-only --out unused --work scratch_dir [--cases 0.012/1 0.012/2 ...]
    python validate_sibc_bored.py --out result.json --work scratch_dir [--cases 0.012/1 0.012/2 ...]
The first form is the CPU stage, with no device; the second runs the production entry on the device and
compares it with the stored reference on the identical right-hand side.
"""
import argparse
import importlib.util
import json
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

FIELD_REL = 1e-6
POWER_REL = 1e-6
BEM_FINEST_W = 73.634332  # bored_workpiece README, finest BEM surface
AIR_RADIUS = .18


def load_case_module():
    spec = importlib.util.spec_from_file_location("bored_workpiece_validate", ROOT / "validation_test/induction_heating/bored_workpiece/validate.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_mesh(case, h, work, air_radius=AIR_RADIUS):
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
    from radia.vol_check import check_consistency
    contract = {"schema": "radia.vol-label-contract.v1", "application": "bored_workpiece", "required": {"materials": ["air"], "boundaries": ["sibc", "outer"]}}
    checked = check_consistency(str(vol), contract=contract)
    if not checked["passed"]:
        raise RuntimeError("air mesh failed check-vol: " + str(checked["warnings"]))
    return vol


def capture_solve(case, kelvin, vol, order, **options):
    """The production ``solve_fem`` with its assembled right-hand side and field captured at return."""
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
                                      half_thickness=.005, nthreads=0, filaments=(paths, currents), **options)
    finally:
        sys.settrace(previous)
    if "error" in result:
        raise RuntimeError(result["error"])
    if not captured:
        raise RuntimeError("the production route returned without its right-hand side being captured")
    return result, captured, time.perf_counter() - started


def parse_case(text):
    return float(text.split("/")[0]), int(text.split("/")[1])


def reference_path(work, text):
    h, order = parse_case(text)
    return work / f"reference_{h:g}_{order}.npz"


def run_reference(case, kelvin, text, work):
    """The CPU stage: sparse Cholesky on the checked mesh. It holds no device and may run outside the device slot."""
    h, order = parse_case(text)
    vol = build_mesh(case, h, work)
    result, data, seconds = capture_solve(case, kelvin, vol, order, solver="sparsecholesky")
    np.savez(reference_path(work, text), rhs=data["rhs"], field=data["field"], P_total=float(result["P_total"]), seconds=seconds, ndof=int(data["fes"].ndof),
             true_residual=float(result["linear_true_relative_residual"]))
    print(f"reference {text}: ndof {data['fes'].ndof} power {result['P_total']:.6f} W, {seconds:.1f} s", flush=True)


def surface_power(g, mesh, sibc, omega, Z_s):
    """P = Re(Z_s)/2 * (omega/|Z_s|)^2 * int |A_t|^2 dS, the expression the CPU route evaluates."""
    normal = ng.specialcf.normal(3)
    squared = sum(g[i].real * g[i].real + g[i].imag * g[i].imag for i in range(3))
    normal_re = sum(g[i].real * normal[i] for i in range(3))
    normal_im = sum(g[i].imag * normal[i] for i in range(3))
    tangential = ng.Integrate(squared - normal_re * normal_re - normal_im * normal_im, mesh, ng.BND, definedon=mesh.Boundaries(sibc)).real
    return 0.5 * Z_s.real * (omega / abs(Z_s)) ** 2 * tangential


def gauge_invariant_differences(fes, mesh, sibc, actual, reference):
    """norm(curl(A1 - A2)) over norm(curl A2) in the air, and the same ratio of the tangential trace on the surface."""
    g, r = ng.GridFunction(fes), ng.GridFunction(fes)
    g.vec.FV().NumPy()[:] = actual
    r.vec.FV().NumPy()[:] = reference
    order = 2 * fes.globalorder + 2

    def volume(cf):
        return ng.Integrate(ng.InnerProduct(cf, ng.Conj(cf)), mesh, order=order).real

    def surface(cf):
        tangent = cf - ng.InnerProduct(cf, ng.specialcf.normal(3)) * ng.specialcf.normal(3)
        return ng.Integrate(ng.InnerProduct(tangent, ng.Conj(tangent)), mesh, ng.BND, definedon=mesh.Boundaries(sibc), order=order).real

    field = float(np.sqrt(volume(ng.curl(g) - ng.curl(r)) / volume(ng.curl(r))))
    trace = float(np.sqrt(surface(g - r) / surface(r)))
    return field, trace


def run_case(case, kelvin, text, work):
    h, order = parse_case(text)
    reference_file = np.load(reference_path(work, text))
    vol = work / f"air_{h:g}_{AIR_RADIUS:g}.vol"
    solver = "ams" if order == 1 else "bddc"
    entry, data, device_seconds = capture_solve(case, kelvin, vol, order, solver=solver, Backend="gpu")
    if entry["backend"] != "gpu" or not entry["device_report"]:
        raise RuntimeError(f"the production entry did not report the device route: {entry['backend']}")
    device = entry["device_report"]
    fes, mesh = data["fes"], data["mesh"]
    omega, Z_s, sibc = data["omega"], data["Z_s"], data["sibc"]
    rhs, reference = reference_file["rhs"], reference_file["field"]
    if rhs.shape != data["rhs"].shape or np.linalg.norm(rhs - data["rhs"]) > 1e-12 * np.linalg.norm(rhs):
        raise RuntimeError("the device and CPU right-hand sides differ, so the inputs are not identical")
    free = np.array(list(fes.FreeDofs()), dtype=bool)
    difference = float(np.linalg.norm(data["field"][free] - reference[free]) / np.linalg.norm(reference[free]))
    field_difference, trace_difference = gauge_invariant_differences(fes, mesh, sibc, data["field"], reference)
    g = ng.GridFunction(fes)
    g.vec.FV().NumPy()[:] = data["field"]
    power = surface_power(g, mesh, sibc, omega, Z_s)
    g_cpu = ng.GridFunction(fes)
    g_cpu.vec.FV().NumPy()[:] = reference
    power_cpu_same_formula = surface_power(g_cpu, mesh, sibc, omega, Z_s)
    source_work = float(.5 * omega * np.vdot(data["field"], rhs).imag)
    source_work_cpu = float(.5 * omega * np.vdot(reference, rhs).imag)
    power_entry, power_cpu = float(entry["P_total"]), float(reference_file["P_total"])
    limits = case.LIMITS
    gates = {"true_residual": {"actual": device["true_residual"], "bound": limits["true_residual"], "pass": device["true_residual"] <= limits["true_residual"]},
             "bem_fem_power": {"actual": abs(power_entry - BEM_FINEST_W) / BEM_FINEST_W, "bound": limits["bem_fem_power"], "pass": abs(power_entry - BEM_FINEST_W) / BEM_FINEST_W <= limits["bem_fem_power"]},
             "curl_field_vs_sparsecholesky_proposed": {"actual": field_difference, "bound": FIELD_REL, "pass": field_difference <= FIELD_REL},
             "surface_trace_vs_sparsecholesky_proposed": {"actual": trace_difference, "bound": FIELD_REL, "pass": trace_difference <= FIELD_REL},
             "power_vs_cpu_route_proposed": {"actual": abs(power_entry - power_cpu) / power_cpu, "bound": POWER_REL, "pass": abs(power_entry - power_cpu) / power_cpu <= POWER_REL},
             "source_work_vs_sparsecholesky_proposed": {"actual": abs(source_work / source_work_cpu - 1), "bound": POWER_REL, "pass": abs(source_work / source_work_cpu - 1) <= POWER_REL}}
    if order == 2:
        gates["fem_power_balance"] = {"actual": abs(source_work / power - 1), "bound": limits["fem_power_balance"], "pass": abs(source_work / power - 1) <= limits["fem_power_balance"]}
    return {"case": text, "maxh_sibc_m": h, "order": order, "ndof": int(fes.ndof), "elements": int(mesh.ne), "device_iterations": device["iterations"], "true_residual": device["true_residual"],
            "preconditioner": device["preconditioner"], "krylov": "COCR", "coarse_cycles": kelvin.BDDC_AMS_COARSE_CYCLES if order > 1 else None, "raw_vector_difference_vs_sparsecholesky_recorded": difference,
            "curl_field_difference_vs_sparsecholesky": field_difference, "surface_trace_difference_vs_sparsecholesky": trace_difference,
            "power_W": power_entry, "power_field_formula_W": power, "power_cpu_same_formula_W": power_cpu_same_formula, "power_cpu_route_W": power_cpu, "source_work_W": source_work, "source_work_cpu_W": source_work_cpu, "fem_power_balance_recorded": abs(source_work / power - 1), "fem_power_balance_sparsecholesky_recorded": abs(source_work_cpu / power_cpu_same_formula - 1), "bem_finest_W": BEM_FINEST_W,
            "gates": gates, "hierarchy": device["hierarchy"], "nonzeros": device["nonzeros"], "device": device["device"], "device_peak_bytes": device["device_peak_bytes"], "ledger": device["ledger"],
            "timing_diagnostic": {**device["timing_diagnostic"], "production_entry_seconds": device_seconds, "cpu_route_seconds": float(reference_file["seconds"])}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--cases", nargs="*", default=["0.012/1", "0.012/2"])
    parser.add_argument("--reference-only", action="store_true", help="the CPU stage: write the sparse Cholesky reference into --work")
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    case = load_case_module()
    import calc_fem_kelvin as kelvin
    if args.reference_only:
        for text in args.cases:
            run_reference(case, kelvin, text, args.work)
        return
    cases = []
    for text in args.cases:
        row = run_case(case, kelvin, text, args.work)
        cases.append(row)
        print(f"{text}: ndof {row['ndof']} iterations {row['device_iterations']} true residual {row['true_residual']:.3e} curl {row['curl_field_difference_vs_sparsecholesky']:.3e} trace {row['surface_trace_difference_vs_sparsecholesky']:.3e} raw {row['raw_vector_difference_vs_sparsecholesky_recorded']:.3e} "
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
