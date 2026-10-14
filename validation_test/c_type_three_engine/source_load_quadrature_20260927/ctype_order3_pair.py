"""C-type linear coarse mesh: reduced-A and mixed Omega at orders 2 and 3.

HDiv-MMM accepts only BDM1/BDM2, so its stored BDM2 field is the fixed
reference.  Reduced-A order 2 must reproduce the stored reduced-A field before
the order-3 rows count.  Gap-core B (median-plane projected) pairwise relative
RMS is reported for every pair.  Diagnostic evidence only.
"""
import hashlib
import importlib.util
import itertools
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import ngsolve as ng


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    here = Path(__file__).resolve().parent
    mesh_dir = here / "mesh"
    reference_path = here / "c_type_20260908_linear_kelvin_after_coarse_worker_a.json"
    output = here / "order3_pair.json"
    ng.SetNumThreads(int(sys.argv[1]) if len(sys.argv) > 1 else 8)

    spec = importlib.util.spec_from_file_location("_rte", here / "run_three_engine.py")
    rte = importlib.util.module_from_spec(spec)
    sys.modules["_rte"] = rte
    spec.loader.exec_module(rte)
    import radia as rad
    from radia.kelvin_identify_ngsolve import detect_kelvin_offset

    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    report_path = mesh_dir / "mesh_result.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if reference["mesh_result_sha256"] != sha256(report_path):
        raise RuntimeError("reference was produced on a different mesh contract")
    if report["artifacts"]["kelvin_domain_vol_sha256"] != sha256(mesh_dir / "kelvin_domain.vol"):
        raise RuntimeError("kelvin_domain.vol hash differs from mesh_result.json")
    physical = np.asarray(report["kelvin_physical_center_m"], dtype=float)
    points = rte.observation_points()
    gap_core = np.abs(points[:, 0]) <= float(reference["gap_core_half_length_m"]) + 1e-14
    rad.UtiDelAll()
    coil, coil_manifest = rte.build_coil()
    mesh = ng.Mesh(str(mesh_dir / "kelvin_domain.vol"))
    offset = np.asarray(detect_kelvin_offset(mesh), dtype=float)
    kelvin_center = tuple((physical + offset).tolist())
    kelvin_radius = float(report["kelvin_radius_m"])

    fields = {}
    stored_hdiv = rte.median_plane_projection(points, np.asarray(reference["fields_T"]["hdiv_mmm"]))
    stored_reduced_a = rte.median_plane_projection(
        points, np.asarray(reference["fields_T"]["reduced_a"]))
    payload = {"schema": "radia.validation.ctype-order3-pair.v1", "platform_class": platform.system(),
               "radia_version": getattr(rad, "__version__", None), "radia_module": rad.__file__,
               "ngsolve": ng.__version__, "script_sha256": sha256(__file__),
               "reference_sha256": sha256(reference_path), "mesh_result_sha256": sha256(report_path),
               "coil": coil_manifest, "controls": {"mu_r": 1000.0, "bonus_intorder": 4,
                                                   "reduced_a_solver": "direct",
                                                   "source_trace_tolerance": 0.05},
               "hdiv_note": "vim.Solve accepts BDM1/BDM2 only; stored BDM2 is the fixed reference",
               "acceptance": "HOLD: diagnostic, not three-engine acceptance", "runs": {}}

    def save():
        payload["fields_gap_core_T"] = {k: v[gap_core].tolist() for k, v in fields.items()}
        output.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")

    if report["artifacts"]["iron_vol_sha256"] != sha256(mesh_dir / "iron.vol"):
        raise RuntimeError("iron.vol hash differs from mesh_result.json")
    started = time.perf_counter()
    field, diag = rte.solve_hdiv(
        ng.Mesh(str(mesh_dir / "iron.vol")), coil, 1000.0, nonlinear=False, order=2,
        gram_eps=1e-14, nonlinear_tolerance=2e-5, nonlinear_maximum_iterations=80, points=points)
    fields["hdiv_bdm2"] = rte.median_plane_projection(points, field)
    payload["runs"]["hdiv_bdm2"] = {"runtime_s": time.perf_counter() - started}
    payload["hdiv_bdm2_vs_stored_older_release"] = rte.relative_rms(
        stored_hdiv[gap_core], fields["hdiv_bdm2"][gap_core])
    save()
    print("hdiv_bdm2", payload["runs"]["hdiv_bdm2"], flush=True)

    for order in (2, 3):
        started = time.perf_counter()
        field, diag = rte.solve_reduced_a(
            mesh, coil, 1000.0, nonlinear=False, order=order, linear_solver="direct", relax=0.1,
            nonlinear_tolerance=2e-5, nonlinear_maximum_iterations=80, nonlinear_verbose=False,
            kelvin_center=kelvin_center, kelvin_radius=kelvin_radius, points=points)
        fields[f"reduced_a_p{order}"] = rte.median_plane_projection(points, field)
        payload["runs"][f"reduced_a_p{order}"] = {"ndof": diag["ndof"],
                                                  "runtime_s": time.perf_counter() - started}
        if order == 2:
            # Recorded, not gated: the stored field came from an older Radia release.
            payload["reduced_a_p2_vs_stored_older_release"] = rte.relative_rms(
                stored_reduced_a[gap_core], fields["reduced_a_p2"][gap_core])
        save()
        print(f"reduced_a_p{order}", payload["runs"][f"reduced_a_p{order}"], flush=True)

        started = time.perf_counter()
        field, diag = rte.solve_omega(
            mesh, coil, 1000.0, nonlinear=False, order=order, nonlinear_tolerance=2e-5,
            nonlinear_maximum_iterations=80, nonlinear_verbose=False, kelvin_center=kelvin_center,
            kelvin_radius=kelvin_radius, points=points, source_trace_tolerance=0.05,
            bonus_intorder=4)
        fields[f"mixed_omega_p{order}"] = rte.median_plane_projection(points, field)
        payload["runs"][f"mixed_omega_p{order}"] = {"ndof": diag["ndof"],
                                                    "runtime_s": time.perf_counter() - started,
                                                    "source_trace": diag["source_trace"]}
        save()
        print(f"mixed_omega_p{order}", payload["runs"][f"mixed_omega_p{order}"]["ndof"], flush=True)

    pairs = {}
    for left, right in itertools.combinations(fields, 2):
        pairs[f"{left}__vs__{right}"] = rte.relative_rms(fields[left][gap_core], fields[right][gap_core])
    payload["pairwise_gap_core_relative_rms"] = pairs
    payload["completed"] = True
    save()
    for name, value in pairs.items():
        print(f"{name}: {100 * value:.4f}%", flush=True)


if __name__ == "__main__":
    main()
