"""C-type mixed Omega: projected Kelvin source trace vs exact pulled-back exterior source.

Same linear order-2 coarse-mesh solve as run_three_engine.solve_omega, once with
the default projected ``kelvin_int`` trace and once with
``exact_exterior_source=True`` (KelvinRadiaFieldStrength).  Gap-core B is
compared with the stored HDiv-MMM and reduced-A fields of the same mesh.
Diagnostic evidence only, not three-engine acceptance.
"""
import hashlib
import importlib.util
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
    output = here / "exact_exterior.json"
    ng.SetNumThreads(int(sys.argv[1]) if len(sys.argv) > 1 else 6)

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
    if not np.allclose(physical, 0.0, atol=1e-15):
        # run_three_engine passes (0, 0, 0) as the physical Kelvin centre.
        raise RuntimeError(f"physical Kelvin centre is {physical}, not the origin")

    points = rte.observation_points()
    gap_core = np.abs(points[:, 0]) <= float(reference["gap_core_half_length_m"]) + 1e-14
    stored = {name: rte.median_plane_projection(points, np.asarray(values))
              for name, values in reference["fields_T"].items()}
    rad.UtiDelAll()
    coil, coil_manifest = rte.build_coil()
    mesh = ng.Mesh(str(mesh_dir / "kelvin_domain.vol"))
    offset = np.asarray(detect_kelvin_offset(mesh), dtype=float)
    kelvin_center = tuple((physical + offset).tolist())
    kelvin_radius = float(report["kelvin_radius_m"])

    payload = {"schema": "radia.validation.ctype-exact-exterior-source.v1",
               "platform_class": platform.system(), "radia_version": getattr(rad, "__version__", None),
               "radia_module": rad.__file__, "ngsolve": ng.__version__,
               "script_sha256": sha256(__file__), "reference_sha256": sha256(reference_path),
               "mesh_result_sha256": sha256(report_path), "coil": coil_manifest,
               "controls": {"order": 2, "bonus_intorder": 4, "mu_r": 1000.0,
                            "source_trace_tolerance": 0.05},
               "acceptance": "HOLD: diagnostic, not three-engine acceptance", "runs": {}}
    fields = {}
    for label, exact in (("projected_kelvin_trace", False), ("exact_exterior_source", True)):
        started = time.perf_counter()
        field, diagnostics = rte.solve_omega(
            mesh, coil, 1000.0, nonlinear=False, order=2, nonlinear_tolerance=2e-5,
            nonlinear_maximum_iterations=80, nonlinear_verbose=False,
            kelvin_center=kelvin_center, kelvin_radius=kelvin_radius, points=points,
            source_trace_tolerance=0.05, bonus_intorder=4, exact_exterior_source=exact)
        projected = rte.median_plane_projection(points, field)
        fields[label] = projected
        row = {f"vs_{name}": rte.relative_rms(values[gap_core], projected[gap_core])
               for name, values in stored.items()}
        row["runtime_s"] = time.perf_counter() - started
        row["source_trace"] = diagnostics["source_trace"]
        payload["runs"][label] = row
        output.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
        print(label, {k: v for k, v in row.items() if k.startswith("vs_")}, flush=True)
    payload["exact_vs_projected_gap_core_relative_rms"] = rte.relative_rms(
        fields["projected_kelvin_trace"][gap_core], fields["exact_exterior_source"][gap_core])
    payload["completed"] = True
    output.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print("exact vs projected", payload["exact_vs_projected_gap_core_relative_rms"], flush=True)


if __name__ == "__main__":
    main()
