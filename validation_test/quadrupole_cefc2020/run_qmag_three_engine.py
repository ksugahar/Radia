"""Q-mag three-engine comparison: HDiv-MMM against the mixed total/reduced Omega FEM (and reduced-A).

Usage:
    python run_qmag_three_engine.py --fem-mesh <qmag_fem_kelvin.vol> --fem-mesh-report <json>
        --case mu1000 --output <json> [--hdiv-mesh <qmag_h10.vol> | --hdiv-result <run_qmag_hdiv json>]
        [--with-reduced-a] [--fem-order 2] [--hdiv-order 1]

Acceptance (Sugahara, 2026-09-07): HDiv-MMM is validated on this magnet when it agrees with the
repository's own mixed total/reduced Omega formulation on the same coils and the same iron law.
Both are evaluated on the 31 diagonal points of ``qmag_case.observation_points``; the FEM field is
averaged over a small cube around each point (2x2x2 Gauss, half-width 2e-5 m) because the points lie
on the z = 0 mesh seam where a piecewise field is double-valued, while the HDiv-MMM field (a Coulomb
integral) is evaluated pointwise.  The verdict is the relative RMS of the vector field difference over
the 31 points and the B_perp(15 mm) gap, against ``--relative-tolerance`` (default 3 %).
Heavy: run on hibino.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "src" / "radia").exists())
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE))

import ngsolve as ng  # noqa: E402
import numpy as np  # noqa: E402
import radia as rad  # noqa: E402

import qmag_case as Q  # noqa: E402
from run_qmag_hdiv import _load_engines, parse_case  # noqa: E402
from radia.kelvin_identify_ngsolve import detect_kelvin_offset, has_kelvin_identification  # noqa: E402
from radia.vim import mesh_conformity_report  # noqa: E402

SCHEMA = "radia.qmag-cefc2020-three-engine.v1"
GAUSS_2 = (-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _cube_points(centres: np.ndarray, half_width: float) -> np.ndarray:
    offsets = np.array([[gx, gy, gz] for gx in GAUSS_2 for gy in GAUSS_2 for gz in GAUSS_2]) * half_width
    return np.concatenate([c + offsets for c in centres], axis=0)


def _cube_average(values: np.ndarray, centre_count: int) -> np.ndarray:
    return np.asarray(values, dtype=float).reshape(centre_count, 8, 3).mean(axis=1)


def _relative_rms(reference: np.ndarray, candidate: np.ndarray) -> float:
    denominator = float(np.sqrt(np.mean(np.sum(reference * reference, axis=1))))
    if denominator <= 0.0:
        raise RuntimeError("zero reference field")
    return float(np.sqrt(np.mean(np.sum((candidate - reference) ** 2, axis=1))) / denominator)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fem-mesh", type=Path, required=True)
    parser.add_argument("--fem-mesh-report", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--output", type=Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--hdiv-mesh", type=Path, help="solve HDiv-MMM on this conforming HEX mesh")
    group.add_argument("--hdiv-result", type=Path, help="reuse a run_qmag_hdiv.py result of the same case")
    parser.add_argument("--with-reduced-a", action="store_true")
    parser.add_argument("--fem-order", type=int, default=2)
    parser.add_argument("--hdiv-order", type=int, choices=(1, 2), default=1)
    parser.add_argument("--gram-eps", type=float, default=1.0e-10)
    parser.add_argument("--nonlinear-tolerance", type=float, default=2.0e-5)
    parser.add_argument("--nonlinear-maximum-iterations", type=int, default=80)
    parser.add_argument("--source-trace-tolerance", type=float, default=0.05)
    parser.add_argument("--mixed-relaxation", type=float, default=0.3)
    parser.add_argument("--reduced-a-solver", choices=("direct", "bddc", "ams", "auto"), default="direct")
    parser.add_argument("--observation-half-width", type=float, default=2.0e-5)
    parser.add_argument("--relative-tolerance", type=float, default=0.03)
    parser.add_argument("--threads", type=int, default=0)
    options = parser.parse_args(argv)
    if options.threads > 0:
        ng.SetNumThreads(options.threads)
    case = parse_case(options.case)
    fem_path = options.fem_mesh.resolve()
    report_path = options.fem_mesh_report.resolve()
    if not fem_path.is_file() or not report_path.is_file():
        raise FileNotFoundError(f"{fem_path} / {report_path}")
    fem_report = json.loads(report_path.read_text(encoding="utf-8"))
    fem_mesh = ng.Mesh(str(fem_path))
    if not has_kelvin_identification(fem_mesh):
        raise RuntimeError("FEM mesh has no Kelvin point identification")
    kelvin_center = tuple(float(v) for v in detect_kelvin_offset(fem_mesh))
    kelvin_radius = float(fem_report["parameters"]["kelvin_radius_m"])
    points = Q.observation_points()
    fem_points = _cube_points(points, options.observation_half_width)
    for point in fem_points:
        if not fem_mesh(*map(float, point)):
            raise RuntimeError(f"observation point outside the FEM mesh: {point.tolist()}")

    rad.UtiDelAll()
    coil, coil_manifest = Q.build_qmag_coils(case["current_density_A_per_mm2"])
    engines = _load_engines()
    material = Q.load_bh_table() if case["nonlinear"] else float(case["mu_r"])
    fields: dict[str, np.ndarray] = {}
    diagnostics: dict[str, dict] = {}

    # HDiv-MMM: pointwise (the field is a Coulomb integral, smooth at the observation points).
    hdiv_provenance: dict
    if options.hdiv_result is not None:
        prior = json.loads(options.hdiv_result.read_text(encoding="utf-8"))
        if prior["case"]["name"] != case["name"]:
            raise ValueError(f"--hdiv-result is case {prior['case']['name']}, not {case['name']}")
        if not np.allclose(np.asarray(prior["observation_points_m"]), points, atol=1.0e-12):
            raise ValueError("--hdiv-result observation points differ")
        fields["hdiv_mmm"] = np.asarray(prior["B_total_T"], dtype=float)
        diagnostics["hdiv_mmm"] = dict(prior["hdiv"]) | {"reused_result": str(options.hdiv_result.resolve())}
        hdiv_provenance = {"mesh": prior["mesh"], "host": prior.get("host")}
    else:
        hdiv_mesh = ng.Mesh(str(options.hdiv_mesh.resolve()))
        conformity = mesh_conformity_report(hdiv_mesh)
        if not conformity["conforming"]:
            raise RuntimeError(f"non-conforming HDiv mesh: {conformity}")
        field, diag = engines.solve_hdiv(
            hdiv_mesh, coil, material, nonlinear=case["nonlinear"], order=options.hdiv_order,
            gram_eps=options.gram_eps, nonlinear_tolerance=options.nonlinear_tolerance,
            nonlinear_maximum_iterations=options.nonlinear_maximum_iterations, points=points, image=None)
        fields["hdiv_mmm"] = np.asarray(field, dtype=float)
        diagnostics["hdiv_mmm"] = diag
        hdiv_provenance = {"mesh": {"path": str(options.hdiv_mesh.resolve()), "sha256": _sha256(options.hdiv_mesh),
                                    "ne": int(hdiv_mesh.ne), "conformity": conformity}, "host": platform.node()}

    # Mixed total/reduced Omega on the Kelvin mesh (cube-averaged at the seam).
    started = time.perf_counter()
    field, diag = engines.solve_omega(
        fem_mesh, coil, material, nonlinear=case["nonlinear"], order=options.fem_order,
        nonlinear_tolerance=options.nonlinear_tolerance,
        nonlinear_maximum_iterations=options.nonlinear_maximum_iterations, nonlinear_verbose=False,
        kelvin_center=kelvin_center, kelvin_radius=kelvin_radius, points=fem_points,
        source_trace_tolerance=options.source_trace_tolerance, relaxation=options.mixed_relaxation,
        anderson_depth=0, mu_r_initial=1000.0, observation_points=fem_points)
    fields["mixed_total_reduced_omega"] = _cube_average(field, len(points))
    diagnostics["mixed_total_reduced_omega"] = diag | {"wall_s": time.perf_counter() - started}
    if options.with_reduced_a:
        started = time.perf_counter()
        field, diag = engines.solve_reduced_a(
            fem_mesh, coil, material, nonlinear=case["nonlinear"], order=options.fem_order,
            linear_solver=options.reduced_a_solver, relax=0.1, nonlinear_tolerance=options.nonlinear_tolerance,
            nonlinear_maximum_iterations=options.nonlinear_maximum_iterations, nonlinear_verbose=False,
            kelvin_center=kelvin_center, kelvin_radius=kelvin_radius, points=fem_points,
            anderson_depth=0, nu_initial=None, observation_points=fem_points)
        fields["reduced_a"] = _cube_average(field, len(points))
        diagnostics["reduced_a"] = diag | {"wall_s": time.perf_counter() - started}

    names = list(fields)
    pairwise = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            pa, pb = Q.b_perp(fields[a]), Q.b_perp(fields[b])
            pairwise[f"{a}__vs__{b}"] = {
                "relative_rms": _relative_rms(fields[a], fields[b]),
                "B_perp_15mm_T": {a: float(pa[-1]), b: float(pb[-1])},
                "B_perp_15mm_relative_gap": float(abs(pa[-1] - pb[-1]) / abs(pa[-1])),
            }
    key = "hdiv_mmm__vs__mixed_total_reduced_omega"
    failures = []
    if pairwise[key]["relative_rms"] > options.relative_tolerance:
        failures.append(f"HDiv-MMM vs mixed Omega relative RMS {pairwise[key]['relative_rms']:.4f} exceeds "
                        f"{options.relative_tolerance}")
    for name, diag in diagnostics.items():
        if case["nonlinear"] and not diag.get("nonlinear_stats", {}).get("converged", True):
            failures.append(f"{name}: nonlinear loop did not converge")
    report = {
        "schema": SCHEMA, "generated_at_utc": datetime.now(timezone.utc).isoformat(), "host": platform.node(),
        "radia_version": getattr(rad, "__version__", None), "radia_file": rad.__file__, "case": case,
        "fem_mesh": {"path": str(fem_path), "sha256": _sha256(fem_path), "report_sha256": _sha256(report_path),
                     "ne": int(fem_mesh.ne), "kelvin_center": kelvin_center, "kelvin_radius": kelvin_radius,
                     "order": int(options.fem_order)},
        "hdiv": hdiv_provenance, "coils": coil_manifest, "observation_points_m": points.tolist(),
        "observation_half_width_m": options.observation_half_width,
        "fields_T": {k: v.tolist() for k, v in fields.items()},
        "B_perp_T": {k: Q.b_perp(v).tolist() for k, v in fields.items()},
        "engines": diagnostics, "pairwise": pairwise, "failures": failures, "passed": not failures,
    }
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(report, indent=1), encoding="utf-8")
    for k, v in pairwise.items():
        print(f"{case['name']}: {k}: relative RMS {v['relative_rms']:.4f}, B_perp(15 mm) "
              + ", ".join(f"{n} {x:+.5f}" for n, x in v["B_perp_15mm_T"].items())
              + f" (gap {100 * v['B_perp_15mm_relative_gap']:.2f} %)", flush=True)
    print("wrote", options.output, flush=True)
    print("PASSED" if not failures else "FAILED: " + "; ".join(failures), flush=True)
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
