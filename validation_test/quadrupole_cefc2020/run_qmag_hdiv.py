"""HDiv-MMM solve of the CEFC 2020 quadrupole: diagonal field profile and gradient.

Usage:
    python run_qmag_hdiv.py --mesh <qmag_h10.vol> --case mu1000 --output <json>
    python run_qmag_hdiv.py --mesh <vol> --case J3.0 --output <json>          # nonlinear, 3 A/mm^2

Cases: ``mu<mu_r>`` = linear iron of relative permeability mu_r at the 3 A/mm^2 design current,
``J<value>`` = the iron law of ``iron_bh_table.json`` at ``value`` A/mm^2.  The runner reuses the
C-type three-engine ``solve_hdiv`` adapter (the coils are Radia racetracks, no coil mesh),
evaluates B at the 31 diagonal points and reports B_perp(15 mm) and the gradient.  The checks
are internal: finite field, the odd symmetry of B_perp along the diagonal (quadrupole), and a
converged nonlinear loop.  Mesh convergence is judged across the ``qmag_h*.vol`` series.
Heavy: run on hibino (one job at a time).
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
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
from radia.vim import mesh_conformity_report  # noqa: E402

CTYPE_RUNNER_PATH = REPO / "validation_test" / "c_type_three_engine" / "run_three_engine.py"
SCHEMA = "radia.qmag-cefc2020-hdiv-result.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()
def implementation_identity() -> dict:
    """Bind a measurement to the exact implementation that produced it.

    A version string does not separate two candidate builds of the same version, so the record also
    carries the SHA-256 of the loaded native extension and, when the run resolves from a source tree,
    the git commit and worktree state of that tree.
    """
    import subprocess

    module = Path(rad.__file__).resolve()
    native = module.parent / "_radia_pybind.pyd"
    identity = {
        "radia_version": getattr(rad, "__version__", None),
        "radia_module": str(module),
        "native_extension": str(native) if native.is_file() else None,
        "native_sha256": _sha256(native) if native.is_file() else None,
        "native_mtime_utc": (datetime.fromtimestamp(native.stat().st_mtime, timezone.utc).isoformat()
                             if native.is_file() else None),
    }
    for candidate in (module.parents[2], Path(__file__).resolve().parents[2]):
        git = candidate / ".git"
        if not (git.is_dir() or git.is_file()):
            continue
        try:
            head = subprocess.run(["git", "-C", str(candidate), "rev-parse", "HEAD"],
                                  capture_output=True, text=True, check=True).stdout.strip()
            dirty = subprocess.run(["git", "-C", str(candidate), "status", "--porcelain"],
                                   capture_output=True, text=True, check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            continue
        identity["source_tree"] = str(candidate)
        identity["source_commit"] = head
        identity["source_tracked_clean"] = not dirty
        break
    return identity


def _load_engines():
    spec = importlib.util.spec_from_file_location("_radia_c_type_engines", CTYPE_RUNNER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {CTYPE_RUNNER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_case(text: str) -> dict:
    text = text.strip()
    if text.startswith("mu"):
        mu_r = float(text[2:])
        if not mu_r > 1.0:
            raise ValueError("mu_r must exceed 1")
        return {"name": text, "nonlinear": False, "mu_r": mu_r,
                "current_density_A_per_mm2": Q.DESIGN_CURRENT_DENSITY_A_PER_MM2}
    if text.startswith("J"):
        j = float(text[1:])
        if not j > 0.0:
            raise ValueError("J must be positive")
        return {"name": text, "nonlinear": True, "mu_r": None, "current_density_A_per_mm2": j}
    raise ValueError("case must be mu<mu_r> or J<A/mm^2>")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--hdiv-order", type=int, choices=(1, 2), default=1)
    parser.add_argument("--gram-eps", type=float, default=1.0e-10)
    parser.add_argument("--nonlinear-tolerance", type=float, default=2.0e-5)
    parser.add_argument("--nonlinear-maximum-iterations", type=int, default=80)
    parser.add_argument("--symmetry-tolerance", type=float, default=1.0e-3,
                        help="allowed even part of B_perp along the diagonal, relative to its maximum")
    parser.add_argument("--threads", type=int, default=0)
    options = parser.parse_args(argv)
    if options.threads > 0:
        ng.SetNumThreads(options.threads)
    case = parse_case(options.case)
    mesh_path = options.mesh.resolve()
    if not mesh_path.is_file():
        raise FileNotFoundError(mesh_path)
    mesh = ng.Mesh(str(mesh_path))
    conformity = mesh_conformity_report(mesh)
    if not conformity["conforming"]:
        raise RuntimeError(f"non-conforming iron mesh: {conformity}")
    if set(mesh.GetMaterials()) != {"iron"}:
        raise RuntimeError(f"expected the single material 'iron', got {list(mesh.GetMaterials())}")

    rad.UtiDelAll()
    coil, coil_manifest = Q.build_qmag_coils(case["current_density_A_per_mm2"])
    points = Q.observation_points()
    coil_field = np.asarray(rad.Fld(coil, "b", points), dtype=float)
    if not np.isfinite(coil_field).all():
        raise RuntimeError("coil field is non-finite at an observation point")
    engines = _load_engines()
    material = Q.load_bh_table() if case["nonlinear"] else float(case["mu_r"])
    started = time.perf_counter()
    field, diagnostics = engines.solve_hdiv(
        mesh, coil, material, nonlinear=case["nonlinear"], order=options.hdiv_order, gram_eps=options.gram_eps,
        nonlinear_tolerance=options.nonlinear_tolerance,
        nonlinear_maximum_iterations=options.nonlinear_maximum_iterations, points=points, image=None,
    )
    wall = time.perf_counter() - started
    field = np.asarray(field, dtype=float)
    perp = Q.b_perp(field)
    value = float(perp[-1])
    even_part = float(np.max(np.abs(perp + perp[::-1])) / np.max(np.abs(perp)))
    failures = []
    if not np.isfinite(field).all():
        failures.append("non-finite field")
    if case["nonlinear"] and not diagnostics["nonlinear_stats"].get("converged", False):
        failures.append("nonlinear Picard loop did not converge")
    if even_part > options.symmetry_tolerance:
        failures.append(f"B_perp is not odd along the diagonal (even part {even_part:.2e})")
    report = {
        "schema": SCHEMA, "generated_at_utc": datetime.now(timezone.utc).isoformat(), "host": platform.node(),
        "implementation": implementation_identity(),
        "radia_version": getattr(rad, "__version__", None), "radia_file": rad.__file__,
        "case": case,
        "mesh": {"path": str(mesh_path), "sha256": _sha256(mesh_path), "ne": int(mesh.ne), "nv": int(mesh.nv),
                 "nbnd": int(mesh.GetNE(ng.BND)), "conformity": conformity},
        "coils": coil_manifest,
        "hdiv": diagnostics | {"wall_s": wall, "order": int(options.hdiv_order), "gram_eps": float(options.gram_eps)},
        "observation_points_m": points.tolist(), "B_total_T": field.tolist(), "B_coil_T": coil_field.tolist(),
        "B_perp_T": perp.tolist(), "B_perp_15mm_T": value, "gradient_T_per_m": value / Q.OBSERVATION_RADIUS_M,
        "coil_only_gradient_T_per_m": float(Q.b_perp(coil_field)[-1]) / Q.OBSERVATION_RADIUS_M,
        "even_part_of_B_perp": even_part, "failures": failures, "passed": not failures,
    }
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(report, indent=1), encoding="utf-8")
    print(f"{case['name']} on {mesh_path.name}: B_perp(15 mm) {value:+.5f} T, gradient {value / Q.OBSERVATION_RADIUS_M:+.3f} T/m "
          f"(coil only {report['coil_only_gradient_T_per_m']:+.3f}), even part {even_part:.1e}, "
          f"ndof {diagnostics['ndof']}, CG {diagnostics['linear_iterations']}, wall {wall:.0f} s", flush=True)
    print("wrote", options.output, flush=True)
    print("PASSED" if not failures else "FAILED: " + "; ".join(failures), flush=True)
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
