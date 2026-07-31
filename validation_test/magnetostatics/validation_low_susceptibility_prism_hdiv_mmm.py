"""HDiv-MMM convergence for a low-susceptibility rectangular prism."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import ngsolve as ng
import numpy as np
from netgen.csg import CSGeometry, OrthoBrick, Pnt


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from radia import vim  # noqa: E402
from radia.analytical_formulas import (  # noqa: E402
    linear_prism_average_flux_density,
    rectangular_prism_demag_factors,
)


OUT_JSON = Path(__file__).with_name(
    "validation_low_susceptibility_prism_hdiv_mmm_summary.json"
)
SIDES_M = (0.011, 0.009, 0.027)
MU_R = 1.04
H_EXT_A_PER_M = 125000.0
MU_0 = 4.0e-7 * math.pi


def _identity_digest() -> str:
    identity = {
        "geometry": "rectangular_prism",
        "side_lengths_m": SIDES_M,
        "material": {"law": "linear_isotropic", "mu_r": MU_R},
        "applied_field_A_per_m": [H_EXT_A_PER_M, 0.0, 0.0],
    }
    encoded = json.dumps(identity, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("ascii")).hexdigest()


def _solver_source_head() -> str:
    completed = subprocess.run(
        [
            "git",
            "-C",
            str(ROOT),
            "log",
            "-1",
            "--format=%H",
            "--",
            "src/radia/vim/_solve.py",
            "src/radia/vim/_vim.py",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def solve_level(maxh_m: float) -> dict[str, object]:
    sx, sy, sz = SIDES_M
    geometry = CSGeometry()
    geometry.Add(
        OrthoBrick(
            Pnt(-0.5 * sx, -0.5 * sy, -0.5 * sz),
            Pnt(+0.5 * sx, +0.5 * sy, +0.5 * sz),
        ).mat("magnetic_body")
    )
    mesh_started = time.perf_counter()
    mesh = ng.Mesh(geometry.GenerateMesh(maxh=maxh_m))
    mesh_duration = time.perf_counter() - mesh_started

    solve_started = time.perf_counter()
    with ng.TaskManager():
        result = vim.Solve(
            mesh,
            mu_r=MU_R,
            H_ext=ng.CoefficientFunction((H_EXT_A_PER_M, 0.0, 0.0)),
            order=1,
            tol=1.0e-9,
            maxit=4000,
        )
    solve_duration = time.perf_counter() - solve_started
    m_avg = np.asarray(result["M_avg"], dtype=float)
    b_avg = MU_0 * MU_R / (MU_R - 1.0) * m_avg
    return {
        "maxh_m": maxh_m,
        "tetrahedra": int(mesh.ne),
        "vertices": int(mesh.nv),
        "ndof": int(result["ndof"]),
        "iterations": int(result["iters"]),
        "M_avg_A_per_m": m_avg.tolist(),
        "B_avg_T": b_avg.tolist(),
        "timing_s": {"mesh": mesh_duration, "solve": solve_duration},
    }


def main() -> int:
    levels = [solve_level(maxh) for maxh in (0.006, 0.003, 0.0015)]
    nx = rectangular_prism_demag_factors(*SIDES_M)[0]
    analytic_bx = linear_prism_average_flux_density(H_EXT_A_PER_M, MU_R, nx)
    finest_bx = float(levels[-1]["B_avg_T"][0])
    previous_bx = float(levels[-2]["B_avg_T"][0])
    analytic_relative_error = abs(finest_bx - analytic_bx) / abs(analytic_bx)
    last_relative_change = abs(finest_bx - previous_bx) / abs(finest_bx)
    transverse_ratio = float(
        np.linalg.norm(levels[-1]["B_avg_T"][1:]) / abs(finest_bx)
    )
    checks = {
        "tet_mesh_refines_monotonically": all(
            levels[index]["tetrahedra"] < levels[index + 1]["tetrahedra"]
            for index in range(len(levels) - 1)
        ),
        "last_refinement_change_below_1e-4": last_relative_change < 1.0e-4,
        "aharoni_low_susceptibility_error_below_1e-3": (
            analytic_relative_error < 1.0e-3
        ),
        "transverse_average_below_1e-5": transverse_ratio < 1.0e-5,
        "all_linear_solves_converged": all(
            0 < int(level["iterations"]) < 4000 for level in levels
        ),
    }
    formula_path = ROOT / "src" / "radia" / "analytical_formulas" / "rectangular_prism.py"
    artifact = {
        "schema": "radia.validation.magnetostatics.low-chi-prism-hdiv-mmm.v1",
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "execution_version": {
            "producer": Path(__file__).name,
            "python": platform.python_version(),
            "ngsolve": getattr(ng, "__version__", "unknown"),
            "radia_solver_source_head": _solver_source_head(),
            "analytic_formula_sha256": _sha256(formula_path),
        },
        "physics_identity_digest": _identity_digest(),
        "identity": {
            "geometry": "rectangular_prism",
            "side_lengths_m": SIDES_M,
            "material": {"law": "linear_isotropic", "mu_r": MU_R},
            "applied_field_A_per_m": [H_EXT_A_PER_M, 0.0, 0.0],
        },
        "discretization": "BDM1 HDiv-MMM on flat tetrahedra",
        "levels": levels,
        "analytic_reference": {
            "name": "Aharoni magnetometric rectangular-prism factor",
            "Nx": nx,
            "B_avg_x_T": analytic_bx,
            "scope": "low-susceptibility convergence reference",
        },
        "metrics": {
            "last_relative_refinement_change": last_relative_change,
            "analytic_relative_error": analytic_relative_error,
            "transverse_to_axial_average_ratio": transverse_ratio,
        },
        "timing_breakdown_s": {
            f"solve_maxh_{level['maxh_m']}": level["timing_s"]["solve"]
            for level in sorted(
                levels, key=lambda item: item["timing_s"]["solve"], reverse=True
            )
        },
        "checks": checks,
        "pass": all(checks.values()),
    }
    OUT_JSON.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(artifact, indent=2))
    return 0 if artifact["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
