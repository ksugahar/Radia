"""Cross-check finite-domain H1-Hodge BDM against an axisymmetric Q2 disk.

The reference is the static axisymmetric ``radia.axifem`` Q2 solution.  The
three-dimensional lane uses a separately generated TET air/body mesh, a
body-restricted BDM space, and ``H1HodgeDemagOperator``.  It is deliberately a
heavy validation runner rather than a unit test.

No mesh artifact is tracked; every mesh is regenerated in memory.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import ngsolve as ng
import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "src"
DEFAULT_OUTPUT = HERE / "results_h1_hodge_bdm2_disk.json"

for path in (SRC, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import radia  # noqa: E402
from radia import vim  # noqa: E402
import validate_magnetic_conductor_disk as disk  # noqa: E402


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _git_head() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
    ).strip()


def _git_dirty() -> bool:
    return bool(
        subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=REPO, text=True
        ).strip()
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source_fingerprints() -> dict[str, str]:
    paths = (
        Path(__file__).resolve(),
        HERE / "validate_magnetic_conductor_disk.py",
        SRC / "radia" / "vim" / "_vim.py",
        SRC / "radia" / "vim" / "_eddy_hybrid.py",
        SRC / "radia" / "axifem.pyd",
    )
    return {
        path.relative_to(REPO).as_posix(): _sha256(path)
        for path in paths
    }


def _solve_case(
    mesh: ng.Mesh,
    *,
    maxh_m: float,
    hdiv_order: int,
    h1_order: int,
    reference: float,
) -> dict[str, object]:
    started = time.perf_counter()
    body = mesh.Materials("conductor")
    hdiv = ng.HDiv(mesh, order=hdiv_order, definedon=body)
    h1 = ng.H1(mesh, order=h1_order, dirichlet="outer")
    with ng.TaskManager():
        demag = vim.H1HodgeDemagOperator(
            hdiv,
            h1,
            definedon=body,
            boundary_contract="finite-dirichlet-disk-accuracy",
        )
        reduction = vim.NgsolveHDivMMMResponseReduction(
            mesh,
            hdiv,
            mu_r=disk.MU_R,
            external_fields=(ng.CF((0.0, 0.0, disk.H0_A_PER_M)),),
            external_names=("uniform_Hz",),
            max_modes=1,
            materials="conductor",
            demag_operator=demag,
            intorder=4,
            solve_tol=1.0e-10,
            parent_family="BDM",
            parent_order=hdiv_order,
        )

    solution = reduction.solve()
    average_mz = float(solution.average_magnetization[0, 2])
    normalized_bz = (
        disk.MU_R
        / (disk.MU_R - 1.0)
        * average_mz
        / disk.H0_A_PER_M
    )
    generalized_demag = float(
        reduction.demag[0, 0] / reduction.mass[0, 0]
    )
    generation = reduction.basis_generation
    return {
        "maxh_m": float(maxh_m),
        "tetrahedra": int(mesh.ne),
        "hdiv_family": "BDM",
        "hdiv_order": int(hdiv_order),
        "hdiv_dofs": int(hdiv.ndof),
        "hdiv_active_dofs": int(sum(hdiv.FreeDofs())),
        "h1_order": int(h1_order),
        "h1_dofs": int(h1.ndof),
        "h1_free_dofs": int(sum(h1.FreeDofs())),
        "normalized_Bz": float(normalized_bz),
        "reference_relative_error": float(
            abs(normalized_bz - reference) / abs(reference)
        ),
        "reduced_demag_generalized_eigenvalue": generalized_demag,
        "snapshot_backend": generation["snapshot_backend"],
        "snapshot_iterations": generation["snapshot_iterations"],
        "max_snapshot_relative_residual": float(
            generation["max_snapshot_relative_residual"]
        ),
        "response_relative_energy_error": float(
            generation["max_response_relative_energy_error"]
        ),
        "reduced_solve_relative_residual": float(
            solution.residual_relative_norm
        ),
        "operator": demag.Diagnostics(),
        "wall_s": time.perf_counter() - started,
    }


def run() -> dict[str, object]:
    started = time.perf_counter()
    source_start = _source_fingerprints()
    head_start = _git_head()
    reference_row = disk.solve_axisymmetric_q2(0.0, fine=True)
    reference = float(reference_row["normalized_Bz"][0])

    mesh_started = time.perf_counter()
    coarse_mesh = disk._full_3d_mesh(0.002)
    coarse_mesh_s = time.perf_counter() - mesh_started
    coarse_bdm1 = _solve_case(
        coarse_mesh,
        maxh_m=0.002,
        hdiv_order=1,
        h1_order=2,
        reference=reference,
    )
    coarse_bdm2 = _solve_case(
        coarse_mesh,
        maxh_m=0.002,
        hdiv_order=2,
        h1_order=3,
        reference=reference,
    )

    mesh_started = time.perf_counter()
    fine_mesh = disk._full_3d_mesh(0.001)
    fine_mesh_s = time.perf_counter() - mesh_started
    fine_bdm2_h1p3 = _solve_case(
        fine_mesh,
        maxh_m=0.001,
        hdiv_order=2,
        h1_order=3,
        reference=reference,
    )
    fine_bdm2_h1p4 = _solve_case(
        fine_mesh,
        maxh_m=0.001,
        hdiv_order=2,
        h1_order=4,
        reference=reference,
    )
    cases = [coarse_bdm1, coarse_bdm2, fine_bdm2_h1p3, fine_bdm2_h1p4]
    source_end = _source_fingerprints()
    bdm2_errors = [
        coarse_bdm2["reference_relative_error"],
        fine_bdm2_h1p3["reference_relative_error"],
        fine_bdm2_h1p4["reference_relative_error"],
    ]
    checks = {
        "source_fingerprints_stable_during_run": source_start == source_end,
        "axisymmetric_static_reference_is_real": (
            abs(float(reference_row["normalized_Bz"][1])) < 1.0e-14
        ),
        "all_response_solves_converged": max(
            row["max_snapshot_relative_residual"] for row in cases
        )
        < 1.0e-8,
        "all_reduced_material_solves_converged": max(
            row["reduced_solve_relative_residual"] for row in cases
        )
        < 1.0e-12,
        "all_hodge_modes_are_contractive": all(
            -1.0e-8 <= row["reduced_demag_generalized_eigenvalue"]
            <= 1.0 + 1.0e-5
            for row in cases
        ),
        "coarse_bdm2_improves_on_coarse_bdm1": (
            coarse_bdm2["reference_relative_error"]
            < coarse_bdm1["reference_relative_error"]
        ),
        "bdm2_h_and_h1_p_ladder_is_strict": bdm2_errors == sorted(
            bdm2_errors, reverse=True
        ),
        "fine_bdm2_h1p4_error_below_one_percent": (
            fine_bdm2_h1p4["reference_relative_error"] < 0.01
        ),
    }
    return {
        "schema": "radia.validation.h1-hodge-bdm2-axisymmetric-disk.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "tool_versions": {
            "radia": getattr(
                radia, "__version__", _package_version("radia-ngsolve")
            ),
            "radia_source_head": head_start,
            "radia_source_head_end": _git_head(),
            "radia_source_dirty": _git_dirty(),
            "python": platform.python_version(),
            "ngsolve": _package_version("ngsolve"),
            "numpy": _package_version("numpy"),
        },
        "source_fingerprints": source_start,
        "identity": {
            "geometry": "solid circular disk in an air box",
            "radius_m": disk.RADIUS_M,
            "thickness_m": disk.THICKNESS_M,
            "mu_r": disk.MU_R,
            "excitation": "uniform axial H",
            "meshes_tracked": False,
        },
        "axisymmetric_q2_static_reference": reference_row,
        "h1_hodge_3d_cases": cases,
        "checks": checks,
        "pass": all(checks.values()),
        "timing_s": {
            "axisymmetric_reference": float(reference_row["wall_s"]),
            "coarse_mesh": coarse_mesh_s,
            "fine_mesh": fine_mesh_s,
            "case_solves": float(sum(row["wall_s"] for row in cases)),
            "total": time.perf_counter() - started,
        },
        "claim_boundary": {
            "established": (
                "on this disk, the finite-domain 3-D BDM2/H1-Hodge response "
                "converges monotonically toward the independent axisymmetric "
                "Q2 static observable and reaches below one-percent error"
            ),
            "not_established": (
                "an exact open-boundary H1 formulation, universal BDM2 "
                "superiority, or transient mixed-solver accuracy"
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {"output": str(output), "pass": result["pass"], "checks": result["checks"]},
            indent=2,
        )
    )
    return 0 if result["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
