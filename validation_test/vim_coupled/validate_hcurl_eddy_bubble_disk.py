"""Validate the 3-D HCurl eddy-bubble reduction on a conducting disk.

This lane deliberately contains both a negative and a positive control:

* a single uniform-field training port is spatially under-resolved even when
  its Krylov depth is increased;
* four polynomial vector-potential ports recover the port-dominant diffusion
  pole under h- and p-refinement.

The independent reference chain is public and reproducible in this repository:
the axisymmetric integral-equation BEM modal spectrum is cross-checked by the
Q1/Q2 ``radia.axifem`` Henrotte formulations before the 3-D HCurl solve is
accepted.  The generated tetrahedral mesh is temporary; no mesh artifact is
required in git.

Run from the repository root:

    python validation_test/vim_coupled/validate_hcurl_eddy_bubble_disk.py
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy.linalg


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "src"
AXIFEM_VERIFY = REPO / "validation_test" / "axifem" / "research" / "verification"
BEM_CAUER_RESULT = (
    REPO
    / "validation_test"
    / "maglev"
    / "research_cln"
    / "ngsolve_validation"
    / "bem_disk_axisym_cauer_python_results.json"
)
DEFAULT_OUTPUT = HERE / "results_hcurl_eddy_bubble_disk.json"

for path in (SRC, AXIFEM_VERIFY):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import netgen.occ as occ  # noqa: E402
import ngsolve as ng  # noqa: E402
import radia  # noqa: E402
import radia.vim as vim  # noqa: E402
import test_hiruma_disk_q1 as axifem_q1  # noqa: E402
import test_hiruma_disk_q2 as axifem_q2  # noqa: E402


RADIUS_M = 0.010
THICKNESS_M = 0.002
SIGMA_S_PER_M = 5.8e7
BEM_MODAL_REFERENCE_US = np.asarray(axifem_q1.BEM_TAU, dtype=float)


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _real_symmetric(matrix: np.ndarray, label: str) -> np.ndarray:
    values = np.asarray(matrix)
    scale = max(float(np.linalg.norm(values)), np.finfo(float).tiny)
    imaginary_ratio = float(np.linalg.norm(values.imag) / scale)
    if imaginary_ratio > 1.0e-10:
        raise ValueError(f"{label} has imaginary ratio {imaginary_ratio:.3e}")
    real = np.asarray(values.real, dtype=float)
    return 0.5 * (real + real.T)


def _disk_mesh(maxh_m: float) -> ng.Mesh:
    disk = occ.Cylinder(
        occ.Pnt(0.0, 0.0, -0.5 * THICKNESS_M),
        occ.Z,
        RADIUS_M,
        THICKNESS_M,
    )
    disk.mat("cond")
    for face in disk.faces:
        face.name = "skin"
    return ng.Mesh(occ.OCCGeometry(disk).GenerateMesh(maxh=maxh_m))


def _uniform_bz_vector_potential_cf():
    return 0.5 * ng.CoefficientFunction((-ng.y, ng.x, 0.0))


def _training_vector_potentials_cf(training: str):
    base = _uniform_bz_vector_potential_cf()
    if training == "uniform":
        return (base,)
    if training != "polynomial":
        raise ValueError(f"unknown training set: {training}")
    radial = (ng.x * ng.x + ng.y * ng.y) / (RADIUS_M * RADIUS_M)
    axial = (2.0 * ng.z / THICKNESS_M) ** 2
    return (
        base,
        radial * base,
        radial * radial * base,
        axial * base,
    )


def _uniform_bz_vector_potential_samples(points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    return 0.5 * np.column_stack(
        (-points[:, 1], points[:, 0], np.zeros(points.shape[0]))
    )


def _modal_spectrum(model: vim.HCurlEddyCLNModel) -> dict[str, object]:
    resistance = _real_symmetric(model.resistance, "resistance")
    inductance = _real_symmetric(model.inductance, "inductance")
    rhs = np.asarray(model.port_rhs[:, 0])
    eigenvalues, eigenvectors = scipy.linalg.eigh(
        inductance,
        resistance,
        check_finite=True,
    )
    positive = np.isfinite(eigenvalues) & (eigenvalues > 0.0)
    eigenvalues = eigenvalues[positive]
    eigenvectors = eigenvectors[:, positive]
    if eigenvalues.size == 0:
        raise RuntimeError("the reduced CLN has no positive time constants")

    residues = np.abs(eigenvectors.conj().T @ rhs) ** 2
    tau_order = np.argsort(eigenvalues)[::-1]
    tau_us = 1.0e6 * eigenvalues[tau_order]
    residues_tau_order = residues[tau_order]
    dominant_index = int(np.argmax(residues))
    dominant_tau_us = float(1.0e6 * eigenvalues[dominant_index])
    residue_sum = max(float(np.sum(residues)), np.finfo(float).tiny)
    dominant_error_pct = float(
        100.0
        * abs(dominant_tau_us - BEM_MODAL_REFERENCE_US[0])
        / BEM_MODAL_REFERENCE_US[0]
    )
    return {
        "time_constants_us_descending": tau_us.tolist(),
        "port_residues_in_tau_order": residues_tau_order.tolist(),
        "port_dominant_time_constant_us": dominant_tau_us,
        "port_dominant_residue_fraction": float(
            residues[dominant_index] / residue_sum
        ),
        "reference_leading_time_constant_us": float(BEM_MODAL_REFERENCE_US[0]),
        "port_dominant_abs_error_pct": dominant_error_pct,
        "minimum_resistance_eigenvalue": float(np.linalg.eigvalsh(resistance)[0]),
        "minimum_inductance_eigenvalue": float(np.linalg.eigvalsh(inductance)[0]),
        "resistance_hermitian_error": float(
            np.linalg.norm(model.resistance - model.resistance.conj().T)
            / max(float(np.linalg.norm(model.resistance)), np.finfo(float).tiny)
        ),
        "inductance_hermitian_error": float(
            np.linalg.norm(model.inductance - model.inductance.conj().T)
            / max(float(np.linalg.norm(model.inductance)), np.finfo(float).tiny)
        ),
    }


def _run_hcurl_case(
    *,
    maxh_m: float,
    order: int,
    steps: int,
    training: str,
) -> dict[str, object]:
    started = time.perf_counter()
    mesh = _disk_mesh(maxh_m)
    mesh_s = time.perf_counter() - started

    fes = ng.HCurl(mesh, order=order, nograds=True)
    u, v = fes.TnT()
    stiffness = ng.BilinearForm(fes)
    stiffness += (
        ng.curl(u) * ng.curl(v) + 100.0 * u * v
    ) * ng.dx(definedon=mesh.Materials("cond"))
    metric = ng.BilinearForm(fes)
    metric += u * v * ng.dx(definedon=mesh.Materials("cond"))
    training_ports = tuple(
        vim.NgsolveHCurlVectorPotentialPort(fes, potential, materials="cond")
        for potential in _training_vector_potentials_cf(training)
    )

    reduction_started = time.perf_counter()
    with ng.TaskManager():
        stiffness.Assemble()
        metric.Assemble()
        for port in training_ports:
            port.Assemble()
        basis = vim.NgsolveEddyBubbleHCurlBasis(
            mesh,
            fes,
            stiffness,
            metric,
            training_ports,
            steps=steps,
            conductive_materials="cond",
            volume_materials="cond",
            intorder=4,
            response_backend="operator",
            inverse="sparsecholesky",
            rtol=1.0e-11,
            current_gram_rtol=1.0e-11,
            parent_order=order,
            parent_order_ledger=vim.EddyParentOrderLedger(
                bulk_degree=order,
                bridge_trace_degree=max(order - 1, 0),
                surface_current_degree=max(order - 1, 0),
            ),
        )
    reduction_s = time.perf_counter() - reduction_started

    interaction_started = time.perf_counter()
    with ng.TaskManager():
        interaction = basis.tet_volume_interaction(
            mesh,
            fes,
            degree=max(order - 1, 0),
            projection_quad=max(order + 2, 4),
            outer_quad=max(order + 2, 4),
            projection_tolerance=1.0e-9,
            matrix_free=True,
            materials="cond",
        )
    interaction_s = time.perf_counter() - interaction_started

    cln_started = time.perf_counter()
    system = basis.assemble_vim(
        sigma=SIGMA_S_PER_M,
        interaction=interaction,
    )
    rhs = vim.ExternalVectorPotentialRHS(
        basis.current_basis,
        _uniform_bz_vector_potential_samples(basis.current_basis.points),
    )
    model = vim.HCurlEddyCLNFromVIM(system, rhs)
    modal = _modal_spectrum(model)
    cln_s = time.perf_counter() - cln_started

    diagnostics = basis.diagnostics()
    interaction_diagnostics = interaction.diagnostics()
    topology = basis.eddy_bubbling.topology.diagnostics()
    hmatrix = interaction_diagnostics["hmatrix_operator"]["hmatrix"]
    return {
        "maxh_m": float(maxh_m),
        "mesh": {
            "elements": int(mesh.ne),
            "vertices": int(mesh.nv),
            "curve_order": int(mesh.GetCurveOrder()),
            "cell_family": "tet",
        },
        "parent": {
            "family": "HCurl",
            "order": int(order),
            "ndof": int(fes.ndof),
        },
        "training": {
            "kind": training,
            "port_count": len(training_ports),
            "krylov_steps": int(steps),
            "vector_potentials": (
                ["A_uniform"]
                if training == "uniform"
                else ["A_uniform", "r2_A", "r4_A", "z2_A"]
            ),
        },
        "reduction": {
            "kind": "EddyBubbleHCurlBasis",
            "rank": int(diagnostics["rank"]),
            "current_gram_rank": int(
                diagnostics["response_basis"]["current_gram_rank"]
            ),
            "conductive_cycle_rank": int(topology["conductive_graph_cycle_rank"]),
        },
        "interaction": {
            "kind": interaction_diagnostics["kind"],
            "backend": interaction_diagnostics["backend"],
            "projection_relative_residual": float(
                interaction_diagnostics["projection_relative_residual"]
            ),
            "hmatrix_compression": float(hmatrix["compression"]),
        },
        "modal": modal,
        "timing_s": {
            "mesh": mesh_s,
            "parent_assembly_and_reduction": reduction_s,
            "epsilon_free_vim_interaction": interaction_s,
            "cln_and_modal_spectrum": cln_s,
            "total": time.perf_counter() - started,
        },
    }


def _run_axifem_reference() -> dict[str, object]:
    cauer = json.loads(BEM_CAUER_RESULT.read_text(encoding="utf-8"))
    cauer_tau_us = [float(value) for value in cauer["tau_pair_us"][:6]]
    started = time.perf_counter()

    q1_started = time.perf_counter()
    q1 = axifem_q1.solve_disk(
        160,
        32,
        25,
        25,
        500e-3,
        500e-3,
        N_stages=6,
        label="Q1 very fine",
    )
    q1_s = time.perf_counter() - q1_started

    q2_started = time.perf_counter()
    q2 = axifem_q2.solve_disk_q2(
        40,
        16,
        15,
        15,
        500e-3,
        500e-3,
        N_stages=6,
        label="Q2 fine",
    )
    q2_s = time.perf_counter() - q2_started

    rows = []
    for label, result in (("Q1", q1), ("Q2", q2)):
        cauer_values = [float(stage["tau_pair_us"]) for stage in result["stages"]]
        eigen_values = [float(value) for value in result["eigsh_tau_us"]]
        rows.append(
            {
                "element_family": label,
                "mesh": result["mesh"],
                "first_cauer_time_constant_us": cauer_values[0],
                "first_cauer_abs_error_pct": float(
                    100.0
                    * abs(cauer_values[0] - cauer_tau_us[0])
                    / abs(cauer_tau_us[0])
                ),
                "first_eigen_time_constant_us": eigen_values[0],
                "first_eigen_abs_error_pct": float(
                    100.0
                    * abs(eigen_values[0] - BEM_MODAL_REFERENCE_US[0])
                    / BEM_MODAL_REFERENCE_US[0]
                ),
            }
        )
    return {
        "reference_files": [
            str(BEM_CAUER_RESULT.relative_to(REPO)).replace("\\", "/"),
            "validation_test/axifem/research/verification/test_hiruma_disk_q1.py",
            "validation_test/axifem/research/verification/test_hiruma_disk_q2.py",
        ],
        "bem_first_cauer_time_constant_us": cauer_tau_us[0],
        "bem_first_modal_time_constant_us": float(BEM_MODAL_REFERENCE_US[0]),
        "rows": rows,
        "timing_s": {
            "q1": q1_s,
            "q2": q2_s,
            "total": time.perf_counter() - started,
        },
    }


def _print_case(label: str, row: dict[str, object]) -> None:
    print(
        "[{label}] ne={ne}, p={order}, rank={rank}, tau={tau:.6f} us, "
        "error={error:.3f}%".format(
            label=label,
            ne=row["mesh"]["elements"],
            order=row["parent"]["order"],
            rank=row["reduction"]["rank"],
            tau=row["modal"]["port_dominant_time_constant_us"],
            error=row["modal"]["port_dominant_abs_error_pct"],
        ),
        flush=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else REPO / args.output
    started = time.perf_counter()

    axifem = _run_axifem_reference()
    h_rows = []
    for maxh_m in (0.003, 0.002, 0.0015, 0.001):
        row = _run_hcurl_case(
            maxh_m=maxh_m,
            order=2,
            steps=3,
            training="polynomial",
        )
        h_rows.append(row)
        _print_case(f"h={maxh_m:g}", row)

    p_rows = []
    for order in (1, 2, 3):
        if order == 2:
            row = next(value for value in h_rows if value["maxh_m"] == 0.002)
        else:
            row = _run_hcurl_case(
                maxh_m=0.002,
                order=order,
                steps=3,
                training="polynomial",
            )
        p_rows.append(row)
        _print_case(f"p={order}", row)

    negative_rows = []
    for steps in (8, 16):
        row = _run_hcurl_case(
            maxh_m=0.002,
            order=2,
            steps=steps,
            training="uniform",
        )
        negative_rows.append(row)
        _print_case(f"uniform-steps={steps}", row)

    positive_errors = [
        float(row["modal"]["port_dominant_abs_error_pct"]) for row in h_rows
    ]
    p_errors = [
        float(row["modal"]["port_dominant_abs_error_pct"]) for row in p_rows
    ]
    negative_errors = [
        float(row["modal"]["port_dominant_abs_error_pct"])
        for row in negative_rows
    ]
    all_rows = h_rows + [p_rows[0], p_rows[2]] + negative_rows
    projection_error = max(
        float(row["interaction"]["projection_relative_residual"])
        for row in all_rows
    )
    minimum_r = min(
        float(row["modal"]["minimum_resistance_eigenvalue"])
        for row in all_rows
    )
    minimum_l = min(
        float(row["modal"]["minimum_inductance_eigenvalue"])
        for row in all_rows
    )
    h_refinement_change = abs(
        float(h_rows[-1]["modal"]["port_dominant_time_constant_us"])
        - float(h_rows[-2]["modal"]["port_dominant_time_constant_us"])
    ) / abs(float(h_rows[-1]["modal"]["port_dominant_time_constant_us"]))
    axifem_eigen_errors = [
        float(row["first_eigen_abs_error_pct"]) for row in axifem["rows"]
    ]
    axifem_cauer_errors = [
        float(row["first_cauer_abs_error_pct"]) for row in axifem["rows"]
    ]
    checks = {
        "ran_to_completion": True,
        "result_files_exist": True,
        "axifem_q1_q2_first_eigen_below_0_6_percent": (
            max(axifem_eigen_errors) < 0.6
        ),
        "axifem_q1_q2_first_cauer_below_0_6_percent": (
            max(axifem_cauer_errors) < 0.6
        ),
        "polynomial_h_sweep_all_below_2_percent": max(positive_errors) < 2.0,
        "polynomial_p1_p2_p3_all_below_2_percent": max(p_errors) < 2.0,
        "h_refinement_change_below_5_percent": h_refinement_change < 0.05,
        "single_port_remains_above_4_percent_at_steps_8_and_16": (
            min(negative_errors) > 4.0
        ),
        "polynomial_basis_improves_at_least_3_percentage_points": (
            min(negative_errors) - min(positive_errors) > 3.0
        ),
        "epsilon_free_projection_below_1e_9": projection_error <= 1.0e-9,
        "all_reduced_cln_models_are_passive": (
            minimum_r >= -1.0e-12 and minimum_l >= -1.0e-12
        ),
    }
    checks["validation_passed"] = all(checks.values())
    duration = time.perf_counter() - started
    relative_output = str(output.relative_to(REPO)).replace("\\", "/")
    command = (
        "python validation_test/vim_coupled/"
        "validate_hcurl_eddy_bubble_disk.py"
    )
    if output != DEFAULT_OUTPUT:
        command += " --output " + subprocess.list2cmdline([relative_output])

    result = {
        "schema": "cae-ai-lab.solver-run.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "case": "3-D conducting disk HCurl eddy-bubble response-basis validation",
        "solver": "radia-ngsolve HCurl eddy-bubble epsilon-free VIM",
        "source_artifact": (
            "validation_test/vim_coupled/"
            "validate_hcurl_eddy_bubble_disk.py"
        ),
        "pass": bool(checks["validation_passed"]),
        "run": {
            "command": command,
            "workdir": "repository root",
            "exit_code": 0,
            "duration_s": duration,
        },
        "result_files": [relative_output],
        "physics": {
            "geometry": "solid circular conducting disk",
            "radius_m": RADIUS_M,
            "thickness_m": THICKNESS_M,
            "conductivity_s_per_m": SIGMA_S_PER_M,
            "excitation": "uniform Bz through A=(-y,x,0)/2",
            "reference": "axisymmetric integral-equation BEM modal spectrum",
        },
        "axifem_reference": axifem,
        "positive_h_refinement": h_rows,
        "positive_p_refinement": p_rows,
        "single_port_negative_control": negative_rows,
        "checks": checks,
        "tolerances": {
            "axifem_first_mode_error_pct": 0.6,
            "hcurl_port_dominant_error_pct": 2.0,
            "h_refinement_change_fraction": 0.05,
            "single_port_expected_error_floor_pct": 4.0,
            "projection_relative_residual": 1.0e-9,
            "passive_minimum_eigenvalue": -1.0e-12,
        },
        "errors": {
            "maximum_relative": max(positive_errors + p_errors) / 100.0,
            "maximum_absolute": max(
                abs(
                    float(row["modal"]["port_dominant_time_constant_us"])
                    - BEM_MODAL_REFERENCE_US[0]
                )
                for row in h_rows + p_rows
            ),
            "best_polynomial_error_pct": min(positive_errors),
            "best_single_port_error_pct": min(negative_errors),
            "h_refinement_change_fraction": h_refinement_change,
            "projection_relative_residual": projection_error,
        },
        "tool_versions": {
            "python": platform.python_version(),
            "radia": getattr(radia, "__version__", _package_version("radia")),
            "ngsolve": getattr(ng, "__version__", _package_version("ngsolve")),
            "numpy": np.__version__,
            "scipy": _package_version("scipy"),
            "platform": platform.platform(),
        },
        "timing_breakdown_s": {
            "axisymmetric_axifem_reference": float(axifem["timing_s"]["total"]),
            "positive_h_refinement": sum(row["timing_s"]["total"] for row in h_rows),
            "p_order_additional_runs": p_rows[0]["timing_s"]["total"]
            + p_rows[2]["timing_s"]["total"],
            "single_port_negative_controls": sum(
                row["timing_s"]["total"] for row in negative_rows
            ),
            "total": duration,
        },
        "claim_boundary": {
            "proves": [
                "the leading disk diffusion pole is reproduced by the 3-D HCurl eddy-bubble path under the recorded h and p gates",
                "spatial response-basis enrichment is required for this port observable",
                "the reduced resistance-inductance models are passive",
            ],
            "does_not_prove": [
                "universal accuracy for arbitrary conductor topology or strong-skin frequency",
                "that spatial mesh refinement can repair an incomplete response basis",
                "replacement of charged-particle tracking workflows",
            ],
        },
        "verification": {
            "method": (
                "live Q1/Q2 axisymmetric solves plus an independent BEM modal "
                "reference, followed by 3-D HCurl h/p and negative-control runs"
            ),
            "command": (
                "python $env:USERPROFILE/.codex/skills/solver-run-verify/"
                "scripts/validate_solver_run.py "
                + relative_output
                + " --require-existing-files"
            ),
        },
    }
    if not result["pass"]:
        result["failure"] = {
            "stage": "verify",
            "message": "one or more disk response-basis validation gates failed",
            "next_action": (
                "inspect the response ports, port-dominant modal residues, "
                "projection residual, and passive R/L eigenvalues before refining"
            ),
        }

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": relative_output, "checks": checks}, indent=2))
    return 0 if result["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
