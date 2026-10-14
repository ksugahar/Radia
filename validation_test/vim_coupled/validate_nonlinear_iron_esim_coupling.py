"""Validate a conducting ferromagnet with HDiv-MMM and local ESIM-SIBC.

The same regenerated TET body carries the magnetic and conductive labels.
HDiv-MMM represents a fixed low-field bulk permeability, while HCurl
eddy-bubble plus surface-Omega carries volume/bridge/skin currents.  A local
one-dimensional nonlinear B-H cell updates the SIBC Gram around the complete
HDiv/HCurl mixed solve.

This is a heavy validation runner, not a unit test.  It deliberately separates
the established coupled local-ESIM contract from the still-unimplemented
simultaneous update of an ordinary bulk nonlinear HDiv constitutive operator.
No generated mesh is tracked.
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

import netgen.occ as occ
import ngsolve as ng
import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "src"
DEFAULT_OUTPUT = HERE / "results_nonlinear_iron_esim_coupling.json"
COMPANION_VOLUMETRIC = HERE / "results_magnetic_conductor_disk_adjudication.json"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import radia  # noqa: E402
from radia import vim  # noqa: E402


MU0 = 4.0e-7 * np.pi
SIGMA_S_PER_M = 2.0e6
FREQUENCY_HZ = 50_000.0
FIELD_AMPLITUDES_A_PER_M = (50.0, 1_000.0, 5_000.0)
BH_CURVE = np.array(
    [
        [0.0, 0.0],
        [100.0, 0.1],
        [500.0, 0.45],
        [2_000.0, 1.2],
        [10_000.0, 1.6],
        [50_000.0, 1.8],
    ]
)
LOW_FIELD_MU_R = float(BH_CURVE[1, 1] / (MU0 * BH_CURVE[1, 0]))


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
        SRC / "radia" / "vim" / "_eddy_hybrid.py",
        SRC / "radia" / "vim" / "_vim.py",
        SRC / "radia" / "esim_cell_problem.py",
    )
    return {
        path.relative_to(REPO).as_posix(): _sha256(path)
        for path in paths
    }


def _complex(value) -> dict[str, float]:
    number = complex(value)
    return {"real": float(number.real), "imag": float(number.imag)}


def _make_mesh(maxh: float) -> ng.Mesh:
    body = occ.Box(occ.Pnt(0.0, 0.0, 0.0), occ.Pnt(1.0, 1.0, 1.0))
    body.mat("iron")
    for face in body.faces:
        face.name = "skin"
    return ng.Mesh(occ.OCCGeometry(body).GenerateMesh(maxh=maxh))


def _assemble(maxh: float):
    mesh_started = time.perf_counter()
    mesh = _make_mesh(maxh)
    mesh_s = time.perf_counter() - mesh_started

    hcurl = ng.HCurl(mesh, order=2, nograds=True)
    trial, test = hcurl.TnT()
    stiffness = ng.BilinearForm(hcurl)
    stiffness += (
        ng.curl(trial) * ng.curl(test) + 0.05 * trial * test
    ) * ng.dx
    mass = ng.BilinearForm(hcurl)
    mass += trial * test * ng.dx
    port = ng.LinearForm(hcurl)
    port += ng.CF((0.0, 0.0, ng.y)) * test * ng.dx
    assembly_started = time.perf_counter()
    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        port.Assemble()

    def unit_vector_potential(points):
        points = np.asarray(points)
        return np.column_stack(
            (
                np.zeros(points.shape[0]),
                np.zeros(points.shape[0]),
                MU0 * points[:, 1],
            )
        )

    applicability = vim.EddySIBCApplicability(
        frequency_hz=FREQUENCY_HZ,
        sigma=SIGMA_S_PER_M,
        characteristic_thickness_m=1.0,
        mu=MU0 * LOW_FIELD_MU_R,
        characteristic_curvature_radius_m=0.5,
    )
    material = vim.SharedMeshMaterialModel(
        mesh=mesh,
        magnetic_regions="iron",
        conductive_regions="iron",
        mu=MU0 * LOW_FIELD_MU_R,
        sigma=SIGMA_S_PER_M,
        sibc="local-esim",
    )
    surface_modes = (
        ng.CF((1.0, 0.0, 0.0)),
        ng.CF((0.0, 1.0, 0.0)),
        ng.CF((0.0, 0.0, 1.0)),
    )
    with ng.TaskManager():
        mixed = vim.NgsolveBDMEddyBubbleVIM(
            mesh,
            hcurl,
            stiffness,
            mass,
            (port,),
            surface_modes,
            hdiv_order=1,
            mu_r=LOW_FIELD_MU_R,
            external_fields=(ng.CF((1.0, 0.0, 0.0)),),
            external_names=("unit_uniform_Hx",),
            hdiv_max_modes=1,
            hdiv_solve_tol=1.0e-9,
            magnetic_materials="iron",
            steps=2,
            sigma=SIGMA_S_PER_M,
            conductive_materials="iron",
            surface_boundaries="skin",
            intorder=2,
            kernel_epsilon=0.05,
            response_backend="dense",
            sibc_applicability=applicability,
            port_vector_potentials=(unit_vector_potential,),
            material_model=material,
            coupling_kernel_epsilon=0.05,
        )
    assembly_s = time.perf_counter() - assembly_started
    return mesh, hcurl, mixed, applicability, material, mesh_s, assembly_s


def _solve_amplitude(mixed, model, amplitude: float) -> dict[str, object]:
    started = time.perf_counter()
    s = 1j * 2.0 * np.pi * FREQUENCY_HZ
    magnetic_rhs = amplitude * mixed.magnetic_rhs
    eddy_rhs = -s * amplitude * mixed.eddy_rhs
    solve_options = {
        "magnetic_rhs": magnetic_rhs,
        "eddy_rhs": eddy_rhs,
        "mixed_galerkin_keep_blocks": ("volume1", "surface"),
        "mixed_galerkin_eliminate_blocks": "volume",
    }
    linear_impedance = model.linear_impedance(FREQUENCY_HZ)
    linear = mixed.solve_frequency(
        FREQUENCY_HZ,
        surface_impedance=linear_impedance,
        **solve_options,
    )
    nonlinear = mixed.solve_frequency_local_esim(
        model,
        FREQUENCY_HZ,
        outer_tolerance=5.0e-3,
        outer_max_iterations=20,
        outer_relaxation=0.7,
        **solve_options,
    )
    fixed_gram = mixed.solve_frequency(
        FREQUENCY_HZ,
        surface_impedance=nonlinear.surface_impedance,
        **solve_options,
    )
    local_solution = nonlinear.mixed_solution
    fixed_denominator = max(
        float(np.linalg.norm(local_solution.reduced_solution)),
        np.finfo(float).tiny,
    )
    linear_denominator = max(
        float(np.linalg.norm(linear.reduced_solution)),
        np.finfo(float).tiny,
    )
    port_denominator = max(
        float(np.linalg.norm(linear.port_response)),
        np.finfo(float).tiny,
    )
    samples = nonlinear.surface_impedance.sample_values
    return {
        "applied_H_A_per_m": float(amplitude),
        "linear_surface_impedance_ohm": _complex(linear_impedance),
        "local_surface_impedance": {
            "real_min_ohm": float(np.min(samples.real)),
            "real_mean_ohm": float(np.mean(samples.real)),
            "real_max_ohm": float(np.max(samples.real)),
            "imag_min_ohm": float(np.min(samples.imag)),
            "imag_mean_ohm": float(np.mean(samples.imag)),
            "imag_max_ohm": float(np.max(samples.imag)),
            "passive": bool(nonlinear.surface_impedance.diagnostics()["passive"]),
        },
        "local_esim": {
            "converged": bool(nonlinear.converged),
            "iterations": int(nonlinear.iterations),
            "relative_impedance_change": float(
                nonlinear.relative_impedance_change
            ),
            "h_min_A_per_m": float(
                np.min(nonlinear.tangential_field_magnitude)
            ),
            "h_max_A_per_m": float(
                np.max(nonlinear.tangential_field_magnitude)
            ),
            "cell_solve_count": int(
                sum(row["cell_solve_count"] for row in nonlinear.history)
            ),
        },
        "mixed": {
            "residual_relative_norm": float(
                local_solution.residual_relative_norm
            ),
            "residual_backward_error": float(
                local_solution.residual_backward_error
            ),
            "joule_loss_W": float(local_solution.average_joule_loss[0]),
            "linear_sibc_joule_loss_W": float(linear.average_joule_loss[0]),
            "local_vs_linear_solution_relative_difference": float(
                np.linalg.norm(
                    local_solution.reduced_solution - linear.reduced_solution
                )
                / linear_denominator
            ),
            "local_vs_linear_port_relative_difference": float(
                np.linalg.norm(local_solution.port_response - linear.port_response)
                / port_denominator
            ),
            "fixed_gram_replay_relative_difference": float(
                np.linalg.norm(
                    local_solution.reduced_solution - fixed_gram.reduced_solution
                )
                / fixed_denominator
            ),
            "magnetization_norm": float(
                np.linalg.norm(local_solution.sampled_magnetization)
            ),
            "eddy_current_norm": float(
                sum(
                    np.linalg.norm(values)
                    for values in local_solution.sampled_eddy_currents
                )
            ),
        },
        "wall_s": time.perf_counter() - started,
    }


def run(maxh: float = 0.65) -> dict[str, object]:
    started = time.perf_counter()
    source_start = _source_fingerprints()
    head_start = _git_head()
    (
        mesh,
        hcurl,
        mixed,
        applicability,
        material,
        mesh_s,
        assembly_s,
    ) = _assemble(maxh)
    model = vim.LocalESIMSurfaceModel(
        bh_curve=BH_CURVE,
        sigma=SIGMA_S_PER_M,
        bins=4,
        n_nodes=40,
        cell_tolerance=1.0e-4,
        cell_max_iterations=80,
    )
    rows = [
        _solve_amplitude(mixed, model, amplitude)
        for amplitude in FIELD_AMPLITUDES_A_PER_M
    ]
    source_end = _source_fingerprints()
    real_means = [row["local_surface_impedance"]["real_mean_ohm"] for row in rows]
    nonlinear_differences = [
        row["mixed"]["local_vs_linear_solution_relative_difference"]
        for row in rows
    ]
    companion = json.loads(COMPANION_VOLUMETRIC.read_text(encoding="utf-8"))
    checks = {
        "source_fingerprints_stable_during_run": source_start == source_end,
        "same_region_is_magnetic_and_conductive": (
            material.magnetic_regions == ("iron",)
            and material.conductive_regions == ("iron",)
        ),
        "high_frequency_route_selects_sibc": applicability.sibc_applicable,
        "all_local_esim_outer_loops_converged": all(
            row["local_esim"]["converged"] for row in rows
        ),
        "all_local_surface_grams_are_passive": all(
            row["local_surface_impedance"]["passive"] for row in rows
        ),
        "all_mixed_residuals_below_1e-10": max(
            row["mixed"]["residual_relative_norm"] for row in rows
        )
        < 1.0e-10,
        "all_joule_losses_are_positive": min(
            row["mixed"]["joule_loss_W"] for row in rows
        )
        > 0.0,
        "converged_gram_replay_matches_mixed_solution": max(
            row["mixed"]["fixed_gram_replay_relative_difference"] for row in rows
        )
        < 1.0e-12,
        "surface_impedance_changes_with_field_amplitude": (
            max(real_means) - min(real_means)
        )
        / max(real_means)
        > 0.01,
        "high_field_reduces_mean_surface_resistance": real_means[-1] < real_means[0],
        "nonlinear_departure_grows_from_low_to_high_field": (
            nonlinear_differences[-1] > nonlinear_differences[0]
        ),
        "companion_volumetric_magnetic_conductor_gate_passes": bool(
            companion.get("checks")
            and all(companion["checks"].values())
        ),
    }
    total_s = time.perf_counter() - started
    return {
        "schema": "radia.validation.nonlinear-iron-local-esim-coupling.v1",
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
            "geometry": "solid conducting ferromagnetic cube",
            "mesh": "regenerated in-memory TET; no mesh tracked",
            "material_region": "iron",
            "magnetic_region": "iron",
            "conductive_region": "iron",
            "bulk_magnetic_law": "fixed low-field scalar permeability",
            "bulk_mu_r": LOW_FIELD_MU_R,
            "surface_magnetic_law": "monotone nonlinear B-H local ESIM cell",
            "bh_curve_H_A_per_m_B_T": BH_CURVE.tolist(),
            "sigma_S_per_m": SIGMA_S_PER_M,
            "frequency_Hz": FREQUENCY_HZ,
            "field_amplitudes_A_per_m": list(FIELD_AMPLITUDES_A_PER_M),
        },
        "discretization": {
            "tetrahedra": int(mesh.ne),
            "hcurl_order": 2,
            "hcurl_dofs": int(hcurl.ndof),
            "hdiv_family": "BDM",
            "hdiv_order": 1,
            "hdiv_modes": int(mixed.n_hdiv_modes),
            "eddy_modes": int(mixed.n_hcurl_vim_modes),
            "eddy_block_roles": dict(mixed.eddy_block_roles),
            "sibc_applicability": applicability.diagnostics(),
        },
        "amplitude_ladder": rows,
        "checks": checks,
        "pass": all(checks.values()),
        "timing_s": {
            "mesh": mesh_s,
            "assembly_and_reduction": assembly_s,
            "amplitude_solves": float(sum(row["wall_s"] for row in rows)),
            "total": total_s,
        },
        "companion_evidence": {
            "artifact": COMPANION_VOLUMETRIC.relative_to(REPO).as_posix(),
            "purpose": (
                "axisymmetric Q2 and full 3-D HCurl volumetric accuracy at "
                "frequencies where skin penetration is resolved"
            ),
        },
        "claim_boundary": {
            "established": (
                "a single TET iron region can carry fixed bulk HDiv-MMM "
                "magnetization, HCurl eddy-bubble currents, and a passive "
                "field-dependent local ESIM surface law in one converged "
                "mixed harmonic solve"
            ),
            "not_established": (
                "a simultaneous ordinary bulk nonlinear B-H update, hysteretic "
                "or rotational skin state, multidimensional corner cells, or "
                "universal accuracy and performance superiority"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxh", type=float, default=0.65)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    artifact = run(args.maxh)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"pass": artifact["pass"], "checks": artifact["checks"]}, indent=2))
    if not artifact["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
