"""TEAM 28 3-D HCurl Eddy Bubble + analytic-tet VIM force acceptance.

This is the numerical gate left open by ``team28_hcurl_eddy_bubble.py``.  It
uses a real p=6 HCurl parent, EVRS response reduction, the epsilon-free affine
tetrahedron Newton-potential Gram, and the physical time-average Lorentz force.

The default mesh sweep is validation-class work.  Run it on an idle compute
host; ``--maxh 0.025`` is the smaller developer smoke.
"""

from __future__ import annotations

import argparse
import json
import platform
import socket
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.special import ellipe, ellipk

import radia
import radia.vim as vim


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "team28_hcurl_vim_force_summary.json"
MU0 = 4.0e-7 * np.pi
FREQUENCY_HZ = 50.0
SIGMA_AL = 3.4e7
TARGET_PHYSICAL_FORCE_N = 0.5 * 2.1925321151130186
REFERENCE_COIL_CURRENT_A = 20.0


def _circular_loop_fields(points, radius, z0, current):
    """Vectorized SI ``(A,B)`` for one z-axis circular filament."""

    points = np.asarray(points)
    x, y, z = points.T
    rho = np.hypot(x, y)
    dz = z - z0
    denominator = (radius + rho) ** 2 + dz**2
    parameter = np.clip(4.0 * radius * rho / denominator, 0.0, 1.0 - 1.0e-14)
    elliptic_k = ellipk(parameter)
    elliptic_e = ellipe(parameter)
    root = np.sqrt(denominator)
    distance = (radius - rho) ** 2 + dz**2
    common = MU0 * current / (2.0 * np.pi * root)

    off_axis = rho > 1.0e-14
    radial_b = np.zeros_like(rho)
    radial_b[off_axis] = common[off_axis] * dz[off_axis] / rho[off_axis] * (
        -elliptic_k[off_axis]
        + (radius**2 + rho[off_axis] ** 2 + dz[off_axis] ** 2)
        / distance[off_axis]
        * elliptic_e[off_axis]
    )
    axial_b = common * (
        elliptic_k
        + (radius**2 - rho**2 - dz**2) / distance * elliptic_e
    )

    elliptic_modulus = np.sqrt(parameter)
    valid_a = off_axis & (elliptic_modulus > 1.0e-14)
    a_phi = np.zeros_like(rho)
    a_phi[valid_a] = (
        (MU0 * current / np.pi)
        * np.sqrt(radius / rho[valid_a])
        * (
            (1.0 - parameter[valid_a] / 2.0) * elliptic_k[valid_a]
            - elliptic_e[valid_a]
        )
        / elliptic_modulus[valid_a]
    )

    radial = np.zeros((len(points), 3))
    azimuthal = np.zeros((len(points), 3))
    radial[off_axis, 0] = x[off_axis] / rho[off_axis]
    radial[off_axis, 1] = y[off_axis] / rho[off_axis]
    azimuthal[off_axis, 0] = -y[off_axis] / rho[off_axis]
    azimuthal[off_axis, 1] = x[off_axis] / rho[off_axis]
    magnetic_flux_density = radial_b[:, None] * radial
    magnetic_flux_density[:, 2] += axial_b
    vector_potential = a_phi[:, None] * azimuthal
    return vector_potential, magnetic_flux_density


def _rectangular_coil_fields(
    points,
    mean_radius,
    width,
    height,
    ampere_turns,
    *,
    radial_order=10,
    axial_order=14,
):
    """Gauss integrate a uniform azimuthal winding over its r-z section."""

    vector_potential = np.zeros((len(points), 3))
    magnetic_flux_density = np.zeros((len(points), 3))
    radial_nodes, radial_weights = np.polynomial.legendre.leggauss(radial_order)
    axial_nodes, axial_weights = np.polynomial.legendre.leggauss(axial_order)
    radii = mean_radius + 0.5 * width * radial_nodes
    z_values = -0.5 * height + 0.5 * height * axial_nodes
    for radius, radial_weight in zip(radii, radial_weights):
        for z0, axial_weight in zip(z_values, axial_weights):
            filament_current = ampere_turns * radial_weight * axial_weight / 4.0
            a_value, b_value = _circular_loop_fields(
                points,
                radius,
                z0,
                filament_current,
            )
            vector_potential += a_value
            magnetic_flux_density += b_value
    return vector_potential, magnetic_flux_density


def _external_coil_fields(
    points,
    coil_current=REFERENCE_COIL_CURRENT_A,
    height_offset_m=0.0,
):
    points = np.asarray(points, dtype=float).copy()
    points[:, 2] += float(height_offset_m)
    scale = float(coil_current) / REFERENCE_COIL_CURRENT_A
    a1, b1 = _rectangular_coil_fields(
        points,
        0.041,
        0.028,
        0.052,
        960.0 * REFERENCE_COIL_CURRENT_A * scale,
    )
    a2, b2 = _rectangular_coil_fields(
        points,
        0.0875,
        0.015,
        0.052,
        -576.0 * REFERENCE_COIL_CURRENT_A * scale,
    )
    return a1 + a2, b1 + b2


def _build_eddy_basis(maxh_m):
    import netgen.occ as occ
    import ngsolve as ng

    disk = occ.Cylinder(occ.Pnt(0.0, 0.0, 0.0108), occ.Z, 0.065, 0.003)
    disk.mat("Al")
    for face in disk.faces:
        face.name = "disk_air"
    mesh = ng.Mesh(occ.OCCGeometry(disk).GenerateMesh(maxh=maxh_m))
    fes = ng.HCurl(mesh, order=6, nograds=True)
    u, v = fes.TnT()
    stiffness = ng.BilinearForm(fes)
    stiffness += (ng.curl(u) * ng.curl(v) + 100.0 * u * v) * ng.dx
    metric = ng.BilinearForm(fes)
    metric += u * v * ng.dx
    training = (
        ng.CF((-ng.y, ng.x, 0.0)),
        (ng.x * ng.x + ng.y * ng.y) * ng.CF((-ng.y, ng.x, 0.0)),
        ng.z * ng.CF((-ng.y, ng.x, 0.0)),
    )
    with ng.TaskManager():
        stiffness.Assemble()
        metric.Assemble()
        ports = []
        for vector_potential in training:
            form = ng.LinearForm(fes)
            form += vector_potential * ng.curl(v) * ng.dx
            form.Assemble()
            ports.append(form.vec.FV().NumPy().copy())
        basis = vim.NgsolveEddyBubbleHCurlBasis(
            mesh,
            fes,
            stiffness,
            metric,
            np.column_stack(ports),
            steps=1,
            conductive_materials="Al",
            volume_materials="Al",
            intorder=10,
            parent_order=6,
            current_gram_rtol=1.0e-10,
        )
    return ng, mesh, fes, basis


def run_case(maxh_m, outer_quad, export_model=None):
    started = time.perf_counter()
    ng, mesh, fes, basis = _build_eddy_basis(maxh_m)
    with ng.TaskManager():
        interaction = basis.tet_volume_interaction(
            mesh,
            fes,
            degree=5,
            projection_quad=7,
            outer_quad=outer_quad,
            projection_tolerance=1.0e-10,
            materials="Al",
        )
    external_a, external_b = _external_coil_fields(basis.current_basis.points)
    unit_external_a = external_a / REFERENCE_COIL_CURRENT_A
    unit_external_b = external_b / REFERENCE_COIL_CURRENT_A
    system = basis.assemble_vim(
        sigma=SIGMA_AL,
        interaction=interaction,
    )
    rhs = vim.ExternalVectorPotentialRHS(basis.current_basis, unit_external_a)
    cln_model = vim.HCurlEddyCLNFromVIM(system, rhs)
    s = 2.0j * np.pi * FREQUENCY_HZ
    coefficients = cln_model.solve_vector_potential_drive(s, REFERENCE_COIL_CURRENT_A)
    current = np.einsum("a,aik->ik", coefficients, basis.current_basis.modes)
    force = 0.5 * np.sum(
        basis.current_basis.weights[:, None]
        * np.real(np.cross(current, np.conj(external_b))),
        axis=0,
    )
    force_operator = np.transpose(
        np.sum(
            basis.current_basis.weights[None, :, None]
            * np.cross(
                basis.current_basis.modes,
                np.conj(unit_external_b)[None, :, :],
            ),
            axis=1,
        ),
        (1, 0),
    )[:, :, None]
    if export_model is not None:
        vim.ExportHCurlEddyCLNJSON(
            cln_model,
            export_model,
            force_operator=force_operator,
            metadata={
                "frequency_hz": FREQUENCY_HZ,
                "coil_current_reference_A": REFERENCE_COIL_CURRENT_A,
                "conductivity_S_per_m": SIGMA_AL,
                "force_operator_units": "N per reduced coefficient and ampere port current",
                "parent_space": "HCurl",
                "parent_order": 6,
            },
        )
    relative_error = abs(abs(force[2]) - TARGET_PHYSICAL_FORCE_N) / TARGET_PHYSICAL_FORCE_N
    transverse_ratio = float(np.linalg.norm(force[:2]) / max(abs(force[2]), np.finfo(float).tiny))
    return {
        "maxh_m": float(maxh_m),
        "outer_quad": int(outer_quad),
        "mesh_elements": int(mesh.ne),
        "parent_ndof": int(fes.ndof),
        "evrs_rank": int(basis.rank),
        "runtime_reduction_ratio": float(basis.rank / fes.ndof),
        "current_sample_count": int(basis.current_basis.n_samples),
        "physical_force_N": force.tolist(),
        "target_physical_force_magnitude_N": TARGET_PHYSICAL_FORCE_N,
        "magnitude_relative_error": float(relative_error),
        "transverse_force_ratio": transverse_ratio,
        "interaction": interaction.diagnostics(),
        "cln_handoff": cln_model.diagnostics(),
        "eddy_bubble": basis.eddy_bubbling.diagnostics(),
        "elapsed_seconds": time.perf_counter() - started,
    }


def _height_offsets_mm_from_reference():
    sweep_file = (
        HERE.parent.parent
        / "docs"
        / "maglev"
        / "demos"
        / "team28"
        / "team28_cln_sweep_results.json"
    )
    with sweep_file.open("r", encoding="utf-8") as stream:
        return [float(value) for value in json.load(stream)["dZ_mm"]]


def run_height_sweep(
    height_offsets_mm,
    *,
    maxh_m=0.025,
    outer_quad=4,
    export_family=None,
):
    """Build one p=6 basis and export a common-coordinate height family.

    The disk mesh and reduced state coordinates stay fixed at the reference
    position. Only the incident coil A/B fields and force operator vary with
    height, so stateful moving coupling does not hide a basis transfer.
    """

    started = time.perf_counter()
    ng, mesh, fes, basis = _build_eddy_basis(maxh_m)
    with ng.TaskManager():
        interaction = basis.tet_volume_interaction(
            mesh,
            fes,
            degree=5,
            projection_quad=7,
            outer_quad=outer_quad,
            projection_tolerance=1.0e-10,
            materials="Al",
        )
    system = basis.assemble_vim(sigma=SIGMA_AL, interaction=interaction)
    snapshots = []
    physical_forces = []
    for offset_mm in height_offsets_mm:
        offset_m = float(offset_mm) * 1.0e-3
        external_a, external_b = _external_coil_fields(
            basis.current_basis.points,
            coil_current=REFERENCE_COIL_CURRENT_A,
            height_offset_m=offset_m,
        )
        unit_external_a = external_a / REFERENCE_COIL_CURRENT_A
        unit_external_b = external_b / REFERENCE_COIL_CURRENT_A
        rhs = vim.ExternalVectorPotentialRHS(basis.current_basis, unit_external_a)
        model = vim.HCurlEddyCLNFromVIM(system, rhs)
        s = 2.0j * np.pi * FREQUENCY_HZ
        coefficients = model.solve_vector_potential_drive(
            s,
            REFERENCE_COIL_CURRENT_A,
        )
        current = np.einsum("a,aik->ik", coefficients, basis.current_basis.modes)
        force = 0.5 * np.sum(
            basis.current_basis.weights[:, None]
            * np.real(np.cross(current, np.conj(external_b))),
            axis=0,
        )
        force_operator = np.transpose(
            np.sum(
                basis.current_basis.weights[None, :, None]
                * np.cross(
                    basis.current_basis.modes,
                    np.conj(unit_external_b)[None, :, :],
                ),
                axis=1,
            ),
            (1, 0),
        )[:, :, None]
        snapshots.append(
            {
                "height_m": offset_m,
                "model": model,
                "force_operator": force_operator,
                "metadata": {"height_offset_m": offset_m},
            }
        )
        physical_forces.append(force.tolist())
    if export_family is not None:
        vim.ExportHCurlEddyCLNFamilyJSON(
            snapshots,
            export_family,
            metadata={
                "frequency_hz": FREQUENCY_HZ,
                "coil_current_reference_A": REFERENCE_COIL_CURRENT_A,
                "conductivity_S_per_m": SIGMA_AL,
                "parent_space": "HCurl",
                "parent_order": 6,
                "reference_disk_bottom_m": 0.0108,
                "height_coordinate": "offset from the reference disk position",
                "state_basis": "common p=6 local disk basis",
            },
        )
    return {
        "schema": "radia.team28.hcurl-vim-family.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "runtime": {
            "hostname": socket.gethostname(),
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "radia_version": getattr(radia, "__version__", "unknown"),
            "ngsolve_version": getattr(__import__("ngsolve"), "__version__", "unknown"),
            "elapsed_seconds": time.perf_counter() - started,
        },
        "problem": {
            "frequency_hz": FREQUENCY_HZ,
            "conductivity_S_per_m": SIGMA_AL,
            "parent_space": "HCurl",
            "parent_order": 6,
            "common_state_basis": True,
        },
        "height_offsets_mm": [float(value) for value in height_offsets_mm],
        "physical_force_N": physical_forces,
        "snapshot_count": len(snapshots),
        "export_file": Path(export_family).name if export_family is not None else None,
    }


def run(maxh_values, outer_quad=4, outer_check=None, export_model=None):
    cases = [
        run_case(
            value,
            outer_quad,
            export_model=export_model if index == 0 else None,
        )
        for index, value in enumerate(maxh_values)
    ]
    outer_case = None
    if outer_check is not None:
        outer_case = run_case(maxh_values[0], outer_check)
    force_values = np.asarray([abs(case["physical_force_N"][2]) for case in cases])
    max_force_error = max(case["magnitude_relative_error"] for case in cases)
    outer_relative_change = None
    if outer_case is not None:
        outer_relative_change = float(
            abs(abs(outer_case["physical_force_N"][2]) - force_values[0])
            / force_values[0]
        )
    checks = {
        "all_projection_residuals_below_1e-10": all(
            case["interaction"]["projection_relative_residual"] < 1.0e-10
            for case in cases
        ),
        "all_inductance_blocks_positive": all(
            case["interaction"]["minimum_eigenvalue_H"] > 0.0
            for case in cases
        ),
        "all_force_errors_below_one_percent": max_force_error < 0.01,
        "all_transverse_force_ratios_below_half_percent": all(
            case["transverse_force_ratio"] < 0.005 for case in cases
        ),
        "all_runtime_reductions_below_one_percent": all(
            case["runtime_reduction_ratio"] < 0.01 for case in cases
        ),
        "no_kernel_epsilon": all(
            case["interaction"]["kernel_epsilon_m"] is None for case in cases
        ),
        "all_hcurl_to_cln_handoffs_passive": all(
            case["cln_handoff"]["passive"]
            and case["cln_handoff"]["state_order"] == case["evrs_rank"]
            and case["cln_handoff"]["port_count"] == 1
            for case in cases
        ),
    }
    if outer_relative_change is not None:
        checks["outer_quadrature_change_below_0p1_percent"] = outer_relative_change < 0.001
    checks = {name: bool(passed) for name, passed in checks.items()}
    return {
        "schema": "radia.team28.hcurl-vim-force.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "runtime": {
            "hostname": socket.gethostname(),
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "radia_version": getattr(radia, "__version__", "unknown"),
            "ngsolve_version": getattr(__import__("ngsolve"), "__version__", "unknown"),
        },
        "problem": {
            "frequency_hz": FREQUENCY_HZ,
            "conductivity_S_per_m": SIGMA_AL,
            "disk_radius_m": 0.065,
            "disk_thickness_m": 0.003,
            "disk_bottom_m": 0.0108,
            "parent_space": "HCurl",
            "parent_order": 6,
            "force_convention": "positive z is upward; compare magnitude to stored TEAM physical force",
            "coil_cross_section_rule": "Gauss-Legendre uniform winding density",
        },
        "cases": cases,
        "outer_quadrature_check": outer_case,
        "outer_quadrature_relative_force_change": outer_relative_change,
        "force_magnitude_range_N": [float(force_values.min()), float(force_values.max())],
        "maximum_force_relative_error": float(max_force_error),
        "checks": checks,
        "hcurl_vim_force_acceptance_complete": bool(all(checks.values())),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxh", type=float, nargs="+", default=[0.025, 0.020, 0.015])
    parser.add_argument("--outer-quad", type=int, default=4)
    parser.add_argument("--outer-check", type=int)
    parser.add_argument(
        "--export-model",
        type=Path,
        help="write the first p=6 HCurl-VIM case as a MATLAB exchange JSON",
    )
    parser.add_argument(
        "--export-family",
        type=Path,
        help="write a common-basis height-indexed p=6 HCurl-VIM family JSON",
    )
    parser.add_argument(
        "--height-offsets-mm",
        type=str,
        nargs="+",
        help="height offsets for --export-family (space or comma separated; default is the 25-point TEAM 28 sweep)",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.export_family is not None:
        offsets = (
            [
                float(part)
                for token in args.height_offsets_mm
                for part in token.split(",")
                if part.strip()
            ]
            if args.height_offsets_mm is not None
            else _height_offsets_mm_from_reference()
        )
        result = run_height_sweep(
            offsets,
            maxh_m=args.maxh[0],
            outer_quad=args.outer_quad,
            export_family=args.export_family,
        )
    else:
        result = run(
            args.maxh,
            outer_quad=args.outer_quad,
            outer_check=args.outer_check,
            export_model=args.export_model,
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    if not result.get("hcurl_vim_force_acceptance_complete", True):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
