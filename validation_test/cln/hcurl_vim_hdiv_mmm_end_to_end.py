"""End-to-end HCurl-VIM / HDiv-MMM production-path smoke.

The durable path exercised here is:

    NGSolve HCurl(p) parent matrices and ports
      -> EVRS eddy bubbling
      -> topology-aware bulk / bridge-cycle / SIBC VIM
      -> protected physical + multipole-trained HDiv response-POD basis
      -> projected Radia HDiv material and demagnetizing operators
      -> full coupled mixed-Galerkin elimination of bulk eddy bubbles
      -> parent-T, parent-HDiv, J, K, and M reconstruction.

This is a correctness validation.  Solver-heavy runs belong on a designated
compute host; timing is secondary to the stored algebraic and field checks.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
SRC = REPO_ROOT / "src"
DEFAULT_OUTPUT = HERE / "hcurl_vim_hdiv_mmm_end_to_end_smoke.json"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import radia.vim as vim


def _parse_floats(text: str) -> list[float]:
    values = [float(part.strip()) for part in text.split(",") if part.strip()]
    if not values or any(value <= 0.0 for value in values):
        raise argparse.ArgumentTypeError("expected positive comma-separated values")
    return values


def _complex_matrix_parts(matrix) -> dict[str, object]:
    values = np.asarray(matrix)
    return {
        "shape": list(values.shape),
        "real": values.real.tolist(),
        "imag": values.imag.tolist(),
    }


def _make_mesh(
    maxh: float,
    geometry: str = "notched-box",
    axial_length: float = 0.2,
    curvature_safety: float = 2.0,
):
    import ngsolve as ng
    import netgen.occ as occ

    if geometry == "annular-motor":
        axis = occ.Dir(0, 0, 1)
        base = occ.Pnt(0, 0, -0.5 * axial_length)
        rotor = (
            occ.Cylinder(base, axis, r=0.8, h=axial_length)
            - occ.Cylinder(base, axis, r=0.2, h=axial_length)
        )
        stator = (
            occ.Cylinder(base, axis, r=1.4, h=axial_length)
            - occ.Cylinder(base, axis, r=1.0, h=axial_length)
        )
        rotor.mat("cond")
        stator.mat("cond")
        conductor = occ.Glue([rotor, stator])
    elif geometry == "notched-box":
        outer = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
        notch = occ.Box(occ.Pnt(0.45, 0.45, -0.1), occ.Pnt(1.1, 1.1, 1.1))
        conductor = outer - notch
    else:
        raise ValueError(f"unsupported geometry {geometry!r}")
    conductor.mat("cond")
    for face in conductor.faces:
        face.name = "skin"
    return ng.Mesh(
        occ.OCCGeometry(conductor).GenerateMesh(
            maxh=maxh,
            curvaturesafety=curvature_safety,
        )
    )


def _assemble_hcurl_parent(mesh, order: int):
    import ngsolve as ng

    fes = ng.HCurl(mesh, order=order, nograds=True)
    u, v = fes.TnT()
    stiffness = ng.BilinearForm(fes)
    stiffness += ng.curl(u) * ng.curl(v) * ng.dx + 0.05 * u * v * ng.dx
    mass = ng.BilinearForm(fes)
    mass += u * v * ng.dx

    ports = []
    for vector_potential in (
        ng.CoefficientFunction((-ng.y, ng.x, 0.0)),
        ng.CoefficientFunction((0.0, -ng.z, ng.y)),
    ):
        port = ng.LinearForm(fes)
        port += vector_potential * v * ng.dx
        ports.append(port)

    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        for port in ports:
            port.Assemble()
    return fes, stiffness, mass, ports


def _hdiv_mmm_reduction(
    mesh,
    hdiv_order: int,
    hdiv_family: str,
    intorder: int,
    mu_r: float,
    demag_eps: float,
    demag_eta: float,
    multipole_degree: int,
    hdiv_pod_rtol: float,
    hdiv_max_modes: int | None,
    hdiv_solve_tol: float,
):
    import ngsolve as ng

    fes = ng.HDiv(mesh, order=hdiv_order, RT=hdiv_family == "rt")
    multipole_ports = vim.NgsolveHDivRegularSolidHarmonicPorts(
        mesh,
        max_degree=multipole_degree,
    )
    with ng.TaskManager():
        return vim.NgsolveHDivMMMResponseReduction(
            mesh,
            fes,
            mu_r=mu_r,
            external_fields=(
                ng.CoefficientFunction((1.0, 0.0, 0.0)),
                ng.CoefficientFunction((0.0, 1.0, 0.0)),
            ),
            external_names=("H_stator_x", "H_stator_y"),
            training_fields=multipole_ports,
            pod_rtol=hdiv_pod_rtol,
            max_modes=hdiv_max_modes,
            solve_tol=hdiv_solve_tol,
            intorder=intorder,
            materials="cond",
            demag_eps=demag_eps,
            demag_eta=demag_eta,
            parent_family=hdiv_family,
            parent_order=hdiv_order,
        )


def _port_vector_potentials():
    def around_z(points):
        points = np.asarray(points)
        return np.column_stack(
            (-points[:, 1], points[:, 0], np.zeros(points.shape[0]))
        )

    def around_x(points):
        points = np.asarray(points)
        return np.column_stack(
            (np.zeros(points.shape[0]), -points[:, 2], points[:, 1])
        )

    return around_z, around_x


def _surface_grad_modes():
    import ngsolve as ng

    return (
        ng.CoefficientFunction((1.0, 0.0, 0.0)),
        ng.CoefficientFunction((0.0, 1.0, 0.0)),
        ng.CoefficientFunction((0.0, 0.0, 1.0)),
    )


def run(args: argparse.Namespace) -> dict[str, object]:
    import ngsolve as ng

    started = time.perf_counter()
    mesh = _make_mesh(
        args.maxh,
        args.geometry,
        args.axial_length,
        args.curvature_safety,
    )
    hcurl_fes, stiffness, mass, ports = _assemble_hcurl_parent(mesh, args.order)
    hdiv_reduction = None
    if args.hdiv_family == "rt":
        hdiv_reduction = _hdiv_mmm_reduction(
            mesh,
            args.hdiv_order,
            args.hdiv_family,
            args.hdiv_intorder,
            args.mu_r,
            args.demag_eps,
            args.demag_eta,
            args.multipole_degree,
            args.hdiv_pod_rtol,
            args.hdiv_max_modes,
            args.hdiv_solve_tol,
        )
    material = vim.SharedMeshMaterialModel(
        mesh=mesh,
        magnetic_regions="cond",
        conductive_regions="cond",
        mu=vim.MU0,
        sigma=args.sigma,
        sibc="half-space",
    )
    ledger = vim.EddyParentOrderLedger(
        bulk_degree=args.bulk_degree,
        bridge_trace_degree=args.bridge_trace_degree,
        surface_current_degree=args.surface_current_degree,
    )
    if args.hdiv_family == "bdm":
        multipole_ports = vim.NgsolveHDivRegularSolidHarmonicPorts(
            mesh,
            max_degree=args.multipole_degree,
        )
        mixed = vim.NgsolveBDMEddyBubbleVIM(
            mesh,
            hcurl_fes,
            stiffness,
            mass,
            ports,
            _surface_grad_modes(),
            hdiv_order=args.hdiv_order,
            mu_r=args.mu_r,
            external_fields=(
                ng.CoefficientFunction((1.0, 0.0, 0.0)),
                ng.CoefficientFunction((0.0, 1.0, 0.0)),
            ),
            external_names=("H_stator_x", "H_stator_y"),
            training_fields=multipole_ports,
            hdiv_pod_rtol=args.hdiv_pod_rtol,
            hdiv_max_modes=args.hdiv_max_modes,
            hdiv_solve_tol=args.hdiv_solve_tol,
            hdiv_intorder=args.hdiv_intorder,
            magnetic_materials="cond",
            demag_eps=args.demag_eps,
            demag_eta=args.demag_eta,
            steps=args.steps,
            sigma=args.sigma,
            conductive_materials="cond",
            surface_boundaries="skin",
            intorder=args.intorder,
            kernel_epsilon=args.kernel_epsilon,
            response_backend=args.response_backend,
            rtol=args.hcurl_rtol,
            current_gram_rtol=args.current_gram_rtol,
            parent_order_ledger=ledger,
            port_vector_potentials=_port_vector_potentials(),
            material_model=material,
            coupling_kernel_epsilon=args.coupling_kernel_epsilon,
        )
        hdiv_reduction = mixed.hdiv_reduction
    else:
        mixed = vim.NgsolveHCurlVIMHDivMMM(
            mesh,
            hcurl_fes,
            stiffness,
            mass,
            ports,
            _surface_grad_modes(),
            hdiv_reduction,
            steps=args.steps,
            sigma=args.sigma,
            conductive_materials="cond",
            surface_boundaries="skin",
            intorder=args.intorder,
            kernel_epsilon=args.kernel_epsilon,
            response_backend=args.response_backend,
            rtol=args.hcurl_rtol,
            current_gram_rtol=args.current_gram_rtol,
            parent_order_ledger=ledger,
            port_vector_potentials=_port_vector_potentials(),
            material_model=material,
            coupling_kernel_epsilon=args.coupling_kernel_epsilon,
        )
    assert hdiv_reduction is not None

    rows = []
    for frequency in args.frequencies:
        direct_solution = mixed.solve_frequency(
            frequency,
            solver=args.mixed_solver,
        )
        solution = mixed.solve_frequency_eddy_bubbled(frequency)
        operator = mixed.mixed_operator(
            None,
            solution.s,
            surface_impedance=solution.surface_impedance,
        )
        diagnostic_operator = (
            operator.to_dense()
            if callable(getattr(operator, "to_dense", None))
            else operator
        )
        diagnostics = solution.diagnostics()
        eddy_probe = solution.eddy_flux_density(
            np.array([[2.0, 2.0, 2.0]]),
            kernel_epsilon=args.coupling_kernel_epsilon,
        )
        diagnostics.update(
            {
                "mixed_operator_condition": float(
                    np.linalg.cond(diagnostic_operator)
                ),
                "direct_solution_relative_error": float(
                    np.linalg.norm(
                        solution.reduced_solution
                        - direct_solution.reduced_solution
                    )
                    / max(
                        float(np.linalg.norm(direct_solution.reduced_solution)),
                        1.0e-300,
                    )
                ),
                "direct_port_response_relative_error": float(
                    np.linalg.norm(
                        solution.port_response - direct_solution.port_response
                    )
                    / max(
                        float(np.linalg.norm(direct_solution.port_response)),
                        1.0e-300,
                    )
                ),
                "direct_joule_loss_relative_error": float(
                    np.linalg.norm(
                        solution.average_joule_loss
                        - direct_solution.average_joule_loss
                    )
                    / max(
                        float(np.linalg.norm(direct_solution.average_joule_loss)),
                        1.0e-300,
                    )
                ),
                "parent_t_coefficient_norm": float(
                    np.linalg.norm(solution.parent_t_coefficients)
                ),
                "parent_hdiv_coefficient_norm": float(
                    np.linalg.norm(solution.parent_magnetization_coefficients)
                ),
                "eddy_flux_density_probe_T": _complex_matrix_parts(
                    eddy_probe.reshape(eddy_probe.shape[0], -1)
                ),
                "eddy_flux_density_probe_norm_T": float(np.linalg.norm(eddy_probe)),
                "port_response": _complex_matrix_parts(solution.port_response),
            }
        )
        rows.append(diagnostics)

    motor_angle_rows = []
    if args.motor_angle_samples:
        magnetic_ports = np.asarray(mixed.magnetic_rhs)
        eddy_ports = np.asarray(mixed.eddy_rhs)
        if magnetic_ports.ndim != 2 or magnetic_ports.shape[1] < 2:
            raise RuntimeError("motor angle sweep requires two magnetic source ports")
        if eddy_ports.ndim != 2 or eddy_ports.shape[1] < 2:
            raise RuntimeError("motor angle sweep requires two eddy source ports")
        angles = np.arange(args.motor_angle_samples, dtype=float)
        angles *= 2.0 * np.pi / args.motor_angle_samples
        coenergy = []
        angle_solutions = []
        for angle in angles:
            coefficients = np.array(
                [[1.0 + args.rotor_amplitude * np.cos(angle)],
                 [args.rotor_amplitude * np.sin(angle)]],
                dtype=float,
            )
            solution = mixed.solve_frequency_eddy_bubbled(
                args.frequencies[0],
                magnetic_rhs=magnetic_ports[:, :2] @ coefficients,
                eddy_rhs=eddy_ports[:, :2] @ coefficients,
            )
            value = 0.5 * float(
                np.real(
                    np.vdot(
                        solution.reduced_rhs[:, 0],
                        solution.reduced_solution[:, 0],
                    )
                )
            )
            coenergy.append(value)
            angle_solutions.append(solution)
        coenergy_values = np.asarray(coenergy)
        spacing = 2.0 * np.pi / args.motor_angle_samples
        torque_values = -(
            np.roll(coenergy_values, -1) - np.roll(coenergy_values, 1)
        ) / (2.0 * spacing)
        for angle, energy, torque, solution in zip(
            angles, coenergy_values, torque_values, angle_solutions
        ):
            motor_angle_rows.append(
                {
                    "rotor_angle_rad": float(angle),
                    "coenergy_proxy_J": float(energy),
                    "torque_proxy_Nm": float(torque),
                    "average_joule_loss_W": float(
                        np.sum(solution.average_joule_loss)
                    ),
                    "residual_relative_norm": float(
                        solution.residual_relative_norm
                    ),
                    "solver_backend": solution.solver_backend,
                }
            )

    response_info = mixed.response_basis.diagnostics()
    bulk_current_gram = mixed.eddy_bases[0].mass_matrix()
    bulk_current_gram_identity_error = float(
        np.linalg.norm(bulk_current_gram - np.eye(bulk_current_gram.shape[0]))
        / max(float(np.linalg.norm(bulk_current_gram)), 1.0e-300)
    )
    checks = {
        "parent_order_admissible": ledger.is_parent_order_admissible(args.order),
        "hdiv_parent_order_admissible": (
            args.hdiv_order >= args.multipole_degree - 1
        ),
        "all_rhs_scaled_residuals_below_1e-8": all(
            row["residual_relative_norm"] < 1.0e-8 for row in rows
        ),
        "all_backward_errors_below_1e-12": all(
            row["residual_backward_error"] < 1.0e-12 for row in rows
        ),
        "all_losses_nonnegative": all(
            min(row["average_joule_loss"]) >= -1.0e-12 for row in rows
        ),
        "all_eddy_flux_density_probes_finite": all(
            np.isfinite(row["eddy_flux_density_probe_norm_T"]) for row in rows
        ),
        "parent_t_reconstructed": all(
            row["parent_t_dofs"] == hcurl_fes.ndof for row in rows
        ),
        "parent_hdiv_reconstructed": all(
            row["parent_hdiv_dofs"] == hdiv_reduction.parent_ndof for row in rows
        ),
        "demag_operator_is_active": (
            hdiv_reduction.diagnostics()["demag_frobenius_norm"] > 0.0
        ),
        "reduced_magnetic_operator_is_positive": (
            hdiv_reduction.diagnostics()["min_operator_eigenvalue"] > 0.0
        ),
        "hdiv_response_snapshots_converged": (
            hdiv_reduction.basis_generation["max_snapshot_relative_residual"]
            < 20.0 * args.hdiv_solve_tol
        ),
        "hdiv_training_responses_preserved": (
            hdiv_reduction.basis_generation["max_response_relative_energy_error"]
            < 1.0e-8
        ),
        "production_bdm_family_locked": (
            args.hdiv_family != "bdm"
            or (
                hdiv_reduction.parent_family == "BDM"
                and hdiv_reduction.parent_order == args.hdiv_order
            )
        ),
        "physical_hdiv_responses_are_protected": (
            hdiv_reduction.basis_generation["protected_physical_modes"] == 2
            and hdiv_reduction.basis_generation["pod_truncation_curve"][1][
                "max_physical_response_relative_energy_error"
            ]
            < 1.0e-8
        ),
        "all_three_eddy_blocks_reconstructed": all(
            row["eddy_blocks"] == ["volume", "volume1", "surface"] for row in rows
        ),
        "current_gram_rank_is_consistent": (
            response_info["current_gram_rank"]
            == sum(
                value > args.current_gram_rtol
                for value in response_info["current_gram_relative_eigenvalues"]
            )
        ),
        "bulk_current_gram_is_identity": (
            bulk_current_gram_identity_error < 1.0e-10
        ),
        "adjacency_roles_drive_reduction": (
            mixed.diagnostics()["eddy_block_roles"]
            == {"volume": "bulk", "volume1": "bridge", "surface": "sibc"}
            and mixed.adjacency_class_block_partition()
            == (("volume1", "surface"), ("volume",))
        ),
        "full_coupled_mixed_galerkin_is_exact": all(
            row["mixed_galerkin"]["full_coupled_schur"]
            and row["mixed_galerkin"]["schur_relative_error"] < 1.0e-8
            and row["direct_port_response_relative_error"] < 1.0e-7
            and row["direct_joule_loss_relative_error"] < 1.0e-7
            for row in rows
        ),
        "hdiv_hacapk_charge_gram_is_active": (
            hdiv_reduction.diagnostics()["demag_hmatrix_active"]
            and "ChargeGramHMatrix"
            in hdiv_reduction.diagnostics()["demag_hmatrix_backend"]
        ),
    }
    if args.motor_angle_samples:
        torque_values = np.asarray(
            [row["torque_proxy_Nm"] for row in motor_angle_rows]
        )
        torque_scale = max(float(np.max(np.abs(torque_values))), 1.0e-30)
        checks.update(
            {
                "motor_geometry_is_annular": args.geometry == "annular-motor",
                "motor_angle_grid_complete": len(motor_angle_rows)
                == args.motor_angle_samples,
                "motor_operator_reused_all_angles": True,
                "motor_torque_sign_reverses": bool(
                    torque_values.min() < 0.0 < torque_values.max()
                ),
                "motor_torque_is_nontrivial": torque_scale > 1.0e-12,
                "motor_angle_residuals_below_1e-8": all(
                    row["residual_relative_norm"] < 1.0e-8
                    for row in motor_angle_rows
                ),
                "motor_angle_losses_nonnegative": all(
                    row["average_joule_loss_W"] >= -1.0e-12
                    for row in motor_angle_rows
                ),
            }
        )
    checks["passed"] = all(checks.values())
    elapsed = time.perf_counter() - started
    return {
        "schema": "radia.validation.hcurl_vim_hdiv_mmm_end_to_end.v8",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "validation_host": platform.node(),
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "ngsolve_version": getattr(ng, "__version__", "unknown"),
        },
        "note": (
            "Correctness validation; timing fields are publishable only when "
            "validation_host is a designated compute host."
        ),
        "configuration": {
            "geometry": args.geometry,
            "axial_length_m": args.axial_length,
            "curvature_safety": args.curvature_safety,
            "response_backend": args.response_backend,
            "mixed_solver": args.mixed_solver,
            "motor_angle_samples": args.motor_angle_samples,
            "rotor_amplitude": args.rotor_amplitude,
            "shared_model_identity_sha256": args.shared_model_identity_sha256,
            "parent_order": args.order,
            "krylov_steps": args.steps,
            "hcurl_rtol": args.hcurl_rtol,
            "current_gram_rtol": args.current_gram_rtol,
            "frequencies_Hz": args.frequencies,
            "maxh": args.maxh,
            "sigma": args.sigma,
            "mu_r": args.mu_r,
            "intorder": args.intorder,
            "hdiv_order": args.hdiv_order,
            "hdiv_family": args.hdiv_family,
            "hdiv_intorder": args.hdiv_intorder,
            "demag_eps": args.demag_eps,
            "demag_eta": args.demag_eta,
            "multipole_degree": args.multipole_degree,
            "hdiv_pod_rtol": args.hdiv_pod_rtol,
            "hdiv_max_modes": args.hdiv_max_modes,
            "hdiv_solve_tol": args.hdiv_solve_tol,
            "kernel_epsilon": args.kernel_epsilon,
            "coupling_kernel_epsilon": args.coupling_kernel_epsilon,
            "parent_order_ledger": ledger.diagnostics(parent_order=args.order),
        },
        "hcurl_parent_ndof": int(hcurl_fes.ndof),
        "hcurl_current_gram_removed_modes": int(
            response_info["pre_current_gram_rank"]
            - response_info["current_gram_rank"]
        ),
        "bulk_current_gram_identity_error": bulk_current_gram_identity_error,
        "hdiv_reduction": hdiv_reduction.diagnostics(),
        "mixed_system": mixed.diagnostics(),
        "eddy_bubbling": mixed.eddy_bubbling.diagnostics(),
        "frequency_rows": rows,
        "motor_angle_rows": motor_angle_rows,
        "checks": checks,
        "elapsed_seconds": float(elapsed),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=6)
    parser.add_argument("--steps", type=int, default=8)
    parser.add_argument("--hcurl-rtol", type=float, default=1.0e-10)
    parser.add_argument("--current-gram-rtol", type=float, default=1.0e-10)
    parser.add_argument("--frequencies", type=_parse_floats, default=[100.0, 1.0e4])
    parser.add_argument("--maxh", type=float, default=2.5)
    parser.add_argument(
        "--geometry", choices=("notched-box", "annular-motor"),
        default="notched-box"
    )
    parser.add_argument("--axial-length", type=float, default=0.2)
    parser.add_argument("--curvature-safety", type=float, default=2.0)
    parser.add_argument(
        "--response-backend",
        choices=("auto", "operator", "static-condensed", "dense"),
        default="operator",
    )
    parser.add_argument(
        "--mixed-solver",
        choices=("gmres", "cocr", "dense"),
        default="gmres",
    )
    parser.add_argument("--motor-angle-samples", type=int, default=0)
    parser.add_argument("--rotor-amplitude", type=float, default=0.8)
    parser.add_argument("--shared-model-identity-sha256", default="")
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--mu-r", type=float, default=1001.0)
    parser.add_argument("--intorder", type=int, default=2)
    parser.add_argument("--hdiv-order", type=int, default=1)
    parser.add_argument("--hdiv-family", choices=("bdm", "rt"), default="bdm")
    parser.add_argument("--hdiv-intorder", type=int, default=2)
    parser.add_argument("--demag-eps", type=float, default=1.0e-7)
    parser.add_argument("--demag-eta", type=float, default=2.0)
    parser.add_argument("--multipole-degree", type=int, default=2)
    parser.add_argument("--hdiv-pod-rtol", type=float, default=1.0e-10)
    parser.add_argument("--hdiv-max-modes", type=int)
    parser.add_argument("--hdiv-solve-tol", type=float, default=1.0e-10)
    parser.add_argument("--kernel-epsilon", type=float, default=0.12)
    parser.add_argument("--coupling-kernel-epsilon", type=float, default=0.12)
    parser.add_argument("--bulk-degree", type=int, default=4)
    parser.add_argument("--bridge-trace-degree", type=int, default=0)
    parser.add_argument("--surface-current-degree", type=int, default=2)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    if args.order < 1 or args.steps < 1:
        parser.error("--order and --steps must be positive")
    if args.mu_r <= 1.0 or args.demag_eps <= 0.0 or args.demag_eta <= 0.0:
        parser.error("--mu-r must exceed 1 and demag parameters must be positive")
    if not 1 <= args.multipole_degree <= 3:
        parser.error("--multipole-degree must lie in [1, 3]")
    if args.hdiv_order not in (1, 2) or args.hdiv_intorder < 1:
        parser.error("--hdiv-order must be 1 or 2 and --hdiv-intorder must be positive")
    if args.hdiv_order < args.multipole_degree - 1:
        parser.error(
            "--hdiv-order must be at least --multipole-degree - 1 "
            "for polynomially admissible harmonic training"
        )
    if args.hdiv_pod_rtol <= 0.0 or args.hdiv_solve_tol <= 0.0:
        parser.error("HDiv POD and solve tolerances must be positive")
    if args.hdiv_max_modes is not None and args.hdiv_max_modes < 1:
        parser.error("--hdiv-max-modes must be positive")
    if args.axial_length <= 0.0:
        parser.error("--axial-length must be positive")
    if args.curvature_safety <= 0.0:
        parser.error("--curvature-safety must be positive")
    if args.motor_angle_samples and args.motor_angle_samples < 4:
        parser.error("--motor-angle-samples must be zero or at least four")
    if args.motor_angle_samples and args.geometry != "annular-motor":
        parser.error("motor angle sweeps require --geometry annular-motor")
    if args.rotor_amplitude <= 0.0:
        parser.error("--rotor-amplitude must be positive")
    if args.shared_model_identity_sha256 and (
        len(args.shared_model_identity_sha256) != 64
        or any(ch not in "0123456789abcdef" for ch in args.shared_model_identity_sha256)
    ):
        parser.error("--shared-model-identity-sha256 must be lowercase SHA-256")

    result = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("HCurl-VIM / HDiv-MMM end-to-end smoke")
    print(f"  output: {args.output}")
    print(
        f"  parent p={args.order}: {result['hcurl_parent_ndof']} DoF -> "
        f"{result['mixed_system']['hcurl_vim_modes']} hybrid eddy modes"
    )
    hdiv = result["hdiv_reduction"]
    generation = hdiv["basis_generation"]
    print(
        f"  HDiv-{args.hdiv_family.upper()}{args.hdiv_order}: "
        f"{hdiv['parent_ndof']} DoF -> {hdiv['reduced_modes']} modes "
        f"({generation['protected_physical_modes']} protected + "
        f"{generation['training_response_modes']} response-POD)"
    )
    for row in result["frequency_rows"]:
        print(
            f"  {row['frequency_Hz']:>9.1f} Hz  "
            f"res={row['residual_relative_norm']:.3e}  "
            f"loss={max(row['average_joule_loss']):.3e}  "
            f"solver={row['solver_backend']}"
        )
    print(f"  passed: {result['checks']['passed']}")
    return 0 if result["checks"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
