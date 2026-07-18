"""EVRS bulk plus bridge-cycle plus surface-Omega/SIBC mixed Schur validation.

This validation lane connects the IGTE mixed-Galerkin idea to the Radia VIM
API:

    HCurl(p) parent -> EVRS bulk curl(T) modes
    conductor-conductor graph cycles -> bridge current modes
    surface-Omega modes -> Z_s(s) M_Gamma
    mixed VIM matrix -> EVRS/bridge/surface Schur complements.

The run checks that eliminating the bulk EVRS block while keeping bridge-cycle
and surface blocks reproduces the direct mixed-system solve and records how
the port admittance changes with Krylov depth.  With ``--esim-bh-file`` it also
compares linear SIBC, uniform nonlinear ESIM, and solution-shaped local ESIM
on the same mixed-Galerkin operator.  Publish timings only from a designated
compute host.

Run from the repository root:

    python validation_test/cln/evrs_sibc_mixed_schur.py
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
DEFAULT_OUTPUT = HERE / "evrs_sibc_mixed_schur_smoke.json"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import radia.vim as vim

from evrs_pn_convergence import (  # noqa: E402
    _complex_matrix_parts,
    _parse_ints,
    _port_samples,
    _relative_frobenius_error,
    _sampled_basis_diagnostics,
)


def _parse_floats(text: str) -> list[float]:
    values = [float(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected at least one number")
    if any(value <= 0.0 for value in values):
        raise argparse.ArgumentTypeError("values must be positive")
    return values


def _response_basis(fes, stiffness, mass, ports, *, steps, condense, rtol):
    if condense:
        return vim.NgsolveStaticCondensedBlockKrylovBasis(
            stiffness,
            mass,
            ports,
            steps=steps,
            free_dofs=fes.FreeDofs(True),
            rtol=rtol,
        )
    return vim.NgsolveOperatorBlockKrylovBasis(
        stiffness,
        mass,
        ports,
        steps=steps,
        free_dofs=fes.FreeDofs(False),
        rtol=rtol,
    )


def _complex_scalar_parts(value) -> dict[str, float]:
    z = complex(value)
    return {"real": float(z.real), "imag": float(z.imag)}


def _complex_matrix_from_parts(parts: dict[str, object]) -> np.ndarray:
    return np.asarray(parts["real"], dtype=float) + 1j * np.asarray(
        parts["imag"],
        dtype=float,
    )


def _make_skin_shape(geometry: str):
    import netgen.occ as occ

    if geometry == "box":
        shape = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    elif geometry == "bar":
        shape = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(2, 0.5, 0.5))
    elif geometry == "notched-box":
        outer = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
        notch = occ.Box(occ.Pnt(0.45, 0.45, -0.1), occ.Pnt(1.1, 1.1, 1.1))
        shape = outer - notch
    elif geometry == "l-prism":
        leg_x = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 0.45, 1))
        leg_y = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(0.45, 1, 1))
        shape = leg_x + leg_y
    elif geometry == "sphere":
        shape = occ.Sphere(occ.Pnt(0, 0, 0), 1.0)
    else:
        raise ValueError(f"unknown geometry: {geometry}")
    shape.mat("cond")
    for face in shape.faces:
        face.name = "skin"
    return shape


def _make_skin_mesh(
    maxh: float,
    geometry: str,
    curve_order: int = 1,
    corner_edge_maxh: float | None = None,
):
    import ngsolve as ng
    import netgen.occ as occ

    shape = _make_skin_shape(geometry)
    if corner_edge_maxh is not None:
        if geometry not in {"notched-box", "l-prism"}:
            raise ValueError("corner_edge_maxh requires notched-box or l-prism")
        for edge in shape.edges:
            center = np.asarray(tuple(edge.center), dtype=float)
            if np.linalg.norm(center[:2] - (0.45, 0.45)) < 1.0e-8:
                edge.maxh = float(corner_edge_maxh)
    mesh = ng.Mesh(occ.OCCGeometry(shape).GenerateMesh(maxh=maxh))
    if curve_order > 1:
        mesh.Curve(curve_order)
    return mesh


def _assemble_parent(
    order: int,
    maxh: float,
    condense: bool,
    geometry: str,
    curve_order: int = 1,
):
    import ngsolve as ng

    mesh = _make_skin_mesh(maxh, geometry, curve_order)
    fes = ng.HCurl(mesh, order=order, nograds=True)
    u, v = fes.TnT()

    stiffness = ng.BilinearForm(fes, condense=condense)
    stiffness += ng.curl(u) * ng.curl(v) * ng.dx + 0.05 * u * v * ng.dx
    mass = ng.BilinearForm(fes)
    mass += u * v * ng.dx

    ports = []
    for cf in (
        ng.CoefficientFunction((-ng.y, ng.x, 0.0)),
        ng.CoefficientFunction((0.0, -ng.z, ng.y)),
    ):
        port = ng.LinearForm(fes)
        port += cf * v * ng.dx
        ports.append(port)

    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        for port in ports:
            port.Assemble()

    return mesh, fes, stiffness, mass, ports


def _surface_grad_modes():
    import ngsolve as ng

    return (
        ng.CoefficientFunction((1.0, 0.0, 0.0)),
        ng.CoefficientFunction((0.0, 1.0, 0.0)),
        ng.CoefficientFunction((0.0, 0.0, 1.0)),
    )


def _port0(points: np.ndarray) -> np.ndarray:
    return _port_samples(points)[0]


def _port1(points: np.ndarray) -> np.ndarray:
    return _port_samples(points)[1]


def _basis_port_rhs(basis) -> np.ndarray:
    return np.column_stack(
        [
            vim.ExternalVectorPotentialRHS(basis, samples)
            for samples in _port_samples(basis.points)
        ]
    )


def _block_rhs(system: vim.HybridVIMSystem, volume, bridge, surface) -> np.ndarray:
    volume_rhs = np.column_stack(
        [
            vim.ExternalVectorPotentialRHS(volume, samples)
            for samples in _port_samples(volume.points)
        ]
    )
    return system.block_rhs(
        volume=volume_rhs,
        volume1=_basis_port_rhs(bridge),
        surface=_basis_port_rhs(surface),
    )


def _schur_solution_checks(
    system: vim.HybridVIMSystem,
    s,
    rhs: np.ndarray,
    *,
    keep_blocks,
    eliminate_blocks,
    surface_impedance,
) -> dict[str, float]:
    z = system.impedance(s, surface_impedance=surface_impedance)
    keep = system.block_indices(keep_blocks)
    eliminate = system.block_indices(eliminate_blocks)

    z_ee = z[np.ix_(eliminate, eliminate)]
    z_ek = z[np.ix_(eliminate, keep)]
    z_ke = z[np.ix_(keep, eliminate)]
    b_e = rhs[eliminate, :]
    b_k = rhs[keep, :]

    s_keep = system.schur_complement_blocks(
        keep_blocks,
        eliminate_blocks,
        s,
        surface_impedance=surface_impedance,
    )
    b_keep = b_k - z_ke @ np.linalg.solve(z_ee, b_e)
    x_keep = np.linalg.solve(s_keep, b_keep)
    x_eliminate = np.linalg.solve(z_ee, b_e - z_ek @ x_keep)

    x_schur = np.zeros_like(rhs, dtype=np.result_type(rhs, z))
    x_schur[keep, :] = x_keep
    x_schur[eliminate, :] = x_eliminate
    x_direct = np.linalg.solve(z, rhs)
    return {
        "solution_relative_error_to_direct": _relative_frobenius_error(
            x_schur,
            x_direct,
        ),
        "residual_relative_error": _relative_frobenius_error(
            z @ x_schur,
            rhs,
        ),
        "direct_residual_relative_error": _relative_frobenius_error(
            z @ x_direct,
            rhs,
        ),
    }


def _mixed_port_metrics(
    system: vim.HybridVIMSystem,
    rhs: np.ndarray,
    s,
    surface_impedance,
) -> dict[str, object]:
    """Compare the mixed-Galerkin solve with the unreduced hybrid solve."""

    z = system.impedance(s, surface_impedance=surface_impedance)
    t0 = time.perf_counter()
    direct = np.linalg.solve(z, rhs)
    direct_seconds = time.perf_counter() - t0
    direct_admittance = rhs.conj().T @ direct

    reduction = system.mixed_galerkin_orthogonalization(
        ("volume1", "surface"),
        "volume",
        s,
        surface_impedance=surface_impedance,
    )
    t0 = time.perf_counter()
    mixed = reduction.solve(z, rhs)
    mixed_seconds = time.perf_counter() - t0
    mixed_admittance = rhs.conj().T @ mixed
    rhs_scale = max(float(np.linalg.norm(rhs)), np.finfo(float).tiny)

    if isinstance(surface_impedance, vim.SurfaceImpedanceGram):
        surface_loss = surface_impedance.dissipative_matrix
    else:
        surface_loss = complex(surface_impedance).real * system.surface_mass
    dissipative = system.resistance + surface_loss
    average_loss = 0.5 * np.real(
        np.einsum("ip,ij,jp->p", mixed.conj(), dissipative, mixed)
    )
    return {
        "full_modes": int(system.n_modes),
        "retained_modes": int(reduction.rank),
        "direct_solve_seconds": float(direct_seconds),
        "mixed_solve_seconds": float(mixed_seconds),
        "port_admittance": _complex_matrix_parts(mixed_admittance),
        "port_relative_error_to_direct": _relative_frobenius_error(
            mixed_admittance,
            direct_admittance,
        ),
        "solution_relative_error_to_direct": _relative_frobenius_error(
            mixed,
            direct,
        ),
        "residual_relative_norm": float(np.linalg.norm(z @ mixed - rhs) / rhs_scale),
        "average_loss_by_port": [float(value) for value in average_loss],
        "orthogonalization": reduction.diagnostics(),
        "_solution": mixed,
        "_admittance_array": mixed_admittance,
    }


def _normalized_surface_field_profile(
    system: vim.HybridVIMSystem,
    surface_basis,
    solution: np.ndarray,
    target_rms: float,
) -> np.ndarray:
    """Recover a solution-shaped ``|H_t|`` profile at SIBC samples."""

    coefficients = solution[system.block_slice("surface"), :]
    currents = np.einsum("mp,mik->pik", coefficients, surface_basis.modes)
    magnitude = np.sqrt(np.mean(np.sum(np.abs(currents) ** 2, axis=2), axis=0))
    rms = np.sqrt(
        np.sum(surface_basis.weights * magnitude**2)
        / np.sum(surface_basis.weights)
    )
    if not np.isfinite(rms) or rms <= np.finfo(float).tiny:
        return np.full(surface_basis.n_samples, target_rms, dtype=float)
    profile = magnitude * (target_rms / rms)
    return np.maximum(profile, target_rms * 1.0e-3)


def _single_port_surface_field_profile(
    system: vim.HybridVIMSystem,
    surface_basis,
    solution: np.ndarray,
) -> np.ndarray:
    """Recover physical ``|H_t|`` samples from one reduced excitation."""

    coefficients = np.asarray(solution)[system.block_slice("surface")]
    if coefficients.shape != (surface_basis.n_modes,):
        raise ValueError("solution must contain exactly one physical excitation")
    currents = np.einsum("m,mik->ik", coefficients, surface_basis.modes)
    return np.sqrt(np.sum(np.abs(currents) ** 2, axis=1))


def _frequency_result(
    system: vim.HybridVIMSystem,
    rhs: np.ndarray,
    surface_basis,
    *,
    frequency: float,
    sigma: float,
    sibc_mu_r: float,
    surface_measure: float,
    dtn_pole_hz: float,
    esim_configuration: dict[str, object] | None,
) -> dict[str, object]:
    s = 1j * 2.0 * np.pi * frequency
    zs = vim.SkinImpedance(s, sigma, vim.MU0 * sibc_mu_r)
    z = system.impedance(s, surface_impedance=zs)
    y = system.port_admittance(s, rhs, surface_impedance=zs)
    s_bridge_surface = system.schur_complement_blocks(
        ("volume1", "surface"),
        "volume",
        s,
        surface_impedance=zs,
    )
    s_surface = system.schur_complement_blocks(
        "surface",
        ("volume", "volume1"),
        s,
        surface_impedance=zs,
    )
    s_evrs_bridge = system.schur_complement_blocks(
        ("volume", "volume1"),
        "surface",
        s,
        surface_impedance=zs,
    )
    k_sibc = surface_measure * np.sqrt(sigma / vim.MU0)
    dtn_pole = 2.0 * np.pi * dtn_pole_hz
    schur_checks = _schur_solution_checks(
        system,
        s,
        rhs,
        keep_blocks=("volume1", "surface"),
        eliminate_blocks="volume",
        surface_impedance=zs,
    )
    mixed_orthogonalization = system.mixed_galerkin_orthogonalization(
        ("volume1", "surface"),
        "volume",
        s,
        surface_impedance=zs,
    )

    result = {
        "frequency_Hz": float(frequency),
        "s_rad_per_s": _complex_scalar_parts(s),
        "skin_impedance": _complex_scalar_parts(zs),
        "mixed_impedance_condition": float(np.linalg.cond(z)),
        "bridge_surface_schur_norm": float(np.linalg.norm(s_bridge_surface)),
        "surface_schur_norm": float(np.linalg.norm(s_surface)),
        "surface_schur_inverse_norm": float(np.linalg.norm(np.linalg.inv(s_surface))),
        "evrs_bridge_schur_norm": float(np.linalg.norm(s_evrs_bridge)),
        "schur_solution_relative_error_to_direct": schur_checks[
            "solution_relative_error_to_direct"
        ],
        "schur_residual_relative_error": schur_checks["residual_relative_error"],
        "mixed_galerkin_orthogonalization": mixed_orthogonalization.diagnostics(),
        "direct_residual_relative_error": schur_checks[
            "direct_residual_relative_error"
        ],
        "sibc_tail_admittance": _complex_scalar_parts(
            vim.SIBCAdmittanceTail(s, surface_measure, sigma)
        ),
        "dtn_sibc_scalar_admittance": _complex_scalar_parts(
            vim.SIBCSchurTerminationAdmittance(s, k_sibc, d=dtn_pole)
        ),
        "port_admittance": _complex_matrix_parts(y),
        "_admittance_array": y,
    }
    linear_metrics = _mixed_port_metrics(system, rhs, s, zs)
    result["ablation"] = {
        "linear_sibc_mixed_galerkin": {
            key: value
            for key, value in linear_metrics.items()
            if not key.startswith("_")
        },
    }

    if esim_configuration is not None:
        target_h_rms = float(esim_configuration["h_rms"])
        production_model = esim_configuration["_model"]
        h_profile = _normalized_surface_field_profile(
            system,
            surface_basis,
            linear_metrics["_solution"],
            target_h_rms,
        )
        z_samples, cell_diagnostics = production_model.impedance_samples(
            h_profile,
            frequency,
        )
        z_uniform, _ = production_model.impedance_samples(
            np.asarray([float(esim_configuration["h_rms"])]),
            frequency,
        )
        z_uniform = z_uniform[0]
        uniform_metrics = _mixed_port_metrics(system, rhs, s, z_uniform)

        initial_port_profile = _single_port_surface_field_profile(
            system,
            surface_basis,
            linear_metrics["_solution"][:, 0],
        )
        initial_port_rms = float(
            np.sqrt(
                np.sum(surface_basis.weights * initial_port_profile**2)
                / np.sum(surface_basis.weights)
            )
        )
        if not np.isfinite(initial_port_rms) or initial_port_rms <= 0.0:
            raise RuntimeError("production ESIM port has zero tangential field")
        production_solution = vim.SolveLocalESIMSurfaceVIM(
            system,
            surface_basis,
            rhs[:, 0],
            production_model,
            frequency,
            outer_tolerance=float(esim_configuration["outer_tol"]),
            outer_max_iterations=int(esim_configuration["outer_max_iterations"]),
            outer_relaxation=float(esim_configuration["outer_relaxation"]),
            field_amplitude_scale=target_h_rms / initial_port_rms,
            mixed_galerkin_keep_blocks=("volume1", "surface"),
            mixed_galerkin_eliminate_blocks="volume",
        )
        production_operator = system.impedance(
            s,
            surface_impedance=production_solution.surface_impedance,
        )
        production_residual = float(
            np.linalg.norm(
                production_operator @ production_solution.coefficients - rhs[:, 0]
            )
            / max(float(np.linalg.norm(rhs[:, 0])), np.finfo(float).tiny)
        )

        outer_converged = False
        outer_relative_change = float("inf")
        local_metrics = None
        local_gram = None
        for outer_iteration in range(int(esim_configuration["outer_max_iterations"])):
            local_gram = vim.AssembleSurfaceImpedanceGram(
                system,
                surface_basis,
                z_samples,
                label="esim-local-sibc",
            )
            local_metrics = _mixed_port_metrics(system, rhs, s, local_gram)
            updated_h_profile = _normalized_surface_field_profile(
                system,
                surface_basis,
                local_metrics["_solution"],
                float(esim_configuration["h_rms"]),
            )
            updated_z_samples, _ = production_model.impedance_samples(
                updated_h_profile,
                frequency,
            )
            outer_relative_change = float(
                np.linalg.norm(updated_z_samples - z_samples)
                / max(float(np.linalg.norm(z_samples)), np.finfo(float).tiny)
            )
            h_profile = updated_h_profile
            if outer_relative_change < float(esim_configuration["outer_tol"]):
                outer_converged = True
                break
            relaxation = float(esim_configuration["outer_relaxation"])
            z_samples = (
                (1.0 - relaxation) * z_samples
                + relaxation * updated_z_samples
            )
        if local_metrics is None or local_gram is None:
            raise RuntimeError("ESIM outer iteration did not execute")
        linear_admittance = linear_metrics["_admittance_array"]
        uniform_admittance = uniform_metrics["_admittance_array"]
        local_admittance = local_metrics["_admittance_array"]
        result["esim"] = {
            **cell_diagnostics,
            "target_h_rms_A_per_m": float(esim_configuration["h_rms"]),
            "outer_converged": outer_converged,
            "outer_iterations": int(outer_iteration + 1),
            "outer_relative_impedance_change": outer_relative_change,
            "recovered_h_rms_A_per_m": float(
                np.sqrt(
                    np.sum(surface_basis.weights * h_profile**2)
                    / np.sum(surface_basis.weights)
                )
            ),
            "h_min_A_per_m": float(np.min(h_profile)),
            "h_max_A_per_m": float(np.max(h_profile)),
            "uniform_surface_impedance_ohm": _complex_scalar_parts(z_uniform),
            "local_surface_impedance": local_gram.diagnostics(),
            "uniform_to_local_port_relative_difference": _relative_frobenius_error(
                uniform_admittance,
                local_admittance,
            ),
            "linear_to_local_port_relative_difference": _relative_frobenius_error(
                linear_admittance,
                local_admittance,
            ),
            "production_api": {
                **production_solution.diagnostics(),
                "rhs_port": 0,
                "fixed_field_amplitude_scale": float(target_h_rms / initial_port_rms),
                "residual_relative_norm": production_residual,
                "port_admittance": _complex_scalar_parts(
                    np.vdot(rhs[:, 0], production_solution.coefficients)
                ),
            },
        }
        result["ablation"]["uniform_esim_sibc_mixed_galerkin"] = {
            key: value
            for key, value in uniform_metrics.items()
            if not key.startswith("_")
        }
        result["ablation"]["local_esim_sibc_mixed_galerkin"] = {
            key: value
            for key, value in local_metrics.items()
            if not key.startswith("_")
        }
    return result


def _case_result(
    mesh,
    fes,
    stiffness,
    mass,
    ports,
    *,
    geometry: str,
    order: int,
    steps: int,
    condense: bool,
    sigma: float,
    sibc_mu_r: float,
    frequencies: list[float],
    dtn_pole_hz: float,
    intorder: int,
    kernel_epsilon: float,
    rtol: float,
    parent_order_ledger: vim.EddyParentOrderLedger,
    esim_configuration: dict[str, object] | None,
) -> dict[str, object]:
    t0 = time.perf_counter()
    parent_response = _response_basis(
        fes,
        stiffness,
        mass,
        ports,
        steps=steps,
        condense=condense,
        rtol=rtol,
    )
    response_seconds = time.perf_counter() - t0

    names = [
        f"EVRS_p{order}_n{steps}_{i}"
        for i in range(parent_response.rank)
    ]
    t0 = time.perf_counter()
    built = vim.NgsolveTopologyAwareHybridVIM(
        mesh,
        fes,
        parent_response,
        _surface_grad_modes(),
        sigma=sigma,
        conductive_materials="cond",
        surface_boundaries="skin",
        intorder=intorder,
        kernel_epsilon=kernel_epsilon,
        volume_names=names,
        surface_names=("omega_x", "omega_y", "omega_z"),
        parent_order_ledger=parent_order_ledger,
        port_vector_potentials=(_port0, _port1),
    )
    vim_assembly_seconds = time.perf_counter() - t0
    volume = built.volume_basis
    bridge_cycle = built.bridge_cycle_basis
    surface = built.surface_basis
    topology = built.topology
    dof_policy = built.dof_policy
    conductor_graph = built.conductor_graph
    system = built.system
    rhs = built.rhs
    assert rhs is not None
    surface_measure = float(surface.weights.sum())
    response = built.response_basis
    if response is None:
        raise RuntimeError("topology-aware VIM did not retain the response basis")
    info = response.diagnostics()
    bridge_cycle_system = vim.AssembleHybridVIM(
        bridge_cycle,
        sigma=sigma,
        kernel_epsilon=kernel_epsilon,
    )
    reduction_plan = dof_policy.reduction_plan(
        evrs_rank=response.rank,
        surface_modes=surface.n_modes,
        loop_bridge_modes=conductor_graph.cycle_rank,
        bridge_strategy="cycle-basis",
    )
    frequency_rows = [
        _frequency_result(
            system,
            rhs,
            surface,
            frequency=frequency,
            sigma=sigma,
            sibc_mu_r=sibc_mu_r,
            surface_measure=surface_measure,
            dtn_pole_hz=dtn_pole_hz,
            esim_configuration=esim_configuration,
        )
        for frequency in frequencies
    ]

    return {
        "order": order,
        "geometry": geometry,
        "krylov_steps": steps,
        "condensed": condense,
        "ndof": int(info["ndof"]),
        "active_dofs": int(info["active_dofs"]),
        "rank": int(info["rank"]),
        "pre_current_gram_rank": int(
            info.get("pre_current_gram_rank", parent_response.rank)
        ),
        "bridge_cycle_modes": int(bridge_cycle.n_modes),
        "surface_modes": int(surface.n_modes),
        "total_reduced_modes": int(system.n_modes),
        "inactive_dofs": int(info["inactive_dofs"]),
        "eddy_visible_dofs": int(info["eddy_visible_dofs"]),
        "eddy_invisible_dofs": int(info["eddy_invisible_dofs"]),
        "compression_ratio": float(info["compression_ratio"]),
        "surface_measure": surface_measure,
        "response_basis_seconds": float(response_seconds),
        "vim_assembly_seconds": float(vim_assembly_seconds),
        "reduced_dense_storage_bytes": int(
            system.resistance.nbytes
            + system.surface_mass.nbytes
            + (
                0
                if hasattr(system.inductance, "matvec")
                else np.asarray(system.inductance).nbytes
            )
        ),
        "volume_samples": int(volume.n_samples),
        "surface_samples": int(surface.n_samples),
        "blocks": dict(system.blocks),
        "topology_diagnostics": topology.diagnostics(),
        "conductor_graph_diagnostics": conductor_graph.diagnostics(),
        "bridge_cycle_basis_diagnostics": _sampled_basis_diagnostics(bridge_cycle),
        "bridge_cycle_vim_diagnostics": bridge_cycle_system.diagnostics(),
        "dof_policy_diagnostics": dof_policy.diagnostics(),
        "reduction_plan_diagnostics": reduction_plan.diagnostics(),
        "topology_aware_hybrid_vim_diagnostics": built.diagnostics(),
        "system_diagnostics": system.diagnostics(),
        "frequency_rows": frequency_rows,
    }


def run_sweep(args: argparse.Namespace) -> dict[str, object]:
    import ngsolve as ng

    parent_order_ledger = vim.EddyParentOrderLedger(
        bulk_degree=args.bulk_degree,
        bridge_trace_degree=args.bridge_trace_degree,
        surface_current_degree=args.surface_current_degree,
        face_family=args.face_family,
    )
    esim_configuration = None
    if args.esim_bh_file is not None:
        bh_curve = np.loadtxt(args.esim_bh_file)
        lut_h_min = (
            max(args.esim_h_rms * 1.0e-4, 1.0e-6)
            if args.esim_lut_h_min is None
            else float(args.esim_lut_h_min)
        )
        lut_h_max = (
            args.esim_h_rms * 10.0
            if args.esim_lut_h_max is None
            else float(args.esim_lut_h_max)
        )
        cell_model = vim.LocalESIMSurfaceModel(
            bh_curve=bh_curve,
            sigma=args.sigma,
            bins=args.esim_bins,
            n_nodes=args.esim_n_nodes,
            cell_tolerance=args.esim_tol,
            cell_max_iterations=args.esim_max_iterations,
            cell_relaxation=args.esim_relaxation,
            h_floor=lut_h_min,
        )
        lut_start = time.perf_counter()
        if args.esim_lut_input is None:
            lut = vim.BuildLocalESIMSurfaceLUT(
                cell_model,
                args.frequencies,
                np.geomspace(
                    lut_h_min,
                    lut_h_max,
                    args.esim_lut_field_nodes,
                ),
                output_path=args.esim_lut_output,
            )
            lut_source = "built-offline"
        else:
            lut = vim.LocalESIMSurfaceLUT.load(args.esim_lut_input)
            lut.validate_model(cell_model)
            if args.esim_lut_output is not None:
                lut.save(args.esim_lut_output)
            lut_source = "loaded-persistent-npz"
        lut_wall_time = time.perf_counter() - lut_start
        lut_validation = None
        lut_validation_wall_time = None
        if args.esim_lut_validate:
            validation_start = time.perf_counter()
            lut_validation = vim.ValidateLocalESIMSurfaceLUT(cell_model, lut)
            lut_validation_wall_time = time.perf_counter() - validation_start
        esim_configuration = {
            "bh_curve": bh_curve,
            "bh_file": str(args.esim_bh_file),
            "h_rms": args.esim_h_rms,
            "bins": args.esim_bins,
            "n_nodes": args.esim_n_nodes,
            "tol": args.esim_tol,
            "max_iterations": args.esim_max_iterations,
            "relaxation": args.esim_relaxation,
            "outer_tol": args.esim_outer_tol,
            "outer_max_iterations": args.esim_outer_max_iterations,
            "outer_relaxation": args.esim_outer_relaxation,
            "lut": {
                **lut.diagnostics(),
                "source": lut_source,
                "offline_build_or_load_wall_time_s": float(lut_wall_time),
                "validation": lut_validation,
                "validation_wall_time_s": lut_validation_wall_time,
            },
            "_model": cell_model.with_lut(lut),
        }
    orders = args.orders if args.orders is not None else [args.order]
    cases = []
    for order in orders:
        condense = order >= args.condense_from
        mesh, fes, stiffness, mass, ports = _assemble_parent(
            order,
            args.maxh,
            condense,
            args.geometry,
            args.curve_order,
        )
        for steps in args.steps:
            cases.append(
                _case_result(
                    mesh,
                    fes,
                    stiffness,
                    mass,
                    ports,
                    geometry=args.geometry,
                    order=order,
                    steps=steps,
                    condense=condense,
                    sigma=args.sigma,
                    sibc_mu_r=args.sibc_mu_r,
                    frequencies=args.frequencies,
                    dtn_pole_hz=args.dtn_pole_hz,
                    intorder=args.intorder,
                    kernel_epsilon=args.kernel_epsilon,
                    rtol=args.rtol,
                    parent_order_ledger=parent_order_ledger,
                    esim_configuration=esim_configuration,
                )
            )

    reference_order = max(orders)
    reference_steps = max(args.steps)
    reference_case = next(
        case
        for case in cases
        if case["order"] == reference_order and case["krylov_steps"] == reference_steps
    )
    reference_by_frequency = {
        row["frequency_Hz"]: row["_admittance_array"]
        for row in reference_case["frequency_rows"]
    }
    reference_ablation_by_frequency = {
        row["frequency_Hz"]: {
            name: _complex_matrix_from_parts(metrics["port_admittance"])
            for name, metrics in row.get("ablation", {}).items()
        }
        for row in reference_case["frequency_rows"]
    }
    for case in cases:
        for row in case["frequency_rows"]:
            row["relative_error_to_reference"] = _relative_frobenius_error(
                row["_admittance_array"],
                reference_by_frequency[row["frequency_Hz"]],
            )
            for name, metrics in row.get("ablation", {}).items():
                metrics["relative_error_to_reference"] = _relative_frobenius_error(
                    _complex_matrix_from_parts(metrics["port_admittance"]),
                    reference_ablation_by_frequency[row["frequency_Hz"]][name],
                )
            del row["_admittance_array"]

    return {
        "schema": "radia.validation.evrs_sibc_mixed_schur.v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "validation_host": platform.node(),
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "ngsolve_version": getattr(ng, "__version__", "unknown"),
        },
        "note": (
            "Timing fields are wall-clock observations on validation_host; "
            "publish them only when generated on a designated compute host."
        ),
        "configuration": {
            "order": args.order,
            "orders": orders,
            "geometry": args.geometry,
            "krylov_steps": args.steps,
            "reference": {"order": reference_order, "krylov_steps": reference_steps},
            "frequencies_Hz": args.frequencies,
            "maxh": args.maxh,
            "curve_order": args.curve_order,
            "sigma": args.sigma,
            "sibc_mu_r": args.sibc_mu_r,
            "intorder": args.intorder,
            "kernel_epsilon": args.kernel_epsilon,
            "condense_from": args.condense_from,
            "rtol": args.rtol,
            "dtn_pole_hz": args.dtn_pole_hz,
            "parent_order_ledger": parent_order_ledger.diagnostics(),
            "esim": (
                None
                if esim_configuration is None
                else {
                    key: value
                    for key, value in esim_configuration.items()
                    if not key.startswith("_") and key != "bh_curve"
                }
            ),
        },
        "cases": cases,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=6)
    parser.add_argument("--orders", type=_parse_ints, default=None)
    parser.add_argument(
        "--geometry",
        choices=("box", "bar", "notched-box", "l-prism", "sphere"),
        default="box",
    )
    parser.add_argument("--steps", type=_parse_ints, default=[8, 11, 12])
    parser.add_argument("--frequencies", type=_parse_floats, default=[100.0, 1.0e4, 1.0e6])
    parser.add_argument("--maxh", type=float, default=3.0)
    parser.add_argument("--curve-order", type=int, default=1)
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--sibc-mu-r", type=float, default=1.0)
    parser.add_argument("--intorder", type=int, default=2)
    parser.add_argument("--kernel-epsilon", type=float, default=0.12)
    parser.add_argument("--condense-from", type=int, default=7)
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    parser.add_argument("--dtn-pole-hz", type=float, default=0.0)
    parser.add_argument("--bulk-degree", type=int, default=4)
    parser.add_argument("--bridge-trace-degree", type=int, default=0)
    parser.add_argument("--surface-current-degree", type=int, default=2)
    parser.add_argument("--face-family", choices=("simplex", "tensor"), default="simplex")
    parser.add_argument("--esim-bh-file", type=Path, default=None)
    parser.add_argument("--esim-h-rms", type=float, default=1.0e3)
    parser.add_argument("--esim-bins", type=int, default=12)
    parser.add_argument("--esim-n-nodes", type=int, default=100)
    parser.add_argument("--esim-tol", type=float, default=1.0e-5)
    parser.add_argument("--esim-max-iterations", type=int, default=80)
    parser.add_argument("--esim-relaxation", type=float, default=0.5)
    parser.add_argument("--esim-outer-tol", type=float, default=1.0e-3)
    parser.add_argument("--esim-outer-max-iterations", type=int, default=20)
    parser.add_argument("--esim-outer-relaxation", type=float, default=0.5)
    parser.add_argument("--esim-lut-h-min", type=float, default=None)
    parser.add_argument("--esim-lut-h-max", type=float, default=None)
    parser.add_argument("--esim-lut-field-nodes", type=int, default=96)
    parser.add_argument(
        "--esim-lut-validate",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--esim-lut-input", type=Path, default=None)
    parser.add_argument("--esim-lut-output", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    if args.order < 1:
        parser.error("--order must be positive")
    if args.orders is not None and any(order < 1 for order in args.orders):
        parser.error("--orders must contain positive integers")
    if args.dtn_pole_hz < 0.0:
        parser.error("--dtn-pole-hz must be non-negative")
    if args.curve_order < 1:
        parser.error("--curve-order must be positive")
    if args.sibc_mu_r <= 0.0:
        parser.error("--sibc-mu-r must be positive")
    if args.esim_h_rms <= 0.0:
        parser.error("--esim-h-rms must be positive")
    if args.esim_bins < 2 or args.esim_n_nodes < 3:
        parser.error("--esim-bins must be >= 2 and --esim-n-nodes must be >= 3")
    if args.esim_lut_field_nodes < 2:
        parser.error("--esim-lut-field-nodes must be >= 2")
    if args.esim_tol <= 0.0 or args.esim_max_iterations < 1:
        parser.error("ESIM tolerance and iteration limit must be positive")
    if not 0.0 < args.esim_relaxation <= 1.0:
        parser.error("--esim-relaxation must be in (0, 1]")
    if args.esim_outer_tol <= 0.0 or args.esim_outer_max_iterations < 1:
        parser.error("ESIM outer tolerance and iteration limit must be positive")
    if not 0.0 < args.esim_outer_relaxation <= 1.0:
        parser.error("--esim-outer-relaxation must be in (0, 1]")
    if args.esim_bh_file is not None and not args.esim_bh_file.is_file():
        parser.error(f"--esim-bh-file does not exist: {args.esim_bh_file}")
    if args.esim_lut_input is not None and not args.esim_lut_input.is_file():
        parser.error(f"--esim-lut-input does not exist: {args.esim_lut_input}")
    lut_h_min = (
        max(args.esim_h_rms * 1.0e-4, 1.0e-6)
        if args.esim_lut_h_min is None
        else args.esim_lut_h_min
    )
    lut_h_max = (
        args.esim_h_rms * 10.0
        if args.esim_lut_h_max is None
        else args.esim_lut_h_max
    )
    if lut_h_min <= 0.0 or lut_h_max <= lut_h_min:
        parser.error("ESIM LUT field range must satisfy 0 < h_min < h_max")

    result = run_sweep(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    print("EVRS + bridge-cycle + surface-Omega/SIBC mixed Schur smoke")
    print(f"  output: {args.output}")
    print("  p  n  rank/active  total  freq_Hz    relerr      schur_res")
    for case in result["cases"]:
        for row in case["frequency_rows"]:
            print(
                f"  {case['order']:>1}  {case['krylov_steps']:>2}  "
                f"{case['rank']:>4}/{case['active_dofs']:<5}  "
                f"{case['total_reduced_modes']:>5}  "
                f"{row['frequency_Hz']:>8.1f}  "
                f"{row['relative_error_to_reference']:.3e}  "
                f"{row['schur_residual_relative_error']:.3e}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
