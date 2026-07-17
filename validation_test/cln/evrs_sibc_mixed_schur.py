"""EVRS bulk plus bridge-cycle plus surface-Omega/SIBC mixed Schur smoke.

This validation lane connects the IGTE mixed-Galerkin idea to the Radia VIM
API:

    HCurl(p) parent -> EVRS bulk curl(T) modes
    conductor-conductor graph cycles -> bridge current modes
    surface-Omega modes -> Z_s(s) M_Gamma
    mixed VIM matrix -> EVRS/bridge/surface Schur complements.

The run is a desktop smoke, not a benchmark.  It checks that eliminating the
bulk EVRS block while keeping bridge-cycle and surface blocks reproduces the
direct mixed-system solve and records how the port admittance changes with
Krylov depth.

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


def _frequency_result(
    system: vim.HybridVIMSystem,
    rhs: np.ndarray,
    *,
    frequency: float,
    sigma: float,
    surface_measure: float,
    dtn_pole_hz: float,
) -> dict[str, object]:
    s = 1j * 2.0 * np.pi * frequency
    zs = vim.SkinImpedance(s, sigma)
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

    return {
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
    frequencies: list[float],
    dtn_pole_hz: float,
    intorder: int,
    kernel_epsilon: float,
    rtol: float,
    parent_order_ledger: vim.EddyParentOrderLedger,
) -> dict[str, object]:
    t0 = time.perf_counter()
    response = _response_basis(
        fes,
        stiffness,
        mass,
        ports,
        steps=steps,
        condense=condense,
        rtol=rtol,
    )
    response_seconds = time.perf_counter() - t0

    names = [f"EVRS_p{order}_n{steps}_{i}" for i in range(response.rank)]
    built = vim.NgsolveTopologyAwareHybridVIM(
        mesh,
        fes,
        response.vectors,
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
            frequency=frequency,
            sigma=sigma,
            surface_measure=surface_measure,
            dtn_pole_hz=dtn_pole_hz,
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
        "bridge_cycle_modes": int(bridge_cycle.n_modes),
        "surface_modes": int(surface.n_modes),
        "total_reduced_modes": int(system.n_modes),
        "inactive_dofs": int(info["inactive_dofs"]),
        "eddy_visible_dofs": int(info["eddy_visible_dofs"]),
        "eddy_invisible_dofs": int(info["eddy_invisible_dofs"]),
        "compression_ratio": float(info["compression_ratio"]),
        "surface_measure": surface_measure,
        "response_seconds_lab_smoke": float(response_seconds),
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
                    frequencies=args.frequencies,
                    dtn_pole_hz=args.dtn_pole_hz,
                    intorder=args.intorder,
                    kernel_epsilon=args.kernel_epsilon,
                    rtol=args.rtol,
                    parent_order_ledger=parent_order_ledger,
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
    for case in cases:
        for row in case["frequency_rows"]:
            row["relative_error_to_reference"] = _relative_frobenius_error(
                row["_admittance_array"],
                reference_by_frequency[row["frequency_Hz"]],
            )
            del row["_admittance_array"]

    return {
        "schema": "radia.validation.evrs_sibc_mixed_schur.v1",
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
            "intorder": args.intorder,
            "kernel_epsilon": args.kernel_epsilon,
            "condense_from": args.condense_from,
            "rtol": args.rtol,
            "dtn_pole_hz": args.dtn_pole_hz,
            "parent_order_ledger": parent_order_ledger.diagnostics(),
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
    parser.add_argument("--intorder", type=int, default=2)
    parser.add_argument("--kernel-epsilon", type=float, default=0.12)
    parser.add_argument("--condense-from", type=int, default=7)
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    parser.add_argument("--dtn-pole-hz", type=float, default=0.0)
    parser.add_argument("--bulk-degree", type=int, default=4)
    parser.add_argument("--bridge-trace-degree", type=int, default=0)
    parser.add_argument("--surface-current-degree", type=int, default=2)
    parser.add_argument("--face-family", choices=("simplex", "tensor"), default="simplex")
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
