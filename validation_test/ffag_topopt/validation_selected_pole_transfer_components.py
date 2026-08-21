"""Manufactured HEX-pole inverse for selected transfer-map components.

One adjacent BDM1 HEX pole block is added to create a physically realizable
target.  Starting from the connected seed block, HDiv-MMM must recover that
binary pole structure while constraining only the named focusing and
dispersion entries.  The full target map is still required to be symplectic;
unselected entries are not treated as independent design variables.

This is numerical correctness validation, not a performance benchmark.
"""
from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import numpy as np


CONTROLLED_COMPONENTS = (
    "horizontal_focusing",
    "vertical_focusing",
    "horizontal_dispersion",
)


def _one_segment_arc(*, radius=8.0, angle=0.08, rigidity=1.5):
    from radia.accelerator_magnet_topopt import PlanarDesignOrbit

    positions = np.array([
        [0.0, 0.0, 2.0],
        [radius*np.sin(angle), radius*(1.0-np.cos(angle)), 2.0],
    ])
    tangents = np.array([
        [1.0, 0.0, 0.0],
        [np.cos(angle), np.sin(angle), 0.0],
    ])
    return PlanarDesignOrbit(
        positions, tangents, magnetic_rigidity=rigidity,
        bend_axis=np.array([0.0, 0.0, 1.0]))


def run(args):
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh

    from radia.accelerator_magnet_topopt import (
        MultiMomentumTransferMatrixObjective,
        PlanarDesignOrbit,
        build_multi_orbit_field_response_matrix,
        optimize_hdiv_mmm_magnet_from_transfer_matrices,
        solve_transfer_matrix_field_correction,
        static_magnet_symplectic_residual,
        static_magnet_transfer_component_entries,
    )
    from radia.ffag_topopt import (
        build_ffag_fixed_design_orbit_target_family,
    )
    from radia.isochronous_topopt import (
        MU0,
        combined_function_transfer_map_from_field_response,
        uniform_field_load,
    )
    from radia.topology_optimization import solve_hdiv_mmm_active_elements
    from radia.vim._vim import build_charge_gram

    ng.SetNumThreads(args.threads)
    mesh = MakeStructured3DMesh(
        hexes=True, nx=2, ny=1, nz=1,
        mapping=lambda x, y, z: (
            -0.30 + 0.60*x, -0.10 + 0.20*y, 0.05 + 0.20*z))
    fes = ng.HDiv(mesh, order=1, discontinuous=True)
    with ng.TaskManager():
        _, gram, _ = build_charge_gram(
            fes, eps=1.0e-10, leafsize=256, eta=2.0,
            internal_interfaces=True)
        source = uniform_field_load(fes, (0.0, 0.0, args.source_h_a_per_m))
    rhs = np.asarray(source.vec.FV().NumPy(), dtype=float).copy()
    initial_active = np.array([True, False])
    target_active = np.array([True, True])
    zero_response = np.zeros((1, fes.ndof))
    solve_options = dict(
        charge_gram=gram, fes=fes, inv_chi=1.0/(args.mu_r-1.0),
        rhs=rhs, response_matrix=zero_response,
        solve_tolerance=1.0e-11)
    with ng.TaskManager():
        initial_state = solve_hdiv_mmm_active_elements(
            active_elements=initial_active, **solve_options)[0]
        target_state = solve_hdiv_mmm_active_elements(
            active_elements=target_active, **solve_options)[0]

    entries = static_magnet_transfer_component_entries(
        CONTROLLED_COMPONENTS)
    provisional_orbit = _one_segment_arc()
    provisional_objective = MultiMomentumTransferMatrixObjective(
        (provisional_orbit,), np.eye(6)[None, :, :], 1.0, 1.0,
        response_entries=entries)
    with ng.TaskManager():
        response_matrix = build_multi_orbit_field_response_matrix(
            gram, provisional_objective, gradient_offset=args.gradient_offset,
            field_scale=MU0)
    incident = np.array([MU0*args.source_h_a_per_m, 0.0])
    initial_raw = response_matrix@initial_state + incident
    target_raw = response_matrix@target_state + incident
    rigidity = float(
        target_raw[0]/provisional_orbit.signed_curvature[0])
    orbit = PlanarDesignOrbit(
        provisional_orbit.positions, provisional_orbit.tangents,
        magnetic_rigidity=rigidity,
        bend_axis=provisional_orbit.bend_axis)
    target_map = combined_function_transfer_map_from_field_response(
        target_raw, orbit.segment_lengths, orbit.magnetic_rigidity,
        response_entries=entries).matrix
    initial_map = combined_function_transfer_map_from_field_response(
        initial_raw, orbit.segment_lengths, orbit.magnetic_rigidity,
        response_entries=entries).matrix
    selected_change = max(
        abs(target_map[row, column]-initial_map[row, column])
        for row, column in entries)
    bend_change = abs(float(target_raw[0]-initial_raw[0]))
    if selected_change <= 0.0 or bend_change <= 0.0:
        raise RuntimeError("manufactured pole block produced no response")
    transfer_band = args.relative_band*selected_change
    bend_band = args.relative_band*bend_change
    family = build_ffag_fixed_design_orbit_target_family(
        (orbit,), target_map[None, :, :],
        transfer_matrix_band=transfer_band,
        bend_field_band=bend_band,
        controlled_components=CONTROLLED_COMPONENTS)
    objective = family.objective
    field_correction = solve_transfer_matrix_field_correction(
        objective, initial_raw, relative_tolerance=1.0e-10)
    volumes = np.asarray(
        ng.Integrate(1.0, mesh, element_wise=True), dtype=float)
    with ng.TaskManager():
        result = optimize_hdiv_mmm_magnet_from_transfer_matrices(
            (orbit,), target_map[None, :, :],
            transfer_matrix_band=transfer_band,
            bend_field_band=bend_band,
            charge_gram=gram, fes=fes, inv_chi=1.0/(args.mu_r-1.0),
            rhs=rhs, field_response_matrix=response_matrix,
            incident_field_response=incident,
            field_correction=field_correction,
            active_elements=initial_active, element_volumes=volumes,
            volume_max=float(np.sum(volumes))+1.0e-14,
            fixed_active_elements=initial_active,
            response_entries=entries,
            maximum_batch_elements=1,
            graph_front_proposal_limit=0,
            max_iterations=1,
            solve_tolerance=1.0e-11)

    gates = {
        "named_four_component_objective": (
            objective.response_entries
            == ((1, 0), (3, 2), (0, 5), (1, 5))),
        "target_map_is_symplectic": bool(
            static_magnet_symplectic_residual(target_map) < 1.0e-12),
        "known_pole_block_recovered": bool(np.array_equal(
            result.active_elements, target_active)),
        "selected_transfer_components_reached": bool(
            np.max(result.transfer_matrix_max_band_ratios) <= 1.0),
        "design_orbit_bend_field_reached": bool(
            np.max(result.orbit_field_max_band_ratios) <= 1.0),
        "binary_connected_topology": bool(result.topology.valid),
        "realized_map_is_symplectic": bool(
            np.max(result.realized_symplectic_residuals) < 1.0e-12),
    }
    report = {
        "schema": "radia.selected-pole-transfer-components/v1",
        "status": "pass" if all(gates.values()) else "fail",
        "machine": platform.node(),
        "scope": (
            "Manufactured two-HEX BDM1 pole growth with a fixed one-pass "
            "design orbit and four selected transfer-map components."),
        "performance_measurement": "disabled",
        "mesh": {
            "element_family": "HEX",
            "hdiv_family": "BDM1",
            "elements": int(mesh.ne),
            "dofs": int(fes.ndof),
            "initial_active_elements": np.flatnonzero(
                initial_active).tolist(),
            "target_active_elements": np.flatnonzero(target_active).tolist(),
            "realized_active_elements": np.flatnonzero(
                result.active_elements).tolist(),
        },
        "objective": {
            "controlled_components": list(CONTROLLED_COMPONENTS),
            "response_entries": [list(entry) for entry in entries],
            "target_symplectic_residual": (
                static_magnet_symplectic_residual(target_map)),
            "transfer_matrix_band": transfer_band,
            "bend_field_band": bend_band,
            "target_matrix": target_map.tolist(),
            "initial_matrix": initial_map.tolist(),
            "realized_matrix": result.realized_transfer_matrices[0].tolist(),
        },
        "result": {
            "orbit_field_max_band_ratio": float(
                np.max(result.orbit_field_max_band_ratios)),
            "selected_transfer_max_band_ratio": float(
                np.max(result.transfer_matrix_max_band_ratios)),
            "target_symplectic_residual": float(
                np.max(result.target_symplectic_residuals)),
            "realized_symplectic_residual": float(
                np.max(result.realized_symplectic_residuals)),
            "converged": bool(result.converged),
        },
        "gates": gates,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2)+"\n", encoding="utf-8")
    return report


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--source-h-a-per-m", type=float, default=1.0e5)
    parser.add_argument("--mu-r", type=float, default=1000.0)
    parser.add_argument("--gradient-offset", type=float, default=0.02)
    parser.add_argument("--relative-band", type=float, default=0.2)
    args = parser.parse_args(argv)
    if (args.threads < 1 or args.source_h_a_per_m <= 0.0
            or args.mu_r <= 1.0 or args.gradient_offset <= 0.0
            or not 0.0 < args.relative_band < 1.0):
        parser.error("invalid selected-pole validation settings")
    return args


def main():
    report = run(parse_args())
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report["status"] == "pass" else 1)


if __name__ == "__main__":
    main()
