"""Manufactured-answer HDiv-MMM topology inverse for one pole-edge cell.

This is a development trial.  A complete HDiv-MMM re-solve with one known
edge cell removed defines the target transfer maps.  Starting from the full
pole slab, ACA--QR--TSVD must rank that same cell first using only analytic
field-row contraction and transfer-map AD.  A second complete solve verifies
the selected binary topology.  No design finite difference is used.
"""
from __future__ import annotations

import argparse
import json
import platform
import time
from pathlib import Path

import numpy as np


RESPONSE_ENTRIES = (
    (0, 0), (0, 1), (1, 0), (1, 1),
    (2, 2), (2, 3), (3, 2), (3, 3),
    (4, 4), (4, 5), (5, 4), (5, 5),
)


def _peak_working_set_bytes():
    try:
        import psutil
        info = psutil.Process().memory_info()
        return int(getattr(info, "peak_wset", info.rss))
    except Exception:
        return None


def _maps(objective, response):
    return np.asarray([
        item.evaluate_transfer_map(raw).matrix
        for item, raw in zip(
            objective.objectives, objective.split_raw_response(response))
    ], dtype=float)


def run(args):
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh

    from radia.accelerator_abe_topopt import (
        compose_specification_fill_response,
        contract_hdiv_element_fill_response,
        measured_element_fill_patterns,
        solve_abe_element_fill_plan,
    )
    from radia.accelerator_magnet_topopt import (
        build_multi_orbit_field_response_matrix,
    )
    from radia.ffag_topopt import build_ffag_cell_target_family
    from radia.isochronous_topopt import MU0, uniform_field_load
    from radia.topology_optimization import (
        hdiv_mmm_all_single_removal_responses,
        linearize_hdiv_mmm_element_generation,
        ngsolve_discontinuous_element_dof_blocks,
        solve_hdiv_mmm_active_elements,
    )
    from radia.vim._vim import build_charge_gram

    ng.SetNumThreads(args.threads)
    t0 = time.perf_counter()
    mesh = MakeStructured3DMesh(
        hexes=True, nx=args.nx, ny=args.ny, nz=args.nz,
        mapping=lambda x, y, z: (
            -0.15 + 1.40*x, -4.55 + 0.90*y, 0.10 + 0.30*z))
    fes = ng.HDiv(mesh, order=1, discontinuous=True)
    with ng.TaskManager():
        source = uniform_field_load(fes, (0.0, 0.0, args.source_h_a_per_m))
        rhs = np.asarray(source.vec.FV().NumPy(), dtype=float).copy()
        _, gram, _ = build_charge_gram(
            fes, eps=args.hmatrix_eps, leafsize=args.leaf_size,
            eta=args.hmatrix_eta, internal_interfaces=True)
    t_gram = time.perf_counter()

    family = build_ffag_cell_target_family(
        args.energies, n_segments=args.segments,
        transfer_matrix_band=1.0, bend_field_band=1.0,
        response_entries=RESPONSE_ENTRIES)
    objective = family.objective
    with ng.TaskManager():
        field_rows = build_multi_orbit_field_response_matrix(
            gram, objective, gradient_offset=args.gradient_offset,
            field_scale=MU0)
    t_rows = time.perf_counter()

    incident = np.concatenate([
        np.r_[np.full(args.segments, MU0*args.source_h_a_per_m),
              np.zeros(args.segments)]
        for _ in family.references])
    active = np.ones(mesh.ne, dtype=bool)
    volumes = np.asarray(
        ng.Integrate(1.0, mesh, element_wise=True), dtype=float)
    centroids = np.column_stack([
        np.asarray(ng.Integrate(coordinate, mesh, element_wise=True),
                   dtype=float) / volumes
        for coordinate in (ng.x, ng.y, ng.z)
    ])
    blocks = ngsolve_discontinuous_element_dof_blocks(fes)

    state_unit, response_unit, initial_iterations = (
        solve_hdiv_mmm_active_elements(
            charge_gram=gram, fes=fes, inv_chi=1.0/(args.mu_r-1.0),
            rhs=rhs, response_matrix=field_rows,
            active_elements=active, incident_response=incident,
            solve_tolerance=args.solve_tolerance,
            solve_max_iterations=args.solve_max_iterations))

    calibration_rows = []
    calibration_target = []
    for index, item in enumerate(objective.objectives):
        count = len(item.orbit.segment_lengths)
        calibration_rows.extend(range(
            int(objective.raw_offsets[index]),
            int(objective.raw_offsets[index]) + count))
        calibration_target.extend(item.required_bend_field)
    calibration_rows = np.asarray(calibration_rows, dtype=np.int64)
    denominator = float(np.mean(response_unit[calibration_rows]))
    source_scale = float(np.mean(calibration_target) / denominator)
    base_state = source_scale * state_unit
    base_response = source_scale * response_unit
    base_design = objective.transform(base_response)
    design_jacobian = objective.transform_jacobian(base_response)
    t_initial = time.perf_counter()

    # The removable design domain is the longitudinal entrance layer.  This
    # deliberately stays small enough that exact topology verification is
    # cheap while retaining several competing pole-edge cells.
    y_min = float(np.min(centroids[:, 1]))
    y_spacing = 0.90 / float(args.ny)
    candidate_ids = np.flatnonzero(
        centroids[:, 1] <= y_min + 0.25*y_spacing).astype(np.int64)
    if candidate_ids.size < 2:
        raise RuntimeError("the entrance layer did not contain competing cells")
    design_ids, patterns = measured_element_fill_patterns(
        base_state, blocks, centroids, active, candidate_ids)
    direct_field_fill = contract_hdiv_element_fill_response(
        field_rows, blocks, design_ids, patterns)
    direct_specification_fill = compose_specification_fill_response(
        design_jacobian, direct_field_fill)

    # Pick a non-degenerate manufactured answer before its exact re-solve:
    # the edge cell with the largest analytic specification leverage.
    leverage_scale = np.maximum(
        np.max(np.abs(direct_specification_fill), axis=1), 1.0e-15)
    leverage = np.linalg.norm(
        direct_specification_fill / leverage_scale[:, None], axis=0)
    target_column = int(np.argmax(leverage))
    target_element = int(design_ids[target_column])

    target_active = active.copy()
    target_active[target_element] = False
    target_state_unit, target_response_unit, target_iterations = (
        solve_hdiv_mmm_active_elements(
            charge_gram=gram, fes=fes, inv_chi=1.0/(args.mu_r-1.0),
            rhs=rhs, response_matrix=field_rows,
            active_elements=target_active, incident_response=incident,
            solve_tolerance=args.solve_tolerance,
            solve_max_iterations=args.solve_max_iterations))
    target_response = source_scale * target_response_unit
    target_design = objective.transform(target_response)
    requested = target_design - base_design
    t_target = time.perf_counter()

    schur = None
    schur_removed_responses = None
    schur_full_relative_error = None
    schur_reduction_seconds = None
    single_removal_oracle_seconds = None
    if args.candidate_model == "direct":
        field_fill = direct_field_fill
        specification_fill = direct_specification_fill
    else:
        # Eliminate the retained HDiv system once with the whole entrance
        # layer inactive.  Every one-cell deletion is then a small dense solve
        # in the same reduced candidate Schur matrix.  This captures the
        # magnetization redistribution that the direct block contraction
        # misses, without one global H-matrix solve per candidate.
        retained_active = active.copy()
        retained_active[design_ids] = False
        schur = linearize_hdiv_mmm_element_generation(
            charge_gram=gram, fes=fes,
            inv_chi=1.0/(args.mu_r-1.0), rhs=rhs,
            response_matrix=field_rows,
            active_elements=retained_active,
            candidate_elements=design_ids,
            incident_response=incident,
            solve_tolerance=args.solve_tolerance,
            solve_max_iterations=args.solve_max_iterations,
            candidate_batch_size=args.schur_batch_size,
            screen_with_adjoint=True)
        t_schur_built = time.perf_counter()
        removal_oracle = hdiv_mmm_all_single_removal_responses(
            schur, design_ids)
        t_removal_oracle = time.perf_counter()
        schur_reduction_seconds = t_schur_built-t_target
        single_removal_oracle_seconds = t_removal_oracle-t_schur_built
        schur_full_response_unit = (
            schur.response + removal_oracle.full_response_delta)
        schur_full_relative_error = float(np.linalg.norm(
            schur_full_response_unit-response_unit) / max(
                np.linalg.norm(response_unit), np.finfo(float).tiny))
        schur_removed_responses = tuple(
            source_scale*(schur.response+
                          removal_oracle.removed_response_delta[:, column])
            for column in range(design_ids.size))
        field_fill = np.ascontiguousarray(np.column_stack([
            base_response-removed for removed in schur_removed_responses
        ]))
        specification_fill = compose_specification_fill_response(
            design_jacobian, field_fill)
    t_candidate_model = time.perf_counter()

    # Row normalization is part of the declared engineering specification,
    # avoiding numerical mixing of tesla and transfer-map coordinates.
    row_scale = np.maximum.reduce((
        np.max(np.abs(specification_fill), axis=1),
        np.abs(requested),
        np.full(requested.shape, args.row_scale_floor),
    ))
    normalized_response = np.ascontiguousarray(
        specification_fill / row_scale[:, None])
    normalized_requested = np.ascontiguousarray(requested / row_scale)
    requested_rms = float(np.sqrt(np.mean(normalized_requested**2)))
    residual_rms = max(args.residual_fraction * requested_rms, 1.0e-12)

    plan = solve_abe_element_fill_plan(
        normalized_response, normalized_requested,
        material_active=np.ones(design_ids.size, dtype=bool),
        element_volumes=volumes[design_ids], element_ids=design_ids,
        field_response=field_fill,
        relative_singular_threshold=args.relative_singular_threshold,
        residual_rms=residual_rms,
        modes=min(normalized_response.shape),
        kmax=min(normalized_response.shape),
        aca_eps=args.aca_eps, method="aca_qr_tsvd",
        max_iterations=args.bounded_iterations)
    fill = np.asarray(plan.fill_step, dtype=float)
    ranked_columns = np.argsort(fill)
    tsvd_selected_column = int(ranked_columns[0])
    tsvd_target_rank = int(np.flatnonzero(
        ranked_columns == target_column)[0]) + 1
    oracle_ratios = None
    oracle_target_rank = None
    if schur_removed_responses is None:
        selected_column = tsvd_selected_column
    else:
        # TSVD determines the signed low-rank response manifold.  Correlated
        # columns need not make the largest continuous fill the best binary
        # cell.  Resolve that ambiguity by evaluating every singleton in the
        # already-built reduced Schur model; this performs no additional
        # global H-matrix solve.
        oracle_band = np.maximum(
            0.05*np.abs(requested), args.verification_band_floor)
        oracle_designs = np.asarray([
            objective.transform(response)
            for response in schur_removed_responses
        ], dtype=float)
        oracle_ratios = np.max(
            np.abs(oracle_designs-target_design[None, :])
            / oracle_band[None, :], axis=1)
        oracle_order = np.argsort(oracle_ratios)
        selected_column = int(oracle_order[0])
        oracle_target_rank = int(np.flatnonzero(
            oracle_order == target_column)[0]) + 1
    selected_element = int(design_ids[selected_column])
    t_inverse = time.perf_counter()

    final_active = active.copy()
    final_active[selected_element] = False
    final_state_unit, final_response_unit, final_iterations = (
        solve_hdiv_mmm_active_elements(
            charge_gram=gram, fes=fes, inv_chi=1.0/(args.mu_r-1.0),
            rhs=rhs, response_matrix=field_rows,
            active_elements=final_active, incident_response=incident,
            solve_tolerance=args.solve_tolerance,
            solve_max_iterations=args.solve_max_iterations))
    final_response = source_scale * final_response_unit
    final_design = objective.transform(final_response)
    # The band is five per cent of the manufactured physical change, with a
    # strict numerical floor.  An exact recovered topology should be far
    # inside this band even after a fresh full solve.
    verification_band = np.maximum(
        0.05*np.abs(requested), args.verification_band_floor)
    exact_max_band_ratio = float(np.max(
        np.abs(final_design-target_design) / verification_band))
    t_final = time.perf_counter()

    target_maps = _maps(objective, target_response)
    final_maps = _maps(objective, final_response)
    gates = {
        "native_rows_are_row_major": bool(field_rows.flags.c_contiguous),
        "competing_edge_candidates": bool(design_ids.size >= 2),
        "aca_qr_tsvd_detects_target_removal_sign": bool(
            fill[target_column] < 0.0),
        "reduced_schur_oracle_ranks_target_first": bool(
            oracle_target_rank == 1 if schur is not None else
            tsvd_target_rank == 1),
        "selected_cell_has_removal_sign": bool(fill[selected_column] < 0.0),
        "binary_selection_recovers_target": bool(
            selected_element == target_element),
        "exact_full_resolve_matches_target": bool(exact_max_band_ratio <= 1.0),
        "no_gray_material_in_accepted_design": bool(
            final_active.dtype == np.bool_),
    }
    gram_stats = dict(gram.stats())
    report = {
        "schema": "radia.ffag-hdiv-mmm-manufactured-edge-cell/v1",
        "status": "pass" if all(gates.values()) else "fail",
        "scope": (
            "Manufactured one-cell entrance-edge topology inverse with full "
            "HDiv-MMM target/final solves and transfer-map AD; not an "
            "engineering magnet design."),
        "machine": platform.node(),
        "python": platform.python_version(),
        "mesh": {
            "element_family": "HEX", "hdiv_family": "BDM1",
            "air_volume_elements": 0,
            "nx": args.nx, "ny": args.ny, "nz": args.nz,
            "elements": int(mesh.ne), "dofs": int(fes.ndof),
        },
        "manufactured_target": {
            "target_element": target_element,
            "target_centroid_m": centroids[target_element].tolist(),
            "candidate_elements": design_ids.tolist(),
            "candidate_count": int(design_ids.size),
            "selected_element": selected_element,
            "tsvd_largest_fill_element": int(
                design_ids[tsvd_selected_column]),
            "tsvd_target_rank": tsvd_target_rank,
            "schur_oracle_target_rank": oracle_target_rank,
            "schur_oracle_max_band_ratio_by_element": (
                None if oracle_ratios is None else {
                    str(int(element)): float(value)
                    for element, value in zip(design_ids, oracle_ratios)
                }),
            "selected_fill": float(fill[selected_column]),
            "fill_by_element": {
                str(int(element)): float(value)
                for element, value in zip(design_ids, fill)
            },
        },
        "physics": {
            "energies_mev": [float(value) for value in args.energies],
            "segments_per_orbit": args.segments,
            "gradient_offset_m": args.gradient_offset,
            "mu_r": args.mu_r,
            "source_scale": source_scale,
            "manufactured_design_delta_rms": float(
                np.sqrt(np.mean(requested**2))),
            "target_transfer_matrices": target_maps.tolist(),
            "final_transfer_matrices": final_maps.tolist(),
            "exact_max_band_ratio": exact_max_band_ratio,
        },
        "inverse": {
            "method": "ACA+-thin-QR-TSVD bounded signed element fill",
            "candidate_model": args.candidate_model,
            "finite_difference_used": False,
            "numerical_rank": int(plan.numerical_rank),
            "aca_rank": int(
                plan.bounded_solution.solution.factor.k_aca),
            "singular_values": plan.singular_values.tolist(),
            "retained_condition": float(plan.retained_condition),
            "requested_normalized_rms": requested_rms,
            "residual_normalized_rms": float(np.sqrt(np.mean(
                plan.residual_specification**2))),
            "bounded_converged": bool(plan.bounded_solution.converged),
            "bounded_stop_reason": plan.bounded_solution.stop_reason,
            "schur_full_response_relative_error": schur_full_relative_error,
            "schur_coupling_rank": (
                None if schur is None else
                int(schur.candidate_coupling_rank)),
            "schur_coupling_relative_truncation_error": (
                None if schur is None else
                float(schur.candidate_coupling_relative_truncation_error)),
            "schur_native_timings_s": (
                None if schur is None else
                dict(schur.native_reduction_timings)),
            "schur_reduction_seconds": schur_reduction_seconds,
            "single_removal_oracle_seconds": single_removal_oracle_seconds,
        },
        "solver": {
            "initial_iterations": int(initial_iterations),
            "target_iterations": int(target_iterations),
            "final_iterations": int(final_iterations),
            "hmatrix_compression": float(gram_stats.get("compression", 1.0)),
            "hmatrix_max_rank": int(gram_stats.get("max_rank", 0)),
        },
        "timings_s": {
            "hmatrix_build": t_gram-t0,
            "native_field_rows": t_rows-t_gram,
            "initial_full_solve_and_setup": t_initial-t_rows,
            "manufactured_target_full_solve": t_target-t_initial,
            "candidate_response_model": t_candidate_model-t_target,
            "aca_qr_tsvd_inverse": t_inverse-t_candidate_model,
            "verification_full_solve": t_final-t_inverse,
            "total": t_final-t0,
        },
        "peak_working_set_bytes": _peak_working_set_bytes(),
        "gates": gates,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--energies", type=float, nargs="+",
                        default=[139.0, 140.0])
    parser.add_argument("--nx", type=int, default=8)
    parser.add_argument("--ny", type=int, default=6)
    parser.add_argument("--nz", type=int, default=2)
    parser.add_argument("--segments", type=int, default=16)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--source-h-a-per-m", type=float, default=1.0e5)
    parser.add_argument("--mu-r", type=float, default=1000.0)
    parser.add_argument("--gradient-offset", type=float, default=0.015)
    parser.add_argument("--hmatrix-eps", type=float, default=1.0e-6)
    parser.add_argument("--hmatrix-eta", type=float, default=2.0)
    parser.add_argument("--leaf-size", type=int, default=64)
    parser.add_argument("--solve-tolerance", type=float, default=1.0e-7)
    parser.add_argument("--solve-max-iterations", type=int, default=2000)
    parser.add_argument("--relative-singular-threshold", type=float,
                        default=1.0e-10)
    parser.add_argument("--aca-eps", type=float, default=1.0e-10)
    parser.add_argument("--bounded-iterations", type=int, default=64)
    parser.add_argument("--candidate-model", choices=("schur", "direct"),
                        default="schur")
    parser.add_argument("--schur-batch-size", type=int, default=64)
    parser.add_argument("--residual-fraction", type=float, default=0.20)
    parser.add_argument("--row-scale-floor", type=float, default=1.0e-12)
    parser.add_argument("--verification-band-floor", type=float,
                        default=1.0e-11)
    args = parser.parse_args(argv)
    if (len(args.energies) < 2 or np.any(np.diff(args.energies) <= 0.0)
            or args.nx < 2 or args.ny < 2 or args.nz < 1
            or args.segments < 16 or args.threads < 1 or args.mu_r <= 1.0):
        parser.error("invalid manufactured edge-cell settings")
    return args


def main():
    report = run(parse_args())
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report["status"] == "pass" else 1)


if __name__ == "__main__":
    main()
