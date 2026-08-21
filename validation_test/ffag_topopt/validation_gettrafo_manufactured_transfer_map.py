"""Manufactured transfer-map inverse for one smooth GetTrafo pole mode.

A known topology-preserving deformation of a BDM1 TET pole creates the target
through a complete HDiv-MMM solve.  Starting from the undeformed pole, the
production analytic ``dM/dq,dB/dq,dG/dq,dC/dq,drhs/dq`` chain supplies each
relinearized shape LP.  Every accepted coefficient is rebuilt and completely
re-solved; a final independent solve verifies the manufactured transfer maps.
No design finite difference is used.

This is deliberately the straight-TET proof of the shape formulation.  The
configured-field ``dC/dq`` native kernel does not yet expose its HEX analogue;
the earlier whole-element material stage remains BDM1 HEX.
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


def _transfer_matrices(objective, raw_response):
    return np.asarray([
        item.evaluate_transfer_map(values).matrix
        for item, values in zip(
            objective.objectives,
            objective.split_raw_response(raw_response),
        )
    ], dtype=float)


def run(args):
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh

    from radia.accelerator_magnet_topopt import (
        build_multi_orbit_field_response_matrix,
        multi_orbit_field_observations,
    )
    from radia.ffag_topopt import build_ffag_cell_target_family
    from radia.isochronous_topopt import MU0
    from radia.topology_optimization import (
        ShapeLinearization,
        assemble_ngsolve_hdiv_linear_form_shape_tangents,
        production_vim_functional_shape_jacobian_streaming,
        sample_production_gettrafo_displacements,
        solve_shape_lp,
    )
    from radia.vim._vim import _charge_basis, build_charge_gram

    ng.SetNumThreads(args.threads)
    started = time.perf_counter()
    mesh = MakeStructured3DMesh(
        hexes=False, nx=args.nx, ny=args.ny, nz=args.nz,
        mapping=lambda x, y, z: (
            -0.15 + 1.40*x, -4.55 + 0.90*y, 0.10 + 0.30*z))
    fes = ng.HDiv(mesh, order=1)
    mode_space = ng.VectorH1(mesh, order=1)
    mode = ng.GridFunction(mode_space)
    xi = (ng.x + 0.15) / 1.40
    yi = (ng.y + 4.55) / 0.90
    zi = (ng.z - 0.10) / 0.30
    # The upper attachment and four lateral sides remain fixed.  Only the
    # aperture-side pole face moves, with a smooth separable envelope.
    envelope = (4.0*xi*(1.0-xi))*(4.0*yi*(1.0-yi))*(1.0-zi)
    mode.Set(ng.CF((0.0, 0.0, -envelope)))
    deformation = ng.GridFunction(mode_space)

    family = build_ffag_cell_target_family(
        args.energies, n_segments=args.segments,
        transfer_matrix_band=1.0, bend_field_band=1.0,
        response_entries=RESPONSE_ENTRIES)
    objective = family.objective
    observations, response_weights = multi_orbit_field_observations(
        objective, gradient_offset=args.gradient_offset,
        field_scale=MU0)
    incident = np.concatenate([
        np.r_[np.full(args.segments, MU0*args.source_h_a_per_m),
              np.zeros(args.segments)]
        for _ in objective.orbits
    ])
    applied = ng.CF((0.0, 0.0, args.source_h_a_per_m))
    inv_chi = 1.0 / (args.mu_r - 1.0)

    solve_timings = []

    def solve_at(amplitude_m, *, derivative, check_batched_rows=False):
        call_started = time.perf_counter()
        deformation.vec.data = float(amplitude_m) * mode.vec
        mesh.SetDeformation(deformation)
        try:
            with ng.TaskManager():
                # The production geometry sampler consumes this exact
                # GetTrafo-aware charge basis.  It is intentionally not
                # reconstructed from FE shapes in the optimizer.
                basis = (_charge_basis(
                    fes, 3, materialize_mass=False)
                    if derivative else None)
                charge_map, gram, _ = build_charge_gram(
                    fes, eps=args.hmatrix_eps,
                    leafsize=args.leaf_size, eta=args.hmatrix_eta)
                rhs, rhs_jacobian = (
                    assemble_ngsolve_hdiv_linear_form_shape_tangents(
                        fes, applied, (mode,) if derivative else (),
                        bonus_intorder=args.integration_order))
                response_matrix = np.ascontiguousarray(np.asarray(
                    gram.configured_field_functional_rows(
                        observations, response_weights), dtype=float))
                row_parity_error = None
                if check_batched_rows:
                    batched = build_multi_orbit_field_response_matrix(
                        gram, objective,
                        gradient_offset=args.gradient_offset,
                        field_scale=MU0)
                    row_parity_error = float(np.max(np.abs(
                        response_matrix - batched)))
                if derivative:
                    linearized = (
                        production_vim_functional_shape_jacobian_streaming(
                            fes=fes, deformation_modes=(mode,),
                            charge_basis=basis, charge_gram=gram,
                            charge_map=charge_map, inv_chi=inv_chi,
                            rhs=rhs, response_matrix=response_matrix,
                            rhs_jacobian=rhs_jacobian,
                            response_observations=observations,
                            response_weights=response_weights,
                            family="tet", incident_response=incident,
                            solve_tolerance=args.solve_tolerance,
                            solve_max_iterations=args.solve_max_iterations,
                            mass_riesz=True))
                    state = np.asarray(linearized.state, dtype=float)
                    raw_response = np.asarray(
                        linearized.response, dtype=float)
                    raw_jacobian = np.asarray(
                        linearized.response_jacobian, dtype=float)
                    iterations = int(linearized.state_iterations)
                    derivative_timings = dict(linearized.timings_s or {})
                    geometry = sample_production_gettrafo_displacements(
                        fes, (mode,), basis, family="tet")
                    maximum_mode_value = float(max(
                        np.max(np.linalg.norm(values[0], axis=1))
                        for values in geometry.cell))
                else:
                    solved = (
                        gram.solve_configured_linear_material_auto_prec_many(
                            inv_chi, np.ascontiguousarray(rhs[None, :]),
                            tol=args.solve_tolerance,
                            maxit=args.solve_max_iterations,
                            mass_riesz=True))
                    state = np.asarray(solved["m"], dtype=float)[0]
                    raw_response = response_matrix @ state + incident
                    raw_jacobian = None
                    iterations = int(solved["iters"][0])
                    derivative_timings = None
                    maximum_mode_value = None
                record = {
                    "amplitude_m": float(amplitude_m),
                    "derivative": bool(derivative),
                    "wall_s": time.perf_counter() - call_started,
                    "state_iterations": iterations,
                    "analytic_derivative_timings_s": derivative_timings,
                }
                solve_timings.append(record)
                return {
                    "state": state,
                    "raw_response": raw_response,
                    "raw_jacobian": raw_jacobian,
                    "iterations": iterations,
                    "response_matrix_row_major": bool(
                        response_matrix.flags.c_contiguous),
                    "row_parity_error": row_parity_error,
                    "gram_stats": dict(gram.stats()),
                    "maximum_mode_value": maximum_mode_value,
                }
        finally:
            mesh.UnsetDeformation()

    target_amplitude = args.target_amplitude_mm * 1.0e-3
    move_limit = args.move_limit_mm * 1.0e-3
    parameter_lower = args.parameter_min_mm * 1.0e-3
    parameter_upper = args.parameter_max_mm * 1.0e-3

    target_solve = solve_at(target_amplitude, derivative=False)
    target_raw = target_solve["raw_response"]
    target_design = objective.transform(target_raw)
    target_maps = _transfer_matrices(objective, target_raw)
    target_done = time.perf_counter()

    base_solve = solve_at(
        0.0, derivative=False, check_batched_rows=True)
    base_raw = base_solve["raw_response"]
    base_design = objective.transform(base_raw)
    manufactured_change = target_design - base_design
    response_band = np.maximum.reduce((
        args.verification_relative_band*np.abs(manufactured_change),
        args.verification_target_relative_floor*np.abs(target_design),
        np.full(target_design.shape, args.verification_band_floor),
    ))
    base_done = time.perf_counter()

    amplitude = 0.0
    history = []
    converged = False
    initial_linearized_base_error = None
    maximum_mode_value = None
    for iteration in range(args.max_iterations):
        analytic = solve_at(amplitude, derivative=True)
        current_raw = analytic["raw_response"]
        current_design = objective.transform(current_raw)
        design_jacobian = (
            objective.transform_jacobian(current_raw)
            @ analytic["raw_jacobian"])
        normalized_residual = (
            (current_design-target_design) / response_band)
        maximum_ratio = float(np.max(np.abs(normalized_residual)))
        maximum_mode_value = analytic["maximum_mode_value"]
        if iteration == 0:
            initial_linearized_base_error = float(np.max(np.abs(
                current_raw-base_raw)))
        entry = {
            "iteration": iteration,
            "amplitude_m": amplitude,
            "exact_max_band_ratio": maximum_ratio,
            "state_iterations": analytic["iterations"],
            "analytic_derivative_timings_s": solve_timings[-1][
                "analytic_derivative_timings_s"],
        }
        if maximum_ratio <= 1.0:
            entry.update({
                "accepted_delta_m": 0.0,
                "predicted_max_band_ratio": maximum_ratio,
                "lp_status": "already inside manufactured bands",
            })
            history.append(entry)
            converged = True
            break
        objective_value = 0.5*float(
            normalized_residual @ normalized_residual)
        objective_gradient = np.asarray([
            float(design_jacobian[:, 0] @ (
                (current_design-target_design) / response_band**2))
        ])
        linearization = ShapeLinearization(
            objective_value, objective_gradient,
            current_design, design_jacobian,
            target_design, response_band)
        update = solve_shape_lp(
            np.asarray([amplitude]), linearization,
            move_limit=np.asarray([move_limit]),
            parameter_bounds=(
                np.asarray([parameter_lower]),
                np.asarray([parameter_upper])))
        amplitude = float(update.parameters[0])
        entry.update({
            "accepted_delta_m": float(update.delta[0]),
            "predicted_max_band_ratio": float(
                update.predicted_max_band_ratio),
            "lp_status": update.status,
        })
        history.append(entry)
        if abs(float(update.delta[0])) <= args.parameter_tolerance_mm*1.0e-3:
            break
    inverse_done = time.perf_counter()

    # This solve is independent of the state returned by the final analytic
    # call and is the acceptance authority for the manufactured target.
    final_solve = solve_at(
        amplitude, derivative=False, check_batched_rows=True)
    final_raw = final_solve["raw_response"]
    final_design = objective.transform(final_raw)
    final_maps = _transfer_matrices(objective, final_raw)
    final_ratio = float(np.max(np.abs(
        (final_design-target_design) / response_band)))
    amplitude_relative_error = float(
        abs(amplitude-target_amplitude) / abs(target_amplitude))
    finished = time.perf_counter()

    target_change_max = float(np.max(np.abs(manufactured_change)))
    first_delta = history[0]["accepted_delta_m"] if history else 0.0
    row_parity_error = max(
        float(base_solve["row_parity_error"]),
        float(final_solve["row_parity_error"]))
    gates = {
        "bdm1_tet_shape_lane": bool(
            fes.globalorder == 1
            and all(len(element.vertices) == 4
                    for element in mesh.Elements(ng.VOL))),
        "fixed_topology": bool(mesh.ne == 6*args.nx*args.ny*args.nz),
        "manufactured_map_change_is_nontrivial": bool(
            target_change_max >= args.minimum_manufactured_change),
        "packed_native_rows_are_row_major": bool(
            observations.flags.c_contiguous
            and response_weights.flags.c_contiguous
            and base_solve["response_matrix_row_major"]
            and final_solve["response_matrix_row_major"]),
        "packed_rows_match_batched_rows": bool(
            row_parity_error <= args.row_parity_tolerance),
        "analytic_linearization_reproduces_base_solve": bool(
            initial_linearized_base_error is not None
            and initial_linearized_base_error
            <= args.solve_consistency_tolerance),
        "first_lp_step_has_target_sign": bool(
            first_delta*target_amplitude > 0.0),
        "shape_lp_enters_all_transfer_map_bands": bool(
            converged and final_ratio <= 1.0),
        "manufactured_amplitude_is_recovered": bool(
            amplitude_relative_error
            <= args.amplitude_relative_tolerance),
        "finite_difference_not_used": True,
    }
    gram_stats = final_solve["gram_stats"]
    report = {
        "schema": "radia.ffag-hdiv-mmm-manufactured-gettrafo-map/v1",
        "status": "pass" if all(gates.values()) else "fail",
        "scope": (
            "One smooth topology-preserving BDM1-TET pole-face mode, "
            "analytic HDiv-MMM shape derivative, relinearized shape LP, "
            "and full transfer-map re-solves; not an engineering magnet "
            "design and not the pending HEX configured-field derivative."),
        "machine": platform.node(),
        "python": platform.python_version(),
        "mesh": {
            "element_family": "TET",
            "hdiv_family": "BDM1",
            "air_volume_elements": 0,
            "nx": args.nx,
            "ny": args.ny,
            "nz": args.nz,
            "elements": int(mesh.ne),
            "dofs": int(fes.ndof),
            "topology_changed": False,
        },
        "manufactured_target": {
            "mode": (
                "V=(0,0,-4*xi*(1-xi)*4*yi*(1-yi)*(1-zi))"),
            "target_amplitude_m": target_amplitude,
            "recovered_amplitude_m": amplitude,
            "amplitude_relative_error": amplitude_relative_error,
            "maximum_mode_vertex_norm": maximum_mode_value,
            "maximum_target_vertex_displacement_m": (
                None if maximum_mode_value is None
                else abs(target_amplitude)*maximum_mode_value),
            "maximum_design_response_change": target_change_max,
            "target_transfer_matrices": target_maps.tolist(),
            "final_transfer_matrices": final_maps.tolist(),
            "final_exact_max_band_ratio": final_ratio,
        },
        "physics": {
            "energies_mev": [float(value) for value in args.energies],
            "segments_per_orbit": args.segments,
            "gradient_offset_m": args.gradient_offset,
            "mu_r": args.mu_r,
            "source_h_a_per_m": args.source_h_a_per_m,
        },
        "inverse": {
            "method": (
                "relinearized topology-preserving Chebyshev shape LP"),
            "derivative": (
                "analytic GetTrafo dM+dB+dG+dC+drhs, matrix-free "
                "H-matrix adjoint contractions"),
            "finite_difference_used": False,
            "iterations": history,
            "response_band_relative_to_manufactured_change": (
                args.verification_relative_band),
            "row_parity_max_abs_error": row_parity_error,
            "fresh_full_solve_max_band_ratio": final_ratio,
        },
        "solver": {
            "hmatrix_eps": args.hmatrix_eps,
            "hmatrix_compression": float(
                gram_stats.get("compression", 1.0)),
            "hmatrix_max_rank": int(gram_stats.get("max_rank", 0)),
            "target_iterations": target_solve["iterations"],
            "base_iterations": base_solve["iterations"],
            "final_iterations": final_solve["iterations"],
        },
        "timings_s": {
            "manufactured_target_full_solve": target_done-started,
            "base_full_solve_and_row_parity": base_done-target_done,
            "analytic_shape_lp_iterations": inverse_done-base_done,
            "fresh_verification_full_solve": finished-inverse_done,
            "total": finished-started,
            "per_solve": solve_timings,
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
    parser.add_argument("--nx", type=int, default=4)
    parser.add_argument("--ny", type=int, default=3)
    parser.add_argument("--nz", type=int, default=1)
    parser.add_argument("--segments", type=int, default=16)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--source-h-a-per-m", type=float, default=1.0e5)
    parser.add_argument("--mu-r", type=float, default=1000.0)
    parser.add_argument("--gradient-offset", type=float, default=0.015)
    parser.add_argument("--target-amplitude-mm", type=float, default=4.0)
    parser.add_argument("--move-limit-mm", type=float, default=5.0)
    parser.add_argument("--parameter-min-mm", type=float, default=-20.0)
    parser.add_argument("--parameter-max-mm", type=float, default=20.0)
    parser.add_argument("--parameter-tolerance-mm", type=float,
                        default=1.0e-6)
    parser.add_argument("--max-iterations", type=int, default=6)
    parser.add_argument("--integration-order", type=int, default=4)
    parser.add_argument("--hmatrix-eps", type=float, default=1.0e-7)
    parser.add_argument("--hmatrix-eta", type=float, default=2.0)
    parser.add_argument("--leaf-size", type=int, default=64)
    parser.add_argument("--solve-tolerance", type=float, default=1.0e-10)
    parser.add_argument("--solve-max-iterations", type=int, default=5000)
    parser.add_argument("--verification-relative-band", type=float,
                        default=1.0e-4)
    parser.add_argument("--verification-target-relative-floor", type=float,
                        default=1.0e-10)
    parser.add_argument("--verification-band-floor", type=float,
                        default=1.0e-12)
    parser.add_argument("--minimum-manufactured-change", type=float,
                        default=1.0e-6)
    parser.add_argument("--row-parity-tolerance", type=float,
                        default=1.0e-13)
    parser.add_argument("--solve-consistency-tolerance", type=float,
                        default=1.0e-9)
    parser.add_argument("--amplitude-relative-tolerance", type=float,
                        default=2.0e-4)
    args = parser.parse_args(argv)
    target_m = args.target_amplitude_mm*1.0e-3
    if (len(args.energies) < 2 or np.any(np.diff(args.energies) <= 0.0)
            or args.nx < 2 or args.ny < 2 or args.nz < 1
            or args.segments < 16 or args.threads < 1
            or args.mu_r <= 1.0 or args.source_h_a_per_m == 0.0
            or target_m == 0.0 or args.move_limit_mm <= 0.0
            or not (args.parameter_min_mm*1.0e-3 < target_m
                    < args.parameter_max_mm*1.0e-3)
            or args.max_iterations < 1 or args.integration_order < 1
            or args.verification_relative_band <= 0.0
            or args.verification_target_relative_floor <= 0.0
            or args.verification_band_floor <= 0.0
            or args.amplitude_relative_tolerance <= 0.0):
        parser.error("invalid manufactured GetTrafo settings")
    return args


def main():
    report = run(parse_args())
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report["status"] == "pass" else 1)


if __name__ == "__main__":
    main()
