"""Study-scale FFAG multi-momentum BDM1 HDiv-MMM topology PoC.

The candidate domain is a one-sided structured HEX pole slab bordering the
unmeshed orbit aperture.  A uniform applied source excites the iron; native
Laplace field rows sample every orbit in vacuum.  The purpose of this lane is
to demonstrate that one exact whole-element add/remove iteration improves the
Bell--Abell soft-edge multi-momentum bend/map objective.  It is not a final
two-pole/yoke/coil engineering design.
"""
from __future__ import annotations

import argparse
import json
import platform
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np


def _peak_working_set_bytes():
    try:
        import psutil
        info = psutil.Process().memory_info()
        return int(getattr(info, "peak_wset", info.rss))
    except Exception:
        return None


def _history_record(item):
    record = asdict(item)
    for key, value in tuple(record.items()):
        if isinstance(value, np.ndarray):
            if key == "linearized_reachability_residual":
                record[key] = None
            else:
                record[key] = value.tolist()
    return record


def run(args):
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh

    from radia.accelerator_magnet_topopt import (
        build_multi_orbit_field_response_matrix,
        optimize_hdiv_mmm_magnet_from_transfer_matrices,
    )
    from radia.ffag_topopt import build_ffag_cell_target_family
    from radia.isochronous_topopt import MU0, uniform_field_load
    from radia.topology_optimization import solve_hdiv_mmm_active_elements
    from radia.vim._vim import build_charge_gram

    ng.SetNumThreads(args.threads)
    started = time.perf_counter()
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
    hmatrix_done = time.perf_counter()

    family = build_ffag_cell_target_family(
        args.energies, n_segments=args.segments,
        transfer_matrix_band=args.transfer_band,
        bend_field_band=args.bend_band)
    with ng.TaskManager():
        field_rows = build_multi_orbit_field_response_matrix(
            gram, family.objective, gradient_offset=args.gradient_offset,
            field_scale=MU0)
    rows_done = time.perf_counter()

    incident = np.concatenate([
        np.r_[np.full(args.segments, MU0*args.source_h_a_per_m),
              np.zeros(args.segments)]
        for _ in family.references])
    calibration_rows = []
    calibration_target = []
    offsets = family.objective.raw_offsets
    for index, objective in enumerate(family.objective.objectives):
        count = len(objective.orbit.segment_lengths)
        calibration_rows.extend(range(int(offsets[index]),
                                      int(offsets[index]) + count))
        calibration_target.extend(objective.required_bend_field)
    calibration_rows = np.asarray(calibration_rows, dtype=np.int64)
    calibration_target = np.asarray(calibration_target, dtype=float)

    active = np.ones(mesh.ne, dtype=bool)
    fixed_active = np.zeros(mesh.ne, dtype=bool)
    fixed_active[0] = True
    volumes = np.asarray(
        ng.Integrate(1.0, mesh, element_wise=True), dtype=float)
    initial_state, initial_response, initial_solve_iterations = (
        solve_hdiv_mmm_active_elements(
            charge_gram=gram, fes=fes, inv_chi=1.0/(args.mu_r-1.0),
            rhs=rhs, response_matrix=field_rows,
            active_elements=active, incident_response=incident,
            solve_tolerance=args.solve_tolerance,
            solve_max_iterations=args.solve_max_iterations))
    denominator = float(np.mean(initial_response[calibration_rows]))
    source_scale = float(np.mean(calibration_target) / denominator)
    initial_response = initial_response * source_scale
    initial_state = initial_state * source_scale
    initial_objective = family.objective.transform(initial_response)
    initial_ratio = float(np.max(np.abs(
        (initial_objective - family.objective.response_target)
        / family.objective.response_band)))
    initial_done = time.perf_counter()

    result = optimize_hdiv_mmm_magnet_from_transfer_matrices(
        tuple(reference.orbit for reference in family.references),
        np.asarray([reference.transfer.matrix
                    for reference in family.references]),
        transfer_matrix_band=args.transfer_band,
        bend_field_band=args.bend_band,
        charge_gram=gram, fes=fes, inv_chi=1.0/(args.mu_r-1.0), rhs=rhs,
        field_response_matrix=field_rows,
        incident_field_response=incident,
        active_elements=active, element_volumes=volumes,
        volume_max=float(np.sum(volumes)),
        fixed_active_elements=fixed_active,
        max_iterations=args.iterations,
        solve_tolerance=args.solve_tolerance,
        solve_max_iterations=args.solve_max_iterations,
        source_calibration_rows=calibration_rows,
        source_calibration_target=calibration_target,
        initial_material_move_fraction=args.move_fraction,
        maximum_material_move_fraction=args.maximum_move_fraction,
        proposal_adjoint_count=args.proposal_adjoint_count,
        graph_front_proposal_limit=0,
        exact_candidate_limit=args.exact_candidate_limit)
    finished = time.perf_counter()
    final_ratio = float(max(
        np.max(result.orbit_field_max_band_ratios),
        np.max(result.transfer_matrix_max_band_ratios)))
    accepted = len(result.generation.history)
    gates = {
        "native_rows_are_row_major": bool(field_rows.flags.c_contiguous),
        "binary_topology_is_valid": bool(result.topology.valid),
        "accepted_whole_element_move": accepted > 0,
        "exact_resolve_improves_objective": final_ratio < initial_ratio,
        "no_gray_material": bool(result.active_elements.dtype == np.bool_),
    }
    gram_stats = dict(gram.stats())
    report = {
        "schema": "radia.ffag-hdiv-mmm-poc/v1",
        "status": "pass" if all(gates.values()) else "fail",
        "scope": (
            "One-sided HEX pole-slab PoC with a uniform applied source; "
            "not a final two-pole yoke or engineered FFAG coil."),
        "machine": platform.node(),
        "python": platform.python_version(),
        "ngsolve_threads": args.threads,
        "energies_mev": [float(value) for value in args.energies],
        "mesh": {
            "element_family": "HEX",
            "hdiv_family": "BDM1",
            "air_volume_elements": 0,
            "nx": args.nx, "ny": args.ny, "nz": args.nz,
            "elements": int(mesh.ne), "dofs": int(fes.ndof),
            "initial_active_elements": int(active.sum()),
            "final_active_elements": int(result.active_elements.sum()),
        },
        "field_contract": {
            "raw_rows": int(field_rows.shape[0]),
            "design_responses": int(family.objective.response_target.size),
            "segments_per_orbit": args.segments,
            "gradient_offset_m": args.gradient_offset,
            "row_major": bool(field_rows.flags.c_contiguous),
        },
        "solver": {
            "mu_r": args.mu_r,
            "hmatrix_eps": args.hmatrix_eps,
            "hmatrix_eta": args.hmatrix_eta,
            "leaf_size": args.leaf_size,
            "solve_tolerance": args.solve_tolerance,
            "solve_max_iterations": args.solve_max_iterations,
            "initial_solve_iterations": int(initial_solve_iterations),
            "source_scale_initial": source_scale,
            "source_scale_final": result.generation.source_scale,
            "hmatrix_compression": float(gram_stats.get("compression", 1.0)),
            "hmatrix_max_rank": int(gram_stats.get("max_rank", 0)),
            "hmatrix_lowrank_blocks": int(gram_stats.get("n_lowrank", 0)),
        },
        "optimization": {
            "requested_iterations": args.iterations,
            "accepted_iterations": accepted,
            "converged": bool(result.converged),
            "stop_reason": result.generation.stop_reason,
            "initial_max_band_ratio": initial_ratio,
            "final_max_band_ratio": final_ratio,
            "relative_improvement": float(
                (initial_ratio-final_ratio)/initial_ratio),
            "orbit_field_max_band_ratios": (
                result.orbit_field_max_band_ratios.tolist()),
            "transfer_matrix_max_band_ratios": (
                result.transfer_matrix_max_band_ratios.tolist()),
            "history": [_history_record(item)
                        for item in result.generation.history],
        },
        "timings_s": {
            "hmatrix_build": hmatrix_done-started,
            "native_field_rows": rows_done-hmatrix_done,
            "initial_exact_solve": initial_done-rows_done,
            "optimization": finished-initial_done,
            "total": finished-started,
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
                        default=[31.0, 140.0, 250.0])
    parser.add_argument("--nx", type=int, default=12)
    parser.add_argument("--ny", type=int, default=8)
    parser.add_argument("--nz", type=int, default=2)
    parser.add_argument("--segments", type=int, default=32)
    parser.add_argument("--iterations", type=int, default=2)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--source-h-a-per-m", type=float, default=1.0e5)
    parser.add_argument("--mu-r", type=float, default=1000.0)
    parser.add_argument("--gradient-offset", type=float, default=0.015)
    parser.add_argument("--bend-band", type=float, default=0.5)
    parser.add_argument("--transfer-band", type=float, default=0.3)
    parser.add_argument("--hmatrix-eps", type=float, default=1.0e-6)
    parser.add_argument("--hmatrix-eta", type=float, default=2.0)
    parser.add_argument("--leaf-size", type=int, default=64)
    parser.add_argument("--solve-tolerance", type=float, default=1.0e-7)
    parser.add_argument("--solve-max-iterations", type=int, default=2000)
    parser.add_argument("--move-fraction", type=float, default=0.05)
    parser.add_argument("--maximum-move-fraction", type=float, default=0.20)
    parser.add_argument("--proposal-adjoint-count", type=int, default=6)
    parser.add_argument("--exact-candidate-limit", type=int, default=64)
    args = parser.parse_args(argv)
    if (len(args.energies) < 2 or args.nx < 2 or args.ny < 2 or args.nz < 1
            or args.segments < 16 or args.iterations < 1
            or args.threads < 1 or args.mu_r <= 1.0):
        parser.error("invalid FFAG HDiv-MMM PoC settings")
    return args


def main():
    report = run(parse_args())
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report["status"] == "pass" else 1)


if __name__ == "__main__":
    main()
