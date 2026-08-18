"""Multi-orbit three-route test (A-LIE / A-RK / B-RK) on the C-type magnet.

Several parallel-entry design orbits (transverse global-y offsets) sample
different horizontal slices of the same solved field; the certified region
is limited vertically to 80% of the half-gap (|y| <= 4 mm of the 5 mm
half-gap).  For EACH orbit the full pipeline runs in that orbit's own
Bishop frame: graded CanonicalHCurl chain, full-volume B fit with honesty
certificates, element-spoly Lie map, canonical A-RK on the same
polynomials, and the independent mechanical Cartesian B-RK -- so the
three-route agreement is measured per design orbit.

One HDiv-MMM solve is shared by all orbits.  Momentum-family orbits are
nearly degenerate for this short weak magnet (sagitta ~0.1 mm), so the
family is spanned by entrance offsets instead.

Usage (LAB smoke; heavy sweeps go to mdx/hibino):
  python validation_canonical_hcurl_multiorbit.py --orbit-offsets-mm -8 -4 0 4 8
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from datetime import UTC

import ngsolve as ng
from validation_canonical_hcurl_ctype import (
    frame_axes,
    sample_frame_cloud,
    track_chain_a_rk,
    track_mechanical_b_rk,
)
from validation_earlytimes_ctype_ab import (
    build_coil,
    build_iron,
    load_bh_table,
    make_symmetric_b_field,
    track_reference_orbit,
)

import radia as rad
from radia import vim
from radia.accelerator_lie_topopt import (
    _fourth_order_lie_map_from_vector_potential_polynomials,
    apply_dragt_finn_map,
)
from radia.beam_canonical_hcurl import (
    CanonicalHCurlChain,
    graded_breaks,
)

MU0 = 4.0e-7 * np.pi


def make_symmetric_b_batch_tree(solution, coil, tree_options):
    """Symmetrized B batch on the RELAXED Barnes-Hut tree.

    Used for the fit/audit/monitor clouds ONLY -- the independent B-RK
    route and the orbit tracker stay on the exact direct evaluator.  The
    caller must verify the tree against direct on a probe cloud (the
    runner gates it) because the relaxed theta trades ~5e-4 relative
    accuracy for the measured 5x speed.
    """

    def batch(points):
        values = np.ascontiguousarray(
            np.asarray(points, dtype=float).reshape(-1, 3))
        reflected = values.copy()
        reflected[:, 2] *= -1.0
        stacked = np.vstack((values, reflected))
        demag = vim.FieldFromSolution(
            solution, stacked, algorithm="tree", tree_options=tree_options)
        coil_b = np.asarray(rad.Fld(coil, "b", stacked))
        total = MU0 * demag + coil_b
        direct_half = total[:values.shape[0]]
        reflected_half = total[values.shape[0]:] * np.array([-1.0, -1.0, 1.0])
        return 0.5 * (direct_half + reflected_half)

    return batch


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--iron-maxh", type=float, default=0.02)
    result.add_argument("--threads", type=int, default=8)
    result.add_argument("--current", type=float, default=3000.0)
    result.add_argument("--magnetic-rigidity", type=float, default=3.0)
    result.add_argument("--orbit-offsets-mm", type=float, nargs="+",
                        default=[-8.0, -4.0, 0.0, 4.0, 8.0])
    result.add_argument("--half-width", type=float, default=0.010)
    result.add_argument("--half-height", type=float, default=0.004,
                        help="80%% of the 5 mm half-gap by default")
    # Production spec (Sugahara 2026-08-18): decapole content -> p_x = 5
    # (A degree 5, B degree 4); p_s = 2 -- the measured 5-orbit battery is
    # certificate-identical to p_s = 3 (the odd-order fringe advantage was
    # a thick-single-slab artifact), so the C1-quadratic chain suffices.
    result.add_argument("--order-x", type=int, default=5)
    result.add_argument("--order-s", type=int, default=2)
    result.add_argument("--elements", type=int, default=24)
    result.add_argument("--grade-fringe", type=float, default=6.0)
    # Information-minimum budget: the whole magnet is ~p_x*(E+p_s) ~ 200
    # numbers, so a few thousand 3-component samples oversample it several
    # times over -- larger clouds only pay more source evaluation.
    result.add_argument("--fit-stations", type=int, default=24)
    result.add_argument("--fit-points-per-station", type=int, default=40)
    result.add_argument("--lie-degree", type=int, default=5)
    result.add_argument("--reference-orbit-tolerance", type=float,
                        default=5.0e-3)
    result.add_argument("--track-scale", type=float, default=0.5)
    result.add_argument("--rk-tolerance", type=float, default=1.0e-10)
    result.add_argument("--b-rk-tolerance", type=float, default=1.0e-8,
                        help="B-RK integration tolerance; the SOURCE stays "
                             "exact-direct (route purity), only the step "
                             "control relaxes -- 100x under the compared "
                             "5e-6 scale")
    result.add_argument("--tree-theta", type=float, default=0.5)
    result.add_argument("--tree-gate", type=float, default=5.0e-3,
                        help="maximum allowed tree-vs-direct relative "
                             "error on the probe cloud")
    result.add_argument("--output", type=str, default=None)
    return result


def run_orbit(offset_m, options, b_point, b_batch, rng):
    """``b_batch`` here is the verified TREE batch (fit/audit/monitor);
    the exact ``b_point`` keeps the orbit tracker and B-RK source-exact."""
    timings = {}

    def clock(label):
        timings[label] = time.perf_counter()

    def lap(label):
        timings[label] = time.perf_counter() - timings[label]

    clock("orbit_track")
    # The reference curve is a DEFINITION shared by all routes, so relaxed
    # integration controls only move the curve by nanometers (absorbed by
    # the H1 gate).  With the atom-parallel single-point evaluator the
    # EXACT source costs the same as the tree here, so the tracker keeps
    # full source purity.
    orbit = track_reference_orbit(
        b_point, float(options.magnetic_rigidity), station_count=65,
        entrance_x_m=-0.040, exit_x_m=0.040,
        entrance_offset_m=float(offset_m),
        relative_tolerance=1.0e-8, maximum_step_m=1.0e-2)
    lap("orbit_track")
    s_total = float(orbit.arc_length_stations[-1])
    seg_mids = 0.5 * (orbit.arc_length_stations[:-1]
                      + orbit.arc_length_stations[1:])

    def htilde_of_s(s_value):
        return -float(np.interp(s_value, seg_mids, orbit.signed_curvature))

    clock("monitor_chain")
    monitor_s = np.linspace(0.0, s_total, 401)
    positions = orbit.position_at(monitor_s)
    verticals = np.asarray(
        [frame_axes(orbit, sv)[1] for sv in monitor_s])
    by_orbit = np.einsum("ij,ij->i", b_batch(positions), verticals)
    monitor = np.abs(np.gradient(by_orbit, monitor_s))
    breaks = graded_breaks(monitor_s, monitor, int(options.elements),
                           strength=float(options.grade_fringe))
    chain = CanonicalHCurlChain(
        breaks, float(options.half_width), float(options.half_height),
        order_x=int(options.order_x), order_s=int(options.order_s),
        curvature_per_m=htilde_of_s)
    lap("monitor_chain")

    clock("fit_cloud")
    fit_cloud = sample_frame_cloud(
        orbit, b_batch, (0.0, s_total), float(options.half_width),
        float(options.half_height), rng, int(options.fit_stations),
        int(options.fit_points_per_station), breaks=breaks)
    lap("fit_cloud")
    clock("fit_ls")
    fit = chain.fit_frame_samples(*fit_cloud)
    lap("fit_ls")
    clock("audit")
    audit_cloud = sample_frame_cloud(
        orbit, b_batch, (0.0, s_total), float(options.half_width),
        float(options.half_height), rng, int(options.fit_stations) // 2,
        int(options.fit_points_per_station), breaks=breaks)
    audit_b = chain.magnetic_flux_density_frame(*audit_cloud[:3])
    audit_max = float(np.max(np.linalg.norm(
        audit_b - np.column_stack(audit_cloud[3:]), axis=1)))
    lap("audit")

    clock("lie")
    Ay, As, lengths, curvatures = chain.lie_element_spoly_arrays(
        degree=int(options.lie_degree))
    lie = _fourth_order_lie_map_from_vector_potential_polynomials(
        Ay, As, lengths, float(options.magnetic_rigidity),
        reference_curvature_per_m=curvatures,
        longitudinal_component="covariant",
        reference_orbit_tolerance=float(options.reference_orbit_tolerance),
        parameter_jacobians=False,
    )
    lap("lie")
    linear_max = float(np.max(np.abs(lie.hamiltonian_linear), initial=0.0))

    rigidity = float(options.magnetic_rigidity)
    scale = float(options.track_scale)
    mechanical0 = scale * np.array(
        [1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 0.0, 1.0e-3])

    def chain_a_normalized(x_value, y_value, s_value):
        return chain.vector_potential_frame(
            np.array([x_value]), np.array([y_value]),
            np.array([s_value]))[0] / rigidity

    a_entry = chain_a_normalized(mechanical0[0], mechanical0[2], 0.0)
    canonical_entry = mechanical0.copy()
    canonical_entry[1] += a_entry[0]
    canonical_entry[3] += a_entry[1]

    lie_exit = apply_dragt_finn_map(
        lie.transfer.factorization, canonical_entry, generator_substeps=8)
    clock("a_rk")
    a_rk_exit, _ = track_chain_a_rk(
        chain, htilde_of_s, rigidity, canonical_entry, (0.0, s_total),
        float(options.rk_tolerance))
    lap("a_rk")

    def to_mechanical(state):
        a_exit = chain_a_normalized(state[0], state[2], s_total)
        return np.array([state[0], state[1] - a_exit[0],
                         state[2], state[3] - a_exit[1]])

    clock("b_rk")
    b_rk_mech, _ = track_mechanical_b_rk(
        b_point, orbit, rigidity, mechanical0, (0.0, s_total),
        float(options.b_rk_tolerance))
    lap("b_rk")
    lie_vs_a_rk = float(np.max(np.abs(lie_exit - a_rk_exit)))
    a_rk_vs_b_rk = float(np.max(np.abs(to_mechanical(a_rk_exit) - b_rk_mech)))
    lie_vs_b_rk = float(np.max(np.abs(to_mechanical(lie_exit) - b_rk_mech)))

    h_mid = htilde_of_s(0.5 * s_total)
    return {
        "offset_mm": float(offset_m) * 1e3,
        "arc_length_m": s_total,
        "htilde_mid_per_m": h_mid,
        "fit_relative_residual": fit.relative_residual,
        "audit_honesty_ratio": audit_max / fit.maximum_residual_t,
        "hamiltonian_linear_max": linear_max,
        "lie_vs_a_rk": lie_vs_a_rk,
        "a_rk_vs_b_rk": a_rk_vs_b_rk,
        "lie_vs_b_rk": lie_vs_b_rk,
        "stage_timings_s": {k: round(v, 2) for k, v in timings.items()},
    }


def main(argv=None):
    options = parser().parse_args(argv)
    rad.UtiDelAll()
    ng.SetNumThreads(int(options.threads))
    started = time.perf_counter()
    with ng.TaskManager():
        iron_mesh = build_iron(float(options.iron_maxh))
        coil = build_coil(float(options.current))
        iron = vim.MeshSoftIron(iron_mesh, bh_table=load_bh_table(), order=1)
        model = rad.ObjCnt([iron, coil])
        solution = rad.Solve(model)
    print(f"solve {time.perf_counter() - started:.1f} s "
          f"({iron_mesh.ne} elements); vertical region |y| <= "
          f"{1e3 * float(options.half_height):.1f} mm "
          f"({100 * float(options.half_height) / 0.005:.0f}% of the half-gap)")
    b_point, b_batch_exact = make_symmetric_b_field(solution, coil)
    tree_options = {"theta": float(options.tree_theta), "leaf_size": 16,
                    "tree_min_sources": 64}
    b_batch_tree = make_symmetric_b_batch_tree(solution, coil, tree_options)
    rng = np.random.default_rng(20260818)

    # In-run gate: the relaxed tree must match the exact evaluator on a
    # probe cloud before any fit consumes it.
    probe = np.column_stack((rng.uniform(-0.035, 0.035, 500),
                             rng.uniform(-0.012, 0.012, 500),
                             rng.uniform(-float(options.half_height),
                                         float(options.half_height), 500)))
    exact = b_batch_exact(probe)
    tree = b_batch_tree(probe)
    tree_error = float(np.max(np.linalg.norm(tree - exact, axis=1))
                       / np.max(np.linalg.norm(exact, axis=1)))
    print(f"tree gate (theta {float(options.tree_theta):g}): "
          f"max rel {tree_error:.2e} vs exact on 500 probes")
    if tree_error > float(options.tree_gate):
        raise RuntimeError(
            f"relaxed tree failed the accuracy gate: {tree_error:.2e} > "
            f"{float(options.tree_gate):.2e}")

    records = []
    for offset_mm in options.orbit_offsets_mm:
        orbit_started = time.perf_counter()
        record = run_orbit(float(offset_mm) * 1e-3, options, b_point,
                           b_batch_tree, rng)
        record["runtime_s"] = time.perf_counter() - orbit_started
        records.append(record)
        print(f"orbit {record['offset_mm']:+5.1f} mm "
              f"(L={record['arc_length_m']*1e3:6.2f} mm, "
              f"h={record['htilde_mid_per_m']:+.4f}): "
              f"fit rel {record['fit_relative_residual']:.1e} "
              f"honesty {record['audit_honesty_ratio']:.2f} "
              f"H1 {record['hamiltonian_linear_max']:.1e}  "
              f"LIE-ARK {record['lie_vs_a_rk']:.2e}  "
              f"ARK-BRK {record['a_rk_vs_b_rk']:.2e}  "
              f"LIE-BRK {record['lie_vs_b_rk']:.2e}  "
              f"[{record['runtime_s']:.0f} s: "
              + " ".join(f"{k}={v:.1f}"
                         for k, v in record["stage_timings_s"].items())
              + "]")

    result = {
        "half_height_m": float(options.half_height),
        "orbits": records,
    }
    if options.output:
        import platform
        from datetime import datetime
        result["generated_at_utc"] = datetime.now(
            UTC).isoformat()
        result["hostname"] = platform.node()
        result["arguments"] = dict(vars(options))
        Path(options.output).write_text(
            json.dumps(result, indent=2), encoding="utf-8")
        print(f"result written to {options.output}")
    rad.UtiDelAll()
    return result


if __name__ == "__main__":
    main()
