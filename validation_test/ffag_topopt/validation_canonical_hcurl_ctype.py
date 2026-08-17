"""CanonicalHCurl chain on the real C-type: fit certificates + LIE vs A-RK.

Stage-3 first verification of the CanonicalHCurl production space
(:mod:`radia.beam_canonical_hcurl`) against the validation C-type magnet:

  1. solve HDiv-MMM, track the design orbit;
  2. fit a CanonicalHCurl chain by the full-volume 3-component frame B fit
     (correct ``frame_at`` order: horizontal, vertical, tangent);
  3. certificates: honest fit residual vs an independent audit cloud,
     hard interface continuity, vacuum defect spectrum, dimension law;
  4. build the fourth-order Lie map from the chain's per-segment covariant
     transverse coefficient arrays (midpoint s-staging) and compare against
     the NATIVE canonical A-map RK (DOP853 on
     ``canonical_vector_potential_hamiltonian_rhs`` with the chain's exact
     polynomial A and gradients) -- the pair isolates the Lie truncation.

The independent HDiv B-map RK cross-route and the exact-source (ungauged)
A-RK reuse the established three-route harness and are wired in a follow-up
run; the sign convention ``htilde = -signed_curvature`` follows the
EarlyTimes metric contract and is pinned against B-RK there.

Usage (LAB smoke; heavy sweeps go to mdx/hibino):
  python validation_canonical_hcurl_ctype.py --iron-maxh 0.02
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import ngsolve as ng  # noqa: E402
import radia as rad  # noqa: E402
from radia import vim  # noqa: E402
from radia.beam_canonical_hcurl import CanonicalHCurlChain  # noqa: E402
from radia.accelerator_lie_topopt import (  # noqa: E402
    _fourth_order_lie_map_from_vector_potential_polynomials,
    apply_dragt_finn_map,
    canonical_vector_potential_hamiltonian_rhs,
)

from validation_earlytimes_ctype_ab import (  # noqa: E402
    build_iron, build_coil, load_bh_table, make_symmetric_b_field,
    track_reference_orbit,
)


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--iron-maxh", type=float, default=0.02)
    result.add_argument("--threads", type=int, default=8)
    result.add_argument("--current", type=float, default=3000.0)
    result.add_argument("--magnetic-rigidity", type=float, default=3.0)
    result.add_argument("--half-width", type=float, default=0.010)
    result.add_argument("--half-height", type=float, default=0.0035)
    result.add_argument("--order-x", type=int, default=8)
    result.add_argument("--order-s", type=int, default=2)
    result.add_argument("--elements", type=int, default=16)
    result.add_argument("--fit-stations", type=int, default=48)
    result.add_argument("--fit-points-per-station", type=int, default=120)
    result.add_argument("--lie-segments", type=int, default=64)
    result.add_argument("--lie-degree", type=int, default=5)
    result.add_argument("--reference-orbit-tolerance", type=float,
                        default=2.0e-3)
    result.add_argument("--htilde-sign", type=float, default=-1.0,
                        help="htilde = sign * signed_curvature (metric "
                             "g = 1 + htilde*x); the Hamiltonian-linear "
                             "gate pins the correct convention")
    result.add_argument("--track-scales", type=float, nargs="+",
                        default=[0.25, 0.5, 1.0])
    result.add_argument("--rk-tolerance", type=float, default=1.0e-12)
    result.add_argument("--output", type=str, default=None)
    return result


def sample_frame_cloud(orbit, b_batch, s_range, half_width, half_height,
                       rng, n_s, n_xy):
    s_values = rng.uniform(s_range[0], s_range[1], n_s)
    xs, ys, ss, bxl, byl, bsl = [], [], [], [], [], []
    for s_value in s_values:
        sv = np.array([s_value])
        centre = orbit.position_at(sv)[0]
        horizontal, vertical, tangent = (np.asarray(w[0], dtype=float)
                                         for w in orbit.frame_at(sv))
        x = rng.uniform(-half_width, half_width, n_xy)
        y = rng.uniform(-half_height, half_height, n_xy)
        points = centre[None, :] + x[:, None] * horizontal[None, :] \
            + y[:, None] * vertical[None, :]
        b = b_batch(points)
        xs.append(x)
        ys.append(y)
        ss.append(np.full(n_xy, s_value))
        bxl.append(b @ horizontal)
        byl.append(b @ vertical)
        bsl.append(b @ tangent)
    return (np.concatenate(xs), np.concatenate(ys), np.concatenate(ss),
            np.concatenate(bxl), np.concatenate(byl), np.concatenate(bsl))


def track_chain_a_rk(chain, curvature_of_s, rigidity, state, s_span,
                     rk_tolerance):
    """Native canonical A-map RK on the chain's exact polynomials."""

    def rhs(s_value, z):
        a, gradient = chain.vector_potential_and_gradient_frame(
            np.array([float(z[0])]), np.array([float(z[2])]),
            np.array([s_value]))
        return canonical_vector_potential_hamiltonian_rhs(
            z, a[0] / rigidity, gradient[0] / rigidity,
            reference_curvature_per_m=float(curvature_of_s(s_value)),
        )

    solution = solve_ivp(
        rhs, s_span, np.asarray(state, dtype=float), method="DOP853",
        rtol=rk_tolerance, atol=rk_tolerance, dense_output=False,
    )
    if not solution.success:
        raise RuntimeError(f"A-map RK failed: {solution.message}")
    return solution.y[:, -1], int(solution.nfev)


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
    solve_wall = time.perf_counter() - started
    print(f"solve {solve_wall:.1f} s ({iron_mesh.ne} elements)")

    b_point, b_batch = make_symmetric_b_field(solution, coil)
    orbit = track_reference_orbit(
        b_point, float(options.magnetic_rigidity), station_count=65,
        entrance_x_m=-0.040, exit_x_m=0.040)
    s_total = float(orbit.arc_length_stations[-1])
    segment_mids = 0.5 * (orbit.arc_length_stations[:-1]
                          + orbit.arc_length_stations[1:])

    def htilde_of_s(s_value):
        # EarlyTimes metric contract: g = 1 + htilde*x.  The sign relative
        # to orbit.signed_curvature is pinned empirically by the
        # Hamiltonian-linear gate (and against B-RK in the follow-up).
        return float(options.htilde_sign) * float(
            np.interp(s_value, segment_mids, orbit.signed_curvature))

    chain = CanonicalHCurlChain(
        np.linspace(0.0, s_total, int(options.elements) + 1),
        float(options.half_width), float(options.half_height),
        order_x=int(options.order_x), order_s=int(options.order_s),
        curvature_per_m=htilde_of_s)
    print(f"chain: {chain.element_count} elements, "
          f"dim/element={chain.elements[0].dimension}, "
          f"chain dim={chain.chain_dimension} "
          f"(spline law {int(options.order_x)}*"
          f"(E+{int(options.order_s)})="
          f"{int(options.order_x) * (chain.element_count + int(options.order_s))})")
    defect = max(float(np.max(el.vacuum_defects / el.vacuum_defect_scale))
                 for el in chain.elements)
    print(f"vacuum defect (relative, worst element): {defect:.2e}")
    if chain.interface_defects.size:
        print(f"interface defect spectrum (relative): max "
              f"{float(np.max(chain.interface_defects)) / chain.interface_defect_scale:.2e}")

    rng = np.random.default_rng(20260817)
    fit_started = time.perf_counter()
    fit_cloud = sample_frame_cloud(
        orbit, b_batch, (0.0, s_total), float(options.half_width),
        float(options.half_height), rng, int(options.fit_stations),
        int(options.fit_points_per_station))
    fit = chain.fit_frame_samples(*fit_cloud)
    audit_cloud = sample_frame_cloud(
        orbit, b_batch, (0.0, s_total), float(options.half_width),
        float(options.half_height), rng, int(options.fit_stations) // 2,
        int(options.fit_points_per_station))
    audit_b = chain.magnetic_flux_density_frame(*audit_cloud[:3])
    audit_ref = np.column_stack(audit_cloud[3:])
    audit_max = float(np.max(np.linalg.norm(audit_b - audit_ref, axis=1)))
    fit_wall = time.perf_counter() - fit_started
    print(f"fit ({fit.sample_count} samples, {fit_wall:.1f} s): "
          f"max residual {fit.maximum_residual_t:.3e} T "
          f"(rel {fit.relative_residual:.1e})")
    print(f"  independent audit: max |dB| {audit_max:.3e} T "
          f"(honesty ratio {audit_max / fit.maximum_residual_t:.2f})")
    print(f"  interface jumps: a_y {fit.maximum_interface_ay_jump:.2e}, "
          f"b-value {fit.maximum_interface_b_value_jump:.2e}  (hard)")

    # ---- fourth-order Lie map from the chain ---------------------------
    lie_started = time.perf_counter()
    Ay, As, lengths, curvatures = chain.lie_segment_arrays(
        int(options.lie_segments), degree=int(options.lie_degree))
    lie = _fourth_order_lie_map_from_vector_potential_polynomials(
        Ay, As, lengths, float(options.magnetic_rigidity),
        reference_curvature_per_m=curvatures,
        longitudinal_component="covariant",
        reference_orbit_tolerance=float(options.reference_orbit_tolerance),
    )
    lie_wall = time.perf_counter() - lie_started
    linear_max = float(np.max(np.abs(lie.hamiltonian_linear), initial=0.0))
    print(f"LIE map ({int(options.lie_segments)} segments, {lie_wall:.1f} s): "
          f"max |hamiltonian linear| {linear_max:.2e} "
          f"(orbit-consistency gate)")

    # ---- native canonical A-map RK vs LIE ------------------------------
    base = np.array([1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 0.0, 1.0e-3])
    records = []
    for scale in options.track_scales:
        state = float(scale) * base
        lie_state = apply_dragt_finn_map(
            lie.transfer.factorization, state, generator_substeps=8)
        rk_state, nfev = track_chain_a_rk(
            chain, htilde_of_s, float(options.magnetic_rigidity), state,
            (0.0, s_total), float(options.rk_tolerance))
        difference = lie_state - rk_state
        records.append({
            "scale": float(scale),
            "initial_state": state.tolist(),
            "lie_state": np.asarray(lie_state).tolist(),
            "a_rk_state": rk_state.tolist(),
            "lie_minus_a_rk": difference.tolist(),
            "maximum_difference": float(np.max(np.abs(difference))),
            "a_rk_evaluations": nfev,
        })
        print(f"  scale {scale:5.2f}: max|LIE - A-RK| "
              f"{records[-1]['maximum_difference']:.3e}  (nfev {nfev})")
    ratios = [records[i]["maximum_difference"]
              / max(records[i - 1]["maximum_difference"], 1e-300)
              for i in range(1, len(records))]
    if ratios:
        print(f"  amplitude scaling of the truncation error: "
              f"{['%.1f' % value for value in ratios]} "
              f"(degree-5 jet => ~2^5=32 per doubling)")

    import platform
    from datetime import datetime, timezone

    result = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(),
        "arguments": {key: value for key, value in vars(options).items()},
        "iron_maxh": float(options.iron_maxh),
        "elements": chain.element_count,
        "order_x": int(options.order_x),
        "order_s": int(options.order_s),
        "chain_dimension": chain.chain_dimension,
        "vacuum_defect_relative": defect,
        "fit_maximum_residual_t": fit.maximum_residual_t,
        "fit_relative_residual": fit.relative_residual,
        "audit_maximum_residual_t": audit_max,
        "interface_ay_jump": fit.maximum_interface_ay_jump,
        "interface_b_value_jump": fit.maximum_interface_b_value_jump,
        "hamiltonian_linear_max": linear_max,
        "lie_vs_a_rk": records,
        "solve_wall_s": solve_wall,
    }
    if options.output:
        Path(options.output).write_text(
            json.dumps(result, indent=2), encoding="utf-8")
        print(f"result written to {options.output}")
    rad.UtiDelAll()
    return result


if __name__ == "__main__":
    main()
