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

Usage (LAB smoke; heavy sweeps use hibino first or an idle-CI mdx fallback):
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

from datetime import UTC

import ngsolve as ng
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
    canonical_vector_potential_hamiltonian_rhs,
)
from radia.beam_canonical_hcurl import (
    CanonicalHCurlChain,
    graded_breaks,
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
    result.add_argument("--lie-staging",
                        choices=("midpoint", "element-spoly"),
                        default="midpoint",
                        help="midpoint: piecewise-constant segments; "
                             "element-spoly: one Lie segment per chain "
                             "element consuming its zeta-polynomial "
                             "directly (nonautonomous stage-jet RK4)")
    result.add_argument("--reference-orbit-tolerance", type=float,
                        default=2.0e-3)
    result.add_argument("--htilde-sign", type=float, default=-1.0,
                        help="htilde = sign * signed_curvature (metric "
                             "g = 1 + htilde*x); the Hamiltonian-linear "
                             "gate pins the correct convention")
    result.add_argument("--grade-fringe", type=float, default=0.0,
                        help="graded_breaks strength on the |dBy/ds| orbit "
                             "monitor (0 = uniform elements)")
    result.add_argument("--include-cross-routes", action="store_true",
                        help="also run the exact-source canonical A-RK "
                             "(raw analytic A, FD gradients) and the native "
                             "Cartesian mechanical B-RK, compared in frame "
                             "mechanical variables at the exit plane")
    result.add_argument("--cross-scale", type=float, default=0.5)
    result.add_argument("--fd-step", type=float, default=1.0e-6)
    result.add_argument("--track-scales", type=float, nargs="+",
                        default=[0.25, 0.5, 1.0])
    result.add_argument("--rk-tolerance", type=float, default=1.0e-12)
    result.add_argument("--output", type=str, default=None)
    return result


def sample_frame_cloud(orbit, b_batch, s_range, half_width, half_height,
                       rng, n_s, n_xy, breaks=None):
    """Frame-component B cloud; with ``breaks``, stations are allocated per
    element (proportional with a floor) so thin graded elements never
    starve -- uniform-in-s sampling under grading leaves near-empty
    elements whose DOFs the LS drives wild (caught by the fit guard).

    ALL points go through ONE ``b_batch`` call: per-station mini-batches
    paid the Python->C++ marshalling and thread-pool standup ~70 times per
    cloud and dominated the whole pipeline (~500 s/orbit measured); the
    single batched call evaluates the same exact kernels in seconds.
    """
    if breaks is None:
        s_values = rng.uniform(s_range[0], s_range[1], n_s)
    else:
        sizes = np.diff(np.asarray(breaks, dtype=float))
        floor = 3
        allocation = np.maximum(
            floor, np.round(n_s * sizes / sizes.sum()).astype(int))
        s_values = np.concatenate([
            rng.uniform(breaks[e], breaks[e + 1], allocation[e])
            for e in range(len(sizes))
        ])
    xs, ys, ss, points_list, frames = [], [], [], [], []
    for s_value in s_values:
        sv = np.array([s_value])
        centre = orbit.position_at(sv)[0]
        horizontal, vertical, tangent = (np.asarray(w[0], dtype=float)
                                         for w in orbit.frame_at(sv))
        x = rng.uniform(-half_width, half_width, n_xy)
        # The basis parity is exact and the source is symmetrized, so the
        # y<0 half of the cloud carries no independent information; the
        # upper half alone halves the dominant source-evaluation cost.
        y = rng.uniform(0.0, half_height, n_xy)
        points_list.append(
            centre[None, :] + x[:, None] * horizontal[None, :]
            + y[:, None] * vertical[None, :])
        frames.append((horizontal, vertical, tangent))
        xs.append(x)
        ys.append(y)
        ss.append(np.full(n_xy, s_value))
    b = b_batch(np.vstack(points_list))
    bxl, byl, bsl = [], [], []
    for station, (horizontal, vertical, tangent) in enumerate(frames):
        block = b[station * n_xy:(station + 1) * n_xy]
        bxl.append(block @ horizontal)
        byl.append(block @ vertical)
        bsl.append(block @ tangent)
    return (np.concatenate(xs), np.concatenate(ys), np.concatenate(ss),
            np.concatenate(bxl), np.concatenate(byl), np.concatenate(bsl))


def frame_axes(orbit, s_value):
    horizontal, vertical, tangent = (np.asarray(w[0], dtype=float)
                                     for w in orbit.frame_at(
                                         np.array([s_value])))
    return horizontal, vertical, tangent


def make_frame_a(a_of_point, orbit, htilde_of_s):
    """Covariant frame components of a raw global vector potential."""

    def frame_a(x_value, y_value, s_value):
        horizontal, vertical, tangent = frame_axes(orbit, s_value)
        point = (orbit.position_at(np.array([s_value]))[0]
                 + x_value * horizontal + y_value * vertical)
        value = np.asarray(a_of_point(point), dtype=float)
        metric = 1.0 + htilde_of_s(s_value) * x_value
        return np.array([value @ horizontal, value @ vertical,
                         metric * (value @ tangent)])

    return frame_a


def track_exact_source_a_rk(frame_a, htilde_of_s, rigidity, state,
                            s_span, rk_tolerance, fd_step):
    """Canonical A-map RK on the RAW analytic A (no projection, no gauge).

    ``frame_a`` returns the covariant frame components at ``(x, y, s)``;
    transverse gradients use central finite differences (the established
    exact-A route floor).  ``a_x != 0`` here -- the canonical machinery
    accepts it and the mechanical conversion removes the gauge.
    """
    evaluations = 0

    def rhs(s_value, z):
        nonlocal evaluations
        x_value, y_value = float(z[0]), float(z[2])
        centre = frame_a(x_value, y_value, s_value)
        gradient = np.empty((3, 2))
        gradient[:, 0] = (frame_a(x_value + fd_step, y_value, s_value)
                          - frame_a(x_value - fd_step, y_value, s_value)) \
            / (2.0 * fd_step)
        gradient[:, 1] = (frame_a(x_value, y_value + fd_step, s_value)
                          - frame_a(x_value, y_value - fd_step, s_value)) \
            / (2.0 * fd_step)
        evaluations += 5
        return canonical_vector_potential_hamiltonian_rhs(
            z, centre / rigidity, gradient / rigidity,
            reference_curvature_per_m=float(htilde_of_s(s_value)))

    solution = solve_ivp(
        rhs, s_span, np.asarray(state, dtype=float), method="DOP853",
        rtol=rk_tolerance, atol=rk_tolerance)
    if not solution.success:
        raise RuntimeError(f"exact-source A-RK failed: {solution.message}")
    return solution.y[:, -1], evaluations


def track_mechanical_b_rk(b_point, orbit, rigidity, mechanical, s_span,
                          rk_tolerance):
    """Native Cartesian mechanical RK on the symmetrized HDiv-MMM B.

    Integrates ``dr/dtau = u``, ``du/dtau = u x B / (B rho (1+delta))``
    (path-length parameterization, unit direction ``u``) from the entrance
    frame mechanical state and stops at the exit plane through
    ``orbit(s_end)`` perpendicular to the exit tangent.  Fully independent
    of every A representation.
    """
    x0, pxm0, y0, pym0, _, delta = (float(v) for v in mechanical)
    h_entry, v_entry, t_entry = frame_axes(orbit, s_span[0])
    r_entry = (orbit.position_at(np.array([s_span[0]]))[0]
               + x0 * h_entry + y0 * v_entry)
    w0 = np.sqrt((1.0 + delta) ** 2 - pxm0**2 - pym0**2)
    u0 = (pxm0 * h_entry + pym0 * v_entry + w0 * t_entry) / (1.0 + delta)
    h_exit, v_exit, t_exit = frame_axes(orbit, s_span[1])
    r_exit_ref = orbit.position_at(np.array([s_span[1]]))[0]
    inverse_bend = 1.0 / (rigidity * (1.0 + delta))

    def rhs(_tau, state):
        r = state[:3]
        u = state[3:]
        b = np.asarray(b_point(r), dtype=float)
        return np.concatenate((u, np.cross(u, b) * inverse_bend))

    def exit_plane(_tau, state):
        return float((state[:3] - r_exit_ref) @ t_exit)

    exit_plane.terminal = True
    exit_plane.direction = 1.0
    length = s_span[1] - s_span[0]
    solution = solve_ivp(
        rhs, (0.0, 1.5 * length), np.concatenate((r_entry, u0)),
        method="DOP853", rtol=rk_tolerance, atol=rk_tolerance,
        events=exit_plane)
    if not solution.success or not len(solution.t_events[0]):
        raise RuntimeError("mechanical B-RK did not reach the exit plane")
    final = solution.y_events[0][0]
    r_end, u_end = final[:3], final[3:]
    return np.array([
        (r_end - r_exit_ref) @ h_exit,
        (1.0 + delta) * (u_end @ h_exit),
        (r_end - r_exit_ref) @ v_exit,
        (1.0 + delta) * (u_end @ v_exit),
    ]), int(solution.nfev)


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

    if float(options.grade_fringe) > 0.0:
        monitor_s = np.linspace(0.0, s_total, 401)
        by_orbit = np.array([
            float(np.asarray(b_point(orbit.position_at(np.array([sv]))[0]))
                  @ frame_axes(orbit, sv)[1])
            for sv in monitor_s
        ])
        monitor = np.abs(np.gradient(by_orbit, monitor_s))
        breaks = graded_breaks(monitor_s, monitor, int(options.elements),
                               strength=float(options.grade_fringe))
        sizes = np.diff(breaks)
        print(f"graded breaks (strength {float(options.grade_fringe):g}): "
              f"element sizes {1e3*sizes.min():.2f}..{1e3*sizes.max():.2f} mm")
    else:
        breaks = np.linspace(0.0, s_total, int(options.elements) + 1)
    chain = CanonicalHCurlChain(
        breaks,
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
        int(options.fit_points_per_station), breaks=breaks)
    fit = chain.fit_frame_samples(*fit_cloud)
    audit_cloud = sample_frame_cloud(
        orbit, b_batch, (0.0, s_total), float(options.half_width),
        float(options.half_height), rng, int(options.fit_stations) // 2,
        int(options.fit_points_per_station), breaks=breaks)
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
    if options.lie_staging == "element-spoly":
        Ay, As, lengths, curvatures = chain.lie_element_spoly_arrays(
            degree=int(options.lie_degree))
        staging_note = f"{chain.element_count} element-spoly segments"
    else:
        Ay, As, lengths, curvatures = chain.lie_segment_arrays(
            int(options.lie_segments), degree=int(options.lie_degree))
        staging_note = f"{int(options.lie_segments)} midpoint segments"
    lie = _fourth_order_lie_map_from_vector_potential_polynomials(
        Ay, As, lengths, float(options.magnetic_rigidity),
        reference_curvature_per_m=curvatures,
        longitudinal_component="covariant",
        reference_orbit_tolerance=float(options.reference_orbit_tolerance),
        # Verification consumes only the map; the topopt adjoint over
        # 40*n_segment parameters is the measured 264x runtime hog.
        parameter_jacobians=False,
    )
    lie_wall = time.perf_counter() - lie_started
    linear_max = float(np.max(np.abs(lie.hamiltonian_linear), initial=0.0))
    print(f"LIE map ({staging_note}, {lie_wall:.1f} s): "
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
              f"{[f'{value:.1f}' for value in ratios]} "
              f"(degree-5 jet => ~2^5=32 per doubling)")

    import platform
    from datetime import datetime

    cross_record = None
    if options.include_cross_routes:
        cross_started = time.perf_counter()
        from radia.beam import build_curvilinear_beam_mesh
        margin = 0.0005
        extended_orbit = track_reference_orbit(
            b_point, float(options.magnetic_rigidity), station_count=65,
            entrance_x_m=-0.040 - margin, exit_x_m=0.040 + margin)
        with ng.TaskManager():
            tube = build_curvilinear_beam_mesh(
                extended_orbit, half_width_m=float(options.half_width),
                half_height_m=float(options.half_height),
                maxh_m=0.02, vertical_layers=3, curve_order=2)
            demag_a = vim.VectorPotentialCoefficientFromSolution(
                solution, construction="exact",
                reflection_normal=orbit.bend_axis)
            coil_a = rad.RadiaField(coil, "a")
            reflected_coil_a = rad.RadiaField(
                coil, "a", origin=[0.0, 0.0, 0.0],
                u_axis=[1.0, 0.0, 0.0], v_axis=[0.0, 1.0, 0.0],
                w_axis=[0.0, 0.0, -1.0])
            current_a = demag_a + 0.5 * (coil_a + reflected_coil_a)

        def a_of_point(point):
            return np.asarray(current_a(tube.mesh(*point)), dtype=float)

        rigidity = float(options.magnetic_rigidity)
        scale = float(options.cross_scale)
        mechanical0 = scale * np.array(
            [1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 0.0, 1.0e-3])

        def chain_a_normalized(x_value, y_value, s_value):
            value = chain.vector_potential_frame(
                np.array([x_value]), np.array([y_value]),
                np.array([s_value]))[0]
            return value / rigidity

        def exact_a_normalized(frame_a, x_value, y_value, s_value):
            return frame_a(x_value, y_value, s_value) / rigidity

        # Chain-canonical entrance state (chain's own gauge).
        a_entry_chain = chain_a_normalized(mechanical0[0], mechanical0[2], 0.0)
        chain_entry = mechanical0.copy()
        chain_entry[1] += a_entry_chain[0]
        chain_entry[3] += a_entry_chain[1]

        lie_exit = apply_dragt_finn_map(
            lie.transfer.factorization, chain_entry, generator_substeps=8)
        chain_exit, _chain_nfev = track_chain_a_rk(
            chain, htilde_of_s, rigidity, chain_entry, (0.0, s_total),
            float(options.rk_tolerance))

        def chain_mechanical(state):
            a_exit = chain_a_normalized(state[0], state[2], s_total)
            return np.array([state[0], state[1] - a_exit[0],
                             state[2], state[3] - a_exit[1]])

        exact_frame_a = make_frame_a(a_of_point, orbit, htilde_of_s)
        # Entrance conversion for the exact-source route (its own raw A).
        a_entry_exact = exact_a_normalized(
            exact_frame_a, mechanical0[0], mechanical0[2], 0.0)
        exact_entry = mechanical0.copy()
        exact_entry[1] += a_entry_exact[0]
        exact_entry[3] += a_entry_exact[1]
        exact_state, exact_evaluations = track_exact_source_a_rk(
            exact_frame_a, htilde_of_s, rigidity,
            exact_entry, (0.0, s_total), 1.0e-11, float(options.fd_step))
        a_exit_exact = exact_a_normalized(
            exact_frame_a, exact_state[0], exact_state[2], s_total)
        exact_mech = np.array([
            exact_state[0], exact_state[1] - a_exit_exact[0],
            exact_state[2], exact_state[3] - a_exit_exact[1]])

        b_mech, b_nfev = track_mechanical_b_rk(
            b_point, orbit, rigidity, mechanical0, (0.0, s_total),
            float(options.rk_tolerance))

        lie_mech = chain_mechanical(lie_exit)
        chain_mech = chain_mechanical(chain_exit)
        pairs = {
            "lie_vs_chain_a_rk": lie_mech - chain_mech,
            "chain_a_rk_vs_exact_a_rk": chain_mech - exact_mech,
            "chain_a_rk_vs_b_rk": chain_mech - b_mech,
            "exact_a_rk_vs_b_rk": exact_mech - b_mech,
        }
        print(f"\ncross-route mechanical comparison (scale {scale:g}, "
              f"[x, px_mech, y, py_mech], {time.perf_counter() - cross_started:.0f} s):")
        for name, difference in pairs.items():
            print(f"  {name:26s}: max {float(np.max(np.abs(difference))):.3e}"
                  f"  {np.array2string(difference, precision=2)}")
        print(f"  [exact-source A evals {exact_evaluations}, "
              f"B-RK nfev {b_nfev}]")
        cross_record = {
            "scale": scale,
            "mechanical_initial": mechanical0.tolist(),
            "lie_mechanical": lie_mech.tolist(),
            "chain_a_rk_mechanical": chain_mech.tolist(),
            "exact_a_rk_mechanical": exact_mech.tolist(),
            "b_rk_mechanical": b_mech.tolist(),
            "pairwise_max": {name: float(np.max(np.abs(value)))
                             for name, value in pairs.items()},
        }

    result = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
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
        "cross_routes": cross_record,
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
