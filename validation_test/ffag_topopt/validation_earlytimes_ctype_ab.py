"""C-type-magnet convergence study for HCurl-A, Lie, and HDiv-B tracking.

It solves one physical C-type magnet and builds the constrained
``HCurl(order=5)`` A-map.  The optional A-only lane compares the fourth-order
Dragt--Finn Lie factors and direct ``R/T/U/V`` polynomial against the exact
unexpanded canonical A-RK Hamiltonian.  The independent HDiv-MMM B
``HDiv(order=4)`` GridFunction route remains available as a separate check.
The original B CoefficientFunction is retained only for the design-orbit and
source-to-HDiv projection diagnostics.

Heavy sweeps belong on a compute host.  The caller supplies an output JSON path;
the file records the host, runtime, exact settings, field/gauge diagnostics, and
all convergence observables needed to compare successive runs.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import time
import traceback
from pathlib import Path

import ngsolve as ng
import numpy as np
from netgen.occ import Face, OCCGeometry, Pnt, Segment, Vec, Wire
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline

import radia as rad
from radia import vim
from radia.accelerator_lie_topopt import (
    HCurlAMapBCoefficientRKBatchComparison,
    compare_earlytimes_b_representations,
    compare_hcurl_a_map_to_hdiv_b_map_rk,
    compare_hcurl_a_map_to_hdiv_b_map_rk_batch,
    compare_hcurl_lie_map_to_a_rk,
    fourth_order_lie_map_from_hcurl_transverse,
    track_b_coefficient_cartesian_s,
    track_hcurl_vector_potential_canonical_s,
    track_hdiv_b_map_cartesian_s,
)
from radia.accelerator_magnet_topopt import PlanarDesignOrbit
from radia.beam import (
    build_curvilinear_beam_mesh,
    project_design_orbit_gauge,
    sample_transverse_vector_potential,
)

MU0 = 4.0e-7 * math.pi
COMPARABLE_INDICES = (0, 1, 2, 3, 5)
COMPONENT_NAMES = ("x", "px", "y", "py", "ell", "delta")


def build_iron(
    maxh_m: float,
    *,
    pole_half_gap_m: float = 0.005,
    pole_half_width_m: float = 0.020,
):
    """Build the intentionally small C-yoke used for the convergence study."""
    pole_half_gap = float(pole_half_gap_m)
    pole_half_width = float(pole_half_width_m)
    if not 0.0 < pole_half_gap < 0.055:
        raise ValueError("pole_half_gap_m must lie between 0 and 0.055 m")
    if not 0.0 < pole_half_width < 0.100:
        raise ValueError("pole_half_width_m must lie between 0 and 0.100 m")
    cross_section = [
        (-pole_half_width, pole_half_gap),
        (-pole_half_width, 0.105),
        (0.1625, 0.105),
        (0.1625, -0.105),
        (-pole_half_width, -0.105),
        (-pole_half_width, -pole_half_gap),
        (pole_half_width, -pole_half_gap),
        (pole_half_width, -0.055),
        (0.100, -0.055),
        (0.100, 0.055),
        (pole_half_width, 0.055),
        (pole_half_width, pole_half_gap),
    ]
    segments = []
    for index, (y0, z0) in enumerate(cross_section):
        y1, z1 = cross_section[(index + 1) % len(cross_section)]
        segments.append(Segment(Pnt(-0.025, y0, z0), Pnt(-0.025, y1, z1)))
    solid = Face(Wire(segments)).Extrude(Vec(0.05, 0.0, 0.0)).mat("iron")
    return ng.Mesh(OCCGeometry(solid).GenerateMesh(maxh=float(maxh_m)))


def build_coil(current_ampere_turns: float):
    half_x = 0.0475
    points = [
        [-half_x, 0.085, 0.0],
        [half_x, 0.085, 0.0],
        [half_x, 0.1775, 0.0],
        [-half_x, 0.1775, 0.0],
        [-half_x, 0.085, 0.0],
    ]
    return rad.ObjFlmCur(points, float(current_ampere_turns))


def load_bh_table():
    path = (
        Path(rad.__file__).resolve().parent
        / "panels"
        / "samples"
        / "em_sample_bh.txt"
    )
    values = np.loadtxt(path, dtype=float)
    return values[:, :2].tolist()


def make_symmetric_b_field(solution, coil):
    """Return independent median-plane-symmetrized HDiv-MMM B evaluators."""

    def batch(points):
        values = np.ascontiguousarray(np.asarray(points, dtype=float).reshape(-1, 3))
        reflected_points = values.copy()
        reflected_points[:, 2] *= -1.0

        def raw(sample_points):
            demag = vim.FieldFromSolution(solution, sample_points)
            coil_b = np.asarray(rad.Fld(coil, "b", sample_points))
            return MU0 * demag + coil_b

        direct = raw(values)
        reflected = raw(reflected_points)
        reflected_axial = reflected * np.array([-1.0, -1.0, 1.0])
        return 0.5 * (direct + reflected_axial)

    def point(position):
        return batch(np.asarray(position, dtype=float).reshape(1, 3))[0]

    return point, batch


def reflection_axes(normal):
    """Return the three column axes of reflection through the origin plane."""
    unit = np.asarray(normal, dtype=float)
    unit /= np.linalg.norm(unit)
    reflection = np.eye(3) - 2.0 * np.outer(unit, unit)
    return tuple(reflection[:, index].tolist() for index in range(3))


def make_midplane_bz_field(
    b_batch,
    *,
    x_range=(-0.045, 0.045),
    y_range=(-0.014, 0.014),
    x_count=101,
    y_count=31,
):
    """Cubic-spline ``B_z(x, y)`` on the ``z=0`` median plane from ONE batch.

    The median-plane orbit ODE of :func:`track_reference_orbit` consumes
    only ``B_z`` on ``z=0``, so a per-sweep gridded spline replaces the
    ~700 sequential single-point kernel evaluations each orbit costs.  This
    interpolates the MIDPLANE TRACE for the reference-CURVE definition only
    -- it is not a 3D field reconstruction (the map fit still samples the
    full volume) -- and callers must gate the returned field against the
    exact source on corridor probes before consuming it.
    """
    xs = np.linspace(float(x_range[0]), float(x_range[1]), int(x_count))
    ys = np.linspace(float(y_range[0]), float(y_range[1]), int(y_count))
    grid_x, grid_y = np.meshgrid(xs, ys, indexing="ij")
    points = np.column_stack(
        (grid_x.reshape(-1), grid_y.reshape(-1), np.zeros(grid_x.size)))
    bz = np.asarray(b_batch(points))[:, 2].reshape(xs.size, ys.size)
    spline = RectBivariateSpline(xs, ys, bz)

    def field(point):
        return np.array([
            0.0,
            0.0,
            float(spline(float(point[0]), float(point[1]))[0, 0]),
        ])

    return field


def track_reference_orbit(
    magnetic_flux_density,
    magnetic_rigidity,
    *,
    station_count,
    entrance_x_m=-0.040,
    exit_x_m=0.040,
    entrance_offset_m=0.0,
    relative_tolerance=2.0e-11,
    maximum_step_m=1.0e-3,
):
    """Track the median-plane design orbit with the independent B source.

    ``entrance_offset_m`` displaces the entrance point transversely (global
    y) with the same +x heading: a family of parallel-entry design orbits
    sampling different horizontal slices of the gap field.  The orbit is a
    REFERENCE-CURVE definition, so callers may relax the integration
    controls (and even supply an accuracy-gated tree source): every route
    shares the same curve and the Hamiltonian-linear gate absorbs the
    consistency question.
    """
    start = np.array([float(entrance_x_m), float(entrance_offset_m), 0.0])
    direction = np.array([1.0, 0.0, 0.0])

    def rhs(_path, state):
        tangent = state[3:]
        field = magnetic_flux_density(state[:3])
        median_plane_field = np.array([0.0, 0.0, field[2]])
        return np.r_[
            tangent,
            np.cross(tangent, median_plane_field) / float(magnetic_rigidity),
        ]

    def exit_plane(_path, state):
        return state[0] - float(exit_x_m)

    exit_plane.terminal = True
    exit_plane.direction = 1.0
    solved = solve_ivp(
        rhs,
        (0.0, 0.14),
        np.r_[start, direction],
        method="DOP853",
        rtol=float(relative_tolerance),
        atol=float(relative_tolerance) * 1.0e-2,
        max_step=float(maximum_step_m),
        events=exit_plane,
        dense_output=True,
    )
    if not solved.success or len(solved.t_events[0]) != 1:
        raise RuntimeError("reference particle did not cross the C-magnet exit plane")
    length = float(solved.t_events[0][0])
    stations = np.linspace(0.0, length, int(station_count))
    states = solved.sol(stations).T
    states[:, 2] = 0.0
    states[:, 5] = 0.0
    tangents = states[:, 3:]
    tangents /= np.linalg.norm(tangents, axis=1)[:, None]
    geometric_orbit = PlanarDesignOrbit(
        positions=states[:, :3],
        tangents=tangents,
        magnetic_rigidity=float(magnetic_rigidity),
        bend_axis=np.array([0.0, 0.0, 1.0]),
        path_length_stations=stations,
    )
    midpoint_s = 0.5 * (stations[:-1] + stations[1:])
    midpoint_positions = geometric_orbit.position_at(midpoint_s)
    _, midpoint_vertical, _ = geometric_orbit.frame_at(midpoint_s)
    midpoint_b = np.asarray(
        [magnetic_flux_density(point) for point in midpoint_positions]
    )
    collocated_curvature = -np.einsum(
        "ij,ij->i", midpoint_b, midpoint_vertical
    ) / float(magnetic_rigidity)
    return PlanarDesignOrbit(
        positions=states[:, :3],
        tangents=tangents,
        magnetic_rigidity=float(magnetic_rigidity),
        bend_axis=np.array([0.0, 0.0, 1.0]),
        path_length_stations=stations,
        signed_curvature_per_m=collocated_curvature,
    )


def trim_orbit(orbit, margin_stations):
    """Drop ``margin_stations`` stations from each orbit end.

    The loft-chain field maps are polluted in their entrance/exit boundary
    elements: the fitted ``d(a_s)/dx`` there reached 3.3e-4 (normalized) on
    the refined mesh while the interior stayed at the few-1e-6 level, because
    a boundary-face section samples the least-controlled one-sided derivative
    content of the discrete fields.  The standard field-map practice applies:
    the certified map region must lie strictly inside the map support.  The
    tube is therefore built from an EXTENDED tracked orbit and the maps and
    tracks consume this trimmed interior orbit.
    """
    margin = int(margin_stations)
    if margin == 0:
        return orbit
    if margin < 0 or 2 * margin + 2 > len(orbit.positions):
        raise ValueError("orbit_margin_stations leaves too few stations")
    return PlanarDesignOrbit(
        positions=orbit.positions[margin:-margin],
        tangents=orbit.tangents[margin:-margin],
        magnetic_rigidity=float(orbit.magnetic_rigidity),
        bend_axis=orbit.bend_axis,
        path_length_stations=orbit.arc_length_stations[margin:-margin],
        signed_curvature_per_m=orbit.signed_curvature[margin:-margin],
    )


def project_axial_gauge(vector_potential, mesh, orbit, order, *,
                        regularization=1.0e-6):
    """Project the C-type field toward ``A_x=0`` before HCurl interpolation."""
    scalar_space = ng.H1(mesh, order=int(order) + 1)
    trial, test = scalar_space.TnT()
    horizontal_samples, _, _ = orbit.frame_at(orbit.arc_length_stations)
    abscissa = orbit.positions[:, 0]
    degree = min(5, len(abscissa) - 1)
    fitted_components = []
    for component in range(3):
        coefficients = np.polynomial.polynomial.polyfit(
            abscissa, horizontal_samples[:, component], degree
        )
        value = ng.CoefficientFunction(0.0)
        for power, coefficient in enumerate(coefficients):
            value = value + float(coefficient) * ng.x**power
        fitted_components.append(value)
    horizontal = ng.CoefficientFunction(tuple(fitted_components))
    trial_axial = ng.InnerProduct(ng.grad(trial), horizontal)
    test_axial = ng.InnerProduct(ng.grad(test), horizontal)
    # The (y,s) directions of chi are controlled ONLY by this term, so it sets
    # the null-space noise floor: assembly roundoff is amplified by
    # eps/regularization.  1e-12 let solver noise populate y-odd modes at the
    # 1e-4 relative level and broke the median-plane parity of A_s by 2.4e-6
    # T*m (measured stage by stage; the source CF is parity-exact and the
    # HCurl projection preserves parity identically).  1e-6 keeps that noise
    # at ~1e-10 relative while biasing the A_x cancellation only at the
    # 1e-6-relative level.
    regularization = float(regularization)
    aperture_scale = 0.01
    form = ng.BilinearForm(scalar_space, symmetric=True)
    form += (
        trial_axial * test_axial
        + regularization * ng.InnerProduct(ng.grad(trial), ng.grad(test))
        + regularization / aperture_scale**2 * trial * test
    ) * ng.dx
    load = ng.LinearForm(scalar_space)
    load += -ng.InnerProduct(vector_potential, horizontal) * test_axial * ng.dx
    form.Assemble()
    load.Assemble()
    potential = ng.GridFunction(scalar_space, name="chi_axial_ctype_ab")
    potential.vec.data = form.mat.Inverse(
        scalar_space.FreeDofs(), inverse="sparsecholesky"
    ) * load.vec
    return vector_potential + ng.grad(potential)


def centered_states(scale, base_amplitudes):
    states = [np.zeros(6)]
    for index in COMPARABLE_INDICES:
        offset = np.zeros(6)
        offset[index] = float(scale) * float(base_amplitudes[index])
        states.extend((offset, -offset))
    return np.ascontiguousarray(states)


def fit_affine_outputs(states, outputs):
    inputs = states[:, np.asarray(COMPARABLE_INDICES)]
    design = np.column_stack((np.ones(len(states)), inputs))
    coefficients, _, rank, singular_values = np.linalg.lstsq(
        design,
        outputs[:, np.asarray(COMPARABLE_INDICES)],
        rcond=None,
    )
    if rank != design.shape[1]:
        raise RuntimeError("centered tracking basis is rank deficient")
    residuals = outputs[:, np.asarray(COMPARABLE_INDICES)] - design @ coefficients
    return {
        "intercept": coefficients[0],
        "jacobian": coefficients[1:].T,
        "residuals": residuals,
        "rank": int(rank),
        "singular_values": singular_values,
        "maximum_fit_residual": float(np.max(np.abs(residuals), initial=0.0)),
    }


def json_safe(value):
    """Convert NumPy values and non-finite B-route ell entries to strict JSON."""
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (np.floating, float)):
        scalar = float(value)
        return scalar if np.isfinite(scalar) else None
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def evenly_spaced_indices(size, requested):
    count = min(int(size), int(requested))
    if count < 1:
        raise ValueError("three-way B sample counts must be positive")
    return np.unique(
        np.rint(np.linspace(0, int(size) - 1, count)).astype(np.int64)
    )


def _comparison_record(comparison, source_b_map=None):
    record = {
        "initial_state": comparison.a_map.trajectory.canonical_states[0],
        "a_final_state": comparison.a_map.trajectory.final_state,
        "b_final_state": comparison.b_map.final_canonical_state,
        "a_minus_b": comparison.final_state_difference,
        "maximum_a_b_difference": comparison.maximum_final_state_difference,
        "a_field_evaluations": comparison.a_map.field_evaluations,
        "b_field_evaluations": comparison.b_map.field_evaluations,
        "maximum_raw_Ax_t_m": comparison.a_map.maximum_raw_Ax_t_m,
        "b_exit_plane_residual_m": comparison.b_map.exit_plane_residual_m,
        "b_field_representation": comparison.b_map.field_representation,
        "b_independent_variable": comparison.b_map.independent_variable,
    }
    if source_b_map is not None:
        hdiv_minus_source = (
            comparison.b_map.final_canonical_state
            - source_b_map.final_canonical_state
        )
        record.update(
            {
                "source_b_final_state": source_b_map.final_canonical_state,
                "hdiv_b_minus_source_b": hdiv_minus_source,
                "maximum_hdiv_b_source_b_tracking_difference": float(
                    np.max(
                        np.abs(hdiv_minus_source[list(COMPARABLE_INDICES)]),
                        initial=0.0,
                    )
                ),
                "source_b_field_evaluations": source_b_map.field_evaluations,
                "source_b_field_representation": (
                    source_b_map.field_representation
                ),
            }
        )
    return record


def run_case(options, iron_maxh):
    rad.UtiDelAll()
    started = time.perf_counter()
    ng.SetNumThreads(int(options.threads))
    ng.SetHeapSize(int(options.ngsolve_heap_size))
    try:
        solve_started = time.perf_counter()
        with ng.TaskManager():
            iron_mesh = build_iron(
                iron_maxh,
                pole_half_gap_m=options.pole_half_gap,
                pole_half_width_m=options.pole_half_width,
            )
            coil = build_coil(options.current)
            iron = vim.MeshSoftIron(
                iron_mesh,
                bh_table=load_bh_table(),
                order=1,
            )
            model = rad.ObjCnt([iron, coil])
            solution = rad.Solve(model)
        solve_wall = time.perf_counter() - solve_started
        b_point, b_batch = make_symmetric_b_field(solution, coil)
        # Track past both certified planes so the loft tube extends beyond the
        # map region; the boundary-face pollution then lives in the margin.
        margin_stations = int(options.orbit_margin_stations)
        margin_x = margin_stations * 0.080 / (int(options.orbit_stations) - 1)
        extended_orbit = track_reference_orbit(
            b_point,
            options.magnetic_rigidity,
            station_count=int(options.orbit_stations) + 2 * margin_stations,
            entrance_x_m=-0.040 - margin_x,
            exit_x_m=0.040 + margin_x,
        )
        orbit = trim_orbit(extended_orbit, margin_stations)
        midpoint_s = 0.5 * (
            orbit.arc_length_stations[:-1] + orbit.arc_length_stations[1:]
        )
        midpoint_positions = orbit.position_at(midpoint_s)
        _, midpoint_vertical, _ = orbit.frame_at(midpoint_s)
        midpoint_b = b_batch(midpoint_positions)
        curvature_from_b = -np.einsum(
            "ij,ij->i", midpoint_b, midpoint_vertical
        ) / float(options.magnetic_rigidity)
        curvature_discrepancy = orbit.signed_curvature - curvature_from_b
        geometric_curvature_discrepancy = (
            orbit.geometric_signed_curvature - curvature_from_b
        )

        field_started = time.perf_counter()
        with ng.TaskManager():
            tube = build_curvilinear_beam_mesh(
                extended_orbit,
                half_width_m=options.half_width,
                half_height_m=options.half_height,
                maxh_m=options.longitudinal_maxh,
                horizontal_subdivisions_per_macro_strip=(
                    options.horizontal_subdivisions_per_macro_strip
                ),
                vertical_layers=options.vertical_layers,
                curve_order=options.curve_order,
            )
            coil_b = rad.RadiaField(coil, "b")
            u_axis, v_axis, w_axis = reflection_axes(orbit.bend_axis)
            reflected_coil_b = rad.RadiaField(
                coil,
                "b",
                origin=[0.0, 0.0, 0.0],
                u_axis=u_axis,
                v_axis=v_axis,
                w_axis=w_axis,
            )
            raw_magnetic_flux_density_coefficient = (
                MU0 * vim.FieldCoefficientFromSolution(solution)
                + coil_b
            )
            magnetic_flux_density_coefficient = (
                MU0
                * vim.FieldCoefficientFromSolution(
                    solution, reflection_normal=orbit.bend_axis
                )
                # RadiaField's coordinate transform is polar.  Magnetic B is
                # axial, hence the reflected contribution is -R B(R r).
                + 0.5 * (coil_b - reflected_coil_b)
            )
            magnetic_flux_density = ng.GridFunction(
                ng.HDiv(tube.mesh, order=options.hdiv_order),
                name="B_earlytimes_ctype_ab",
            )
            magnetic_flux_density.Set(
                magnetic_flux_density_coefficient, use_simd=False
            )
            if options.a_construction == "exact":
                demag_a = vim.VectorPotentialCoefficientFromSolution(
                    solution,
                    construction="exact",
                    reflection_normal=orbit.bend_axis,
                )
            else:
                demag_a = vim.VectorPotentialCoefficientFromSolution(
                    solution,
                    construction="quadrature",
                    integration_order=options.a_integration_order,
                    reflection_normal=orbit.bend_axis,
                    target_points_m=np.asarray(tube.vertices_m).reshape(-1, 3),
                    maximum_integration_order=options.a_maximum_integration_order,
                )
            coil_a = rad.RadiaField(coil, "a")
            reflected_coil_a = rad.RadiaField(
                coil,
                "a",
                origin=[0.0, 0.0, 0.0],
                u_axis=[1.0, 0.0, 0.0],
                v_axis=[0.0, 1.0, 0.0],
                w_axis=[0.0, 0.0, -1.0],
            )
            current_a = demag_a + 0.5 * (coil_a + reflected_coil_a)
            # Interpolate the expensive analytic source exactly once.  Both
            # gauges add grad(H1(p+1)) corrections, and grad(H1_{p+1}) is a
            # subspace of HCurl_p, so P(A + grad chi) = P(A) + grad chi: the
            # final discrete field is identical whether the projection happens
            # before or after the gauges.  Evaluating the analytic CF in every
            # gauge pass instead of once dominated the whole run (each full
            # pass sums ~5.5e3 analytic kernels at ~5e5 quadrature points).
            projected_once = ng.GridFunction(
                ng.HCurl(tube.mesh, order=options.hcurl_order),
                name="A_projected_once",
            )
            projected_once.Set(current_a, use_simd=False)
            current_a = projected_once
            gauged = None
            for cycle in range(options.gauge_cycles):
                axial_a = project_axial_gauge(
                    current_a, tube.mesh, extended_orbit, options.hcurl_order
                )
                gauged = project_design_orbit_gauge(
                    axial_a,
                    tube,
                    order=options.hcurl_order,
                    constraint_subdivisions=options.gauge_constraint_subdivisions,
                    audit_subdivisions=options.gauge_audit_subdivisions,
                    schur_rcond=options.gauge_schur_rcond,
                    gauge_tolerance=options.gauge_tolerance,
                    name=f"A_orbit_gauge_ctype_ab_{cycle}",
                )
                current_a = gauged.vector_potential
            if gauged is None:
                raise RuntimeError("at least one gauge cycle is required")
            x_offsets = np.linspace(-options.fit_x, options.fit_x, 13)
            y_offsets = np.linspace(-options.fit_y, options.fit_y, 9)
            samples = sample_transverse_vector_potential(
                gauged.vector_potential,
                orbit,
                x_offsets,
                y_offsets,
            )
        field_wall = time.perf_counter() - field_started
        hdiv_midpoint_b = np.asarray(
            [
                magnetic_flux_density(tube.mesh(*point))
                for point in midpoint_positions
            ],
            dtype=float,
        )
        raw_source_midpoint_b = np.asarray(
            [
                raw_magnetic_flux_density_coefficient(tube.mesh(*point))
                for point in midpoint_positions
            ],
            dtype=float,
        )
        symmetric_source_midpoint_b = np.asarray(
            [
                magnetic_flux_density_coefficient(tube.mesh(*point))
                for point in midpoint_positions
            ],
            dtype=float,
        )
        hdiv_source_difference = hdiv_midpoint_b - symmetric_source_midpoint_b
        raw_source_symmetry_difference = (
            raw_source_midpoint_b - symmetric_source_midpoint_b
        )
        symmetric_source_validation_difference = (
            symmetric_source_midpoint_b - midpoint_b
        )
        hdiv_symmetric_source_difference = hdiv_midpoint_b - midpoint_b

        three_way_b_record = None
        if not options.skip_three_way_b:
            three_way_started = time.perf_counter()
            s_indices = evenly_spaced_indices(
                samples.global_points_m.shape[0],
                options.three_way_longitudinal_count,
            )
            x_indices = evenly_spaced_indices(
                samples.global_points_m.shape[1],
                options.three_way_x_count,
            )
            y_indices = evenly_spaced_indices(
                samples.global_points_m.shape[2],
                options.three_way_y_count,
            )
            three_way_points = (
                samples.global_points_m[s_indices]
                [:, x_indices]
                [:, :, y_indices]
                .reshape(-1, 3)
            )
            with ng.TaskManager():
                three_way_b = compare_earlytimes_b_representations(
                    gauged.vector_potential,
                    magnetic_flux_density,
                    b_batch,
                    three_way_points,
                )
            three_way_b_record = {
                "runtime_s": time.perf_counter() - three_way_started,
                "point_count": len(three_way_points),
                "longitudinal_indices": s_indices,
                "x_indices": x_indices,
                "y_indices": y_indices,
                "global_points_m": three_way_b.global_points_m,
                "hdiv_b_t": three_way_b.hdiv_b_t,
                "curl_hcurl_a_t": three_way_b.curl_a_t,
                "rad_fld_source_b_t": three_way_b.source_b_t,
                "hdiv_minus_source_t": three_way_b.hdiv_minus_source_t,
                "curl_a_minus_source_t": (
                    three_way_b.curl_a_minus_source_t
                ),
                "hdiv_minus_curl_a_t": (
                    three_way_b.hdiv_minus_curl_a_t
                ),
                "componentwise_maximum_hdiv_minus_source_t": (
                    three_way_b.componentwise_maximum_hdiv_minus_source_t
                ),
                "componentwise_maximum_curl_a_minus_source_t": (
                    three_way_b.componentwise_maximum_curl_a_minus_source_t
                ),
                "componentwise_maximum_hdiv_minus_curl_a_t": (
                    three_way_b.componentwise_maximum_hdiv_minus_curl_a_t
                ),
                "maximum_hdiv_minus_source_vector_norm_t": (
                    three_way_b.maximum_hdiv_minus_source_vector_norm_t
                ),
                "maximum_curl_a_minus_source_vector_norm_t": (
                    three_way_b.maximum_curl_a_minus_source_vector_norm_t
                ),
                "maximum_hdiv_minus_curl_a_vector_norm_t": (
                    three_way_b.maximum_hdiv_minus_curl_a_vector_norm_t
                ),
                "field_roles": {
                    "hdiv_b_t": "EarlyTimes HDiv(order=4) B GridFunction",
                    "curl_hcurl_a_t": (
                        "diagnostic curl of EarlyTimes HCurl(order=5) A "
                        "GridFunction"
                    ),
                    "rad_fld_source_b_t": (
                        "independent HDiv-MMM plus rad.Fld source truth"
                    ),
                },
            }

        derivative_records = []
        derivative_initial = np.asarray(options.derivative_initial, dtype=float)
        if not options.skip_a_b_track:
            for step in options.transverse_derivative_steps:
                with ng.TaskManager():
                    comparison = compare_hcurl_a_map_to_hdiv_b_map_rk(
                        gauged.vector_potential,
                        magnetic_flux_density,
                        orbit,
                        derivative_initial,
                        curvature_sign=options.curvature_sign,
                        maximum_step_m=options.maximum_track_step,
                        exit_plane_tolerance_m=options.exit_plane_tolerance,
                        axial_gauge_tolerance_t_m=options.axial_gauge_tolerance,
                        transverse_derivative_step_m=step,
                    )
                    source_b_map = track_b_coefficient_cartesian_s(
                        magnetic_flux_density_coefficient,
                        gauged.vector_potential,
                        orbit,
                        derivative_initial,
                        field_mesh=tube.mesh,
                        curvature_sign=options.curvature_sign,
                        maximum_step_m=options.maximum_track_step,
                        exit_plane_tolerance_m=options.exit_plane_tolerance,
                        axial_gauge_tolerance_t_m=options.axial_gauge_tolerance,
                        transverse_derivative_step_m=step,
                    )
                record = _comparison_record(comparison, source_b_map)
                record["transverse_derivative_step_m"] = float(step)
                derivative_records.append(record)

        a_lie_record = None
        if options.include_a_lie:
            lie_started = time.perf_counter()
            a_lie_states = np.ascontiguousarray(
                [
                    float(scale) * derivative_initial
                    for scale in options.a_lie_scales
                ]
            )
            with ng.TaskManager():
                lie_map = fourth_order_lie_map_from_hcurl_transverse(
                    gauged.vector_potential,
                    orbit,
                    x_offsets,
                    y_offsets,
                    field_symmetry="normal",
                    curvature_sign=options.curvature_sign,
                    axial_gauge_tolerance_t_m=(
                        options.axial_gauge_tolerance
                    ),
                    orbit_gauge_tolerance_t_m=options.gauge_tolerance,
                    symmetry_tolerance_t_m=options.gauge_tolerance,
                    reference_orbit_tolerance=(
                        options.lie_reference_orbit_tolerance
                    ),
                    maximum_step_m=options.lie_map_step,
                    factorization_tolerance=(
                        options.lie_factorization_tolerance
                    ),
                )
                a_lie = compare_hcurl_lie_map_to_a_rk(
                    gauged.vector_potential,
                    lie_map,
                    a_lie_states,
                    generator_substeps=options.lie_generator_substeps,
                    curvature_sign=options.curvature_sign,
                    maximum_step_m=options.maximum_track_step,
                    axial_gauge_tolerance_t_m=(
                        options.axial_gauge_tolerance
                    ),
                    transverse_derivative_step_m=(
                        options.basis_derivative_step
                    ),
                )
            transfer = lie_map.lie_map.transfer
            fit = lie_map.polynomial_fit
            a_lie_record = {
                "runtime_s": time.perf_counter() - lie_started,
                "initial_states": a_lie_states,
                "maximum_lie_minus_a_rk": (
                    a_lie.maximum_lie_truncation_error
                ),
                "componentwise_maximum_lie_minus_a_rk": (
                    a_lie.componentwise_maximum_lie_truncation_error
                ),
                "maximum_polynomial_minus_a_rk": (
                    a_lie.maximum_polynomial_truncation_error
                ),
                "maximum_lie_minus_polynomial": (
                    a_lie.maximum_lie_polynomial_difference
                ),
                "cases": [
                    {
                        "initial_state": case.initial_state,
                        "lie_state": case.lie_state,
                        "polynomial_state": case.polynomial_state,
                        "a_rk_state": case.a_rk_state,
                        "lie_minus_a_rk": case.lie_minus_a_rk,
                        "polynomial_minus_a_rk": case.polynomial_minus_a_rk,
                        "lie_minus_polynomial": case.lie_minus_polynomial,
                        "a_field_evaluations": case.a_field_evaluations,
                        "a_rk_accepted_steps": case.a_rk_accepted_steps,
                        "maximum_raw_Ax_t_m": case.maximum_raw_Ax_t_m,
                        "lie_generator_substeps": (
                            case.lie_generator_substeps
                        ),
                    }
                    for case in a_lie.cases
                ],
                "map": {
                    "maximum_hamiltonian_linear": float(
                        np.max(
                            np.abs(lie_map.lie_map.hamiltonian_linear),
                            initial=0.0,
                        )
                    ),
                    "maximum_Ay_fit_residual_t_m": (
                        fit.maximum_Ay_fit_residual_t_m
                    ),
                    "maximum_As_fit_residual_t_m": (
                        fit.maximum_As_fit_residual_t_m
                    ),
                    "maximum_left_right_scaled_coefficient_discrepancy_t_m": (
                        fit.maximum_left_right_scaled_coefficient_discrepancy_t_m
                    ),
                    "maximum_generator_symmetry_defect": (
                        transfer.factorization.maximum_generator_symmetry_defect
                    ),
                    "relative_reconstruction_error": (
                        transfer.factorization.relative_reconstruction_error
                    ),
                    "reconstructed_symplectic_residual": (
                        transfer.factorization.reconstructed_symplectic_residual.maximum
                    ),
                },
            }

        gauge_triangle_record = None
        if options.include_a_lie and a_lie_record is not None:
            # Three-route agreement in gauge-invariant observables:
            #   1. gauged-A LIE (fourth-order Dragt--Finn),
            #   2. UNGAUGED (axial-gauge-only) A exact canonical RK,
            #   3. HDiv B-map Cartesian RK.
            # Canonical momenta are gauge dependent, so every route is
            # converted to mechanical momenta with its OWN vector potential at
            # the entrance and exit; x, y, delta, and ell are gauge invariant
            # (the gauge scalar does not depend on delta).  The gauged-A RK is
            # retained as the bridge: its difference from the ungauged-A RK is
            # a pure gauge-invariance check of the canonical machinery.
            triangle_started = time.perf_counter()
            normalization = float(options.curvature_sign) / float(
                options.magnetic_rigidity
            )
            s_entry = float(orbit.arc_length_stations[0])
            s_exit = float(orbit.arc_length_stations[-1])
            ungauged_a = gauged.ungauged_vector_potential

            def frame_A(field, s_value, x_value, y_value):
                horizontal, vertical, tangent_ = orbit.frame_at(
                    np.array([s_value])
                )
                point = (
                    orbit.position_at(np.array([s_value]))[0]
                    + x_value * horizontal[0]
                    + y_value * vertical[0]
                )
                value = np.asarray(field(tube.mesh(*point)), dtype=float)
                return np.asarray(
                    [
                        float(np.dot(value, horizontal[0])),
                        float(np.dot(value, vertical[0])),
                        float(np.dot(value, tangent_[0])),
                    ]
                )

            def to_mechanical(state, field, s_value):
                local_A = frame_A(field, s_value, float(state[0]), float(state[2]))
                mechanical = np.asarray(state, dtype=float).copy()
                mechanical[1] = float(state[1]) - normalization * local_A[0]
                mechanical[3] = float(state[3]) - normalization * local_A[1]
                return mechanical, local_A

            def to_canonical(mechanical, field, s_value):
                local_A = frame_A(
                    field, s_value, float(mechanical[0]), float(mechanical[2])
                )
                canonical = np.asarray(mechanical, dtype=float).copy()
                canonical[1] = float(mechanical[1]) + normalization * local_A[0]
                canonical[3] = float(mechanical[3]) + normalization * local_A[1]
                return canonical, local_A

            comparable = list(COMPARABLE_INDICES)
            triangle_cases = []
            for case in a_lie.cases:
                initial = np.asarray(case.initial_state, dtype=float)
                mechanical_entry, entry_A_gauged = to_mechanical(
                    initial, gauged.vector_potential, s_entry
                )
                ungauged_initial, entry_A_ungauged = to_canonical(
                    mechanical_entry, ungauged_a, s_entry
                )
                with ng.TaskManager():
                    ungauged_map = track_hcurl_vector_potential_canonical_s(
                        ungauged_a,
                        orbit,
                        ungauged_initial,
                        curvature_sign=options.curvature_sign,
                        maximum_step_m=options.maximum_track_step,
                        axial_gauge_tolerance_t_m=options.axial_gauge_tolerance,
                        transverse_derivative_step_m=(
                            options.basis_derivative_step
                        ),
                    )
                    b_map = track_hdiv_b_map_cartesian_s(
                        magnetic_flux_density,
                        gauged.vector_potential,
                        orbit,
                        initial,
                        curvature_sign=options.curvature_sign,
                        maximum_step_m=options.maximum_track_step,
                        exit_plane_tolerance_m=options.exit_plane_tolerance,
                        axial_gauge_tolerance_t_m=options.axial_gauge_tolerance,
                        transverse_derivative_step_m=(
                            options.basis_derivative_step
                        ),
                    )
                gauged_mech, _ = to_mechanical(
                    case.a_rk_state, gauged.vector_potential, s_exit
                )
                lie_mech, _ = to_mechanical(
                    case.lie_state, gauged.vector_potential, s_exit
                )
                ungauged_mech, _ = to_mechanical(
                    ungauged_map.trajectory.final_state, ungauged_a, s_exit
                )
                b_mech, _ = to_mechanical(
                    b_map.final_canonical_state, gauged.vector_potential, s_exit
                )
                gauge_invariance = ungauged_mech - gauged_mech
                lie_minus_a = lie_mech - gauged_mech
                a_minus_b = gauged_mech[comparable] - b_mech[comparable]
                lie_minus_b = lie_mech[comparable] - b_mech[comparable]
                triangle_cases.append(
                    {
                        "initial_canonical_gauged": initial,
                        "initial_mechanical": mechanical_entry,
                        "initial_canonical_ungauged": ungauged_initial,
                        "entry_A_gauged_Tm": entry_A_gauged,
                        "entry_A_ungauged_Tm": entry_A_ungauged,
                        "gauged_a_rk_mechanical_exit": gauged_mech,
                        "ungauged_a_rk_mechanical_exit": ungauged_mech,
                        "lie_mechanical_exit": lie_mech,
                        "b_rk_mechanical_exit": b_mech,
                        "gauge_invariance_ungauged_minus_gauged": (
                            gauge_invariance
                        ),
                        "maximum_gauge_invariance_defect": float(
                            np.max(np.abs(gauge_invariance[comparable + [4]]))
                        ),
                        "lie_minus_gauged_a_rk": lie_minus_a,
                        "maximum_lie_minus_gauged_a_rk": float(
                            np.max(np.abs(lie_minus_a[comparable + [4]]))
                        ),
                        "gauged_a_rk_minus_b_rk": a_minus_b,
                        "maximum_gauged_a_rk_minus_b_rk": float(
                            np.max(np.abs(a_minus_b))
                        ),
                        "lie_minus_b_rk": lie_minus_b,
                        "maximum_lie_minus_b_rk": float(
                            np.max(np.abs(lie_minus_b))
                        ),
                        "ungauged_maximum_raw_Ax_Tm": float(
                            ungauged_map.maximum_raw_Ax_t_m
                        ),
                        "ungauged_field_evaluations": int(
                            ungauged_map.field_evaluations
                        ),
                        "b_exit_plane_residual_m": float(
                            b_map.exit_plane_residual_m
                        ),
                    }
                )
            gauge_triangle_record = {
                "runtime_s": time.perf_counter() - triangle_started,
                "mechanical_component_order": COMPONENT_NAMES,
                "cases": triangle_cases,
                "maximum_gauge_invariance_defect": float(
                    max(
                        case["maximum_gauge_invariance_defect"]
                        for case in triangle_cases
                    )
                ),
                "maximum_lie_minus_gauged_a_rk": float(
                    max(
                        case["maximum_lie_minus_gauged_a_rk"]
                        for case in triangle_cases
                    )
                ),
                "maximum_gauged_a_rk_minus_b_rk": float(
                    max(
                        case["maximum_gauged_a_rk_minus_b_rk"]
                        for case in triangle_cases
                    )
                ),
                "maximum_lie_minus_b_rk": float(
                    max(
                        case["maximum_lie_minus_b_rk"]
                        for case in triangle_cases
                    )
                ),
            }

        aperture_records = []
        if not options.skip_centered_basis:
            all_states = [np.zeros(6)]
            scale_case_indices = {}
            for scale in options.aperture_scales:
                indices = [0]
                for state in centered_states(scale, options.basis_amplitudes)[1:]:
                    indices.append(len(all_states))
                    all_states.append(state)
                scale_case_indices[float(scale)] = indices
            all_states = np.ascontiguousarray(all_states)
            with ng.TaskManager():
                batch = compare_hcurl_a_map_to_hdiv_b_map_rk_batch(
                    gauged.vector_potential,
                    magnetic_flux_density,
                    orbit,
                    all_states,
                    curvature_sign=options.curvature_sign,
                    maximum_step_m=options.maximum_track_step,
                    exit_plane_tolerance_m=options.exit_plane_tolerance,
                    axial_gauge_tolerance_t_m=options.axial_gauge_tolerance,
                    transverse_derivative_step_m=options.basis_derivative_step,
                )

            for scale in options.aperture_scales:
                indices = scale_case_indices[float(scale)]
                subset = HCurlAMapBCoefficientRKBatchComparison(
                    initial_states=batch.initial_states[indices],
                    cases=tuple(batch.cases[index] for index in indices),
                )
                difference_fit = subset.fit_affine_route_difference()
                a_outputs = np.vstack(
                    [case.a_map.trajectory.final_state for case in subset.cases]
                )
                b_outputs = np.vstack(
                    [case.b_map.final_canonical_state for case in subset.cases]
                )
                a_fit = fit_affine_outputs(subset.initial_states, a_outputs)
                b_fit = fit_affine_outputs(subset.initial_states, b_outputs)
                route_difference_check = (
                    a_fit["jacobian"]
                    - b_fit["jacobian"]
                    - difference_fit.jacobian
                )
                aperture_records.append(
                    {
                        "scale": float(scale),
                        "basis_amplitudes": (
                            float(scale) * np.asarray(options.basis_amplitudes)
                        ),
                        "initial_states": subset.initial_states,
                        "a_minus_b_final_states": subset.final_state_differences,
                        "maximum_a_b_final_state_difference": (
                            subset.maximum_final_state_difference
                        ),
                        "componentwise_maximum_a_b_final_state_difference": (
                            subset.componentwise_maximum_final_state_difference
                        ),
                        "a_minus_b_intercept": difference_fit.intercept,
                        "a_minus_b_jacobian": difference_fit.jacobian,
                        "a_minus_b_maximum_intercept": (
                            difference_fit.maximum_intercept
                        ),
                        "a_minus_b_maximum_jacobian_entry": (
                            difference_fit.maximum_jacobian_entry
                        ),
                        "a_minus_b_maximum_affine_fit_residual": (
                            difference_fit.maximum_fit_residual
                        ),
                        "a_route_intercept": a_fit["intercept"],
                        "a_route_jacobian": a_fit["jacobian"],
                        "a_route_maximum_affine_fit_residual": (
                            a_fit["maximum_fit_residual"]
                        ),
                        "b_route_intercept": b_fit["intercept"],
                        "b_route_jacobian": b_fit["jacobian"],
                        "b_route_maximum_affine_fit_residual": (
                            b_fit["maximum_fit_residual"]
                        ),
                        "maximum_route_difference_identity_residual": float(
                            np.max(np.abs(route_difference_check), initial=0.0)
                        ),
                    }
                )

        orbit_x_index, orbit_y_index = samples.design_orbit_indices
        orbit_Ay = samples.Ay_t_m[:, orbit_x_index, orbit_y_index]
        orbit_As = samples.As_t_m[:, orbit_x_index, orbit_y_index]
        source_stats = list(
            solution["vector_potential_coefficient_stats"].values()
        )[-1]
        return {
            "schema": "radia.validation.earlytimes-ctype-ab.v1",
            "status": "pass",
            "hostname": platform.node(),
            "runtime_s": time.perf_counter() - started,
            "solve_wall_s": solve_wall,
            "field_map_wall_s": field_wall,
            "settings": {
                "iron_maxh_m": float(iron_maxh),
                "current_ampere_turns": float(options.current),
                "magnetic_rigidity_t_m": float(options.magnetic_rigidity),
                "orbit_stations": int(options.orbit_stations),
                "orbit_margin_stations": int(options.orbit_margin_stations),
                "half_width_m": float(options.half_width),
                "half_height_m": float(options.half_height),
                "aperture_profile": (
                    "wide-20mm" if options.wide_20mm else "explicit-or-legacy"
                ),
                "pole_half_gap_m": float(options.pole_half_gap),
                "pole_half_width_m": float(options.pole_half_width),
                "longitudinal_maxh_m": float(options.longitudinal_maxh),
                "vertical_layers": int(options.vertical_layers),
                "horizontal_subdivisions_per_macro_strip": int(
                    options.horizontal_subdivisions_per_macro_strip
                ),
                "curve_order": int(options.curve_order),
                "hcurl_order": int(options.hcurl_order),
                "hdiv_order": int(options.hdiv_order),
                "a_construction": str(options.a_construction),
                "a_integration_order": int(options.a_integration_order),
                "a_maximum_integration_order": int(
                    options.a_maximum_integration_order
                ),
                "basis_derivative_step_m": float(options.basis_derivative_step),
                "maximum_track_step_m": float(options.maximum_track_step),
                "aperture_scales": options.aperture_scales,
                "a_lie_scales": options.a_lie_scales,
                "lie_map_step_m": float(options.lie_map_step),
                "lie_generator_substeps": int(
                    options.lie_generator_substeps
                ),
                "three_way_longitudinal_count": int(
                    options.three_way_longitudinal_count
                ),
                "three_way_x_count": int(options.three_way_x_count),
                "three_way_y_count": int(options.three_way_y_count),
                "component_order": COMPONENT_NAMES,
                "comparable_indices": COMPARABLE_INDICES,
            },
            "field": {
                "iron_elements": int(iron_mesh.ne),
                "hdiv_dofs": int(solution["ndof"]),
                "nonlinear_iterations": int(solution["iters"]),
                "beam_hexes": int(tube.mesh.ne),
                "orbit_length_m": float(orbit.length_m),
                "orbit_exit_m": orbit.positions[-1],
                "B_center_T": b_point(orbit.position_at(0.5 * orbit.length_m)),
                "maximum_curvature_per_m": float(
                    np.max(np.abs(orbit.signed_curvature), initial=0.0)
                ),
                "maximum_orbit_vs_B_midpoint_curvature_discrepancy_per_m": float(
                    np.max(np.abs(curvature_discrepancy), initial=0.0)
                ),
                "maximum_geometric_vs_B_midpoint_curvature_discrepancy_per_m": float(
                    np.max(
                        np.abs(geometric_curvature_discrepancy), initial=0.0
                    )
                ),
                "maximum_HDiv_B_source_projection_difference_T": float(
                    np.max(np.abs(hdiv_source_difference), initial=0.0)
                ),
                "maximum_raw_source_B_symmetrization_difference_T": float(
                    np.max(
                        np.abs(raw_source_symmetry_difference), initial=0.0
                    )
                ),
                "maximum_HDiv_B_symmetric_source_difference_T": float(
                    np.max(
                        np.abs(hdiv_symmetric_source_difference), initial=0.0
                    )
                ),
                "maximum_symmetric_source_B_validation_difference_T": float(
                    np.max(
                        np.abs(symmetric_source_validation_difference),
                        initial=0.0,
                    )
                ),
                "HDiv_B_midpoint_T": hdiv_midpoint_b,
                "raw_source_B_midpoint_T": raw_source_midpoint_b,
                "symmetric_source_B_midpoint_T": symmetric_source_midpoint_b,
                "independent_symmetric_source_B_midpoint_T": midpoint_b,
                "A_source_stats": source_stats,
                "orbit_gauge_residual_Tm": float(
                    gauged.maximum_orbit_gauge_residual_t_m
                ),
                "curl_change_l2_Tm32": float(gauged.curl_change_l2_t_m32),
                "maximum_Ax_Tm": float(
                    np.max(np.abs(samples.Ax_t_m), initial=0.0)
                ),
                "maximum_orbit_Ay_As_Tm": float(
                    max(
                        np.max(np.abs(orbit_Ay), initial=0.0),
                        np.max(np.abs(orbit_As), initial=0.0),
                    )
                ),
            },
            "derivative_step_study": derivative_records,
            "aperture_scale_study": aperture_records,
            "b_representation_comparison": three_way_b_record,
            "a_lie_comparison": a_lie_record,
            "gauge_invariance_triangle": gauge_triangle_record,
        }
    finally:
        rad.UtiDelAll()


def write_result(path, result):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(json_safe(result), indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def parser():
    result = argparse.ArgumentParser()
    result.add_argument("--output-json", type=Path, required=True)
    result.add_argument("--overwrite", action="store_true")
    result.add_argument("--source-revision", default="working-tree")
    result.add_argument("--threads", type=int, default=8)
    result.add_argument("--ngsolve-heap-size", type=int, default=10_000_000)
    result.add_argument("--iron-maxh-values", type=float, nargs="+", default=[0.03])
    result.add_argument("--current", type=float, default=3000.0)
    result.add_argument("--magnetic-rigidity", type=float, default=3.0)
    result.add_argument("--orbit-stations", type=int, default=33)
    result.add_argument(
        "--orbit-margin-stations",
        type=int,
        default=2,
        help=(
            "stations tracked past each certified plane; the loft tube covers "
            "them so the boundary-face fit pollution stays outside the map "
            "region"
        ),
    )
    result.add_argument("--half-width", type=float, default=0.010)
    result.add_argument("--half-height", type=float, default=0.0035)
    result.add_argument("--pole-half-gap", type=float, default=0.005)
    result.add_argument("--pole-half-width", type=float, default=0.020)
    result.add_argument(
        "--wide-20mm",
        action="store_true",
        help=(
            "use the physical +/-20 mm aperture preset: 50 mm pole gap, "
            "60 mm pole width, p6-ready 8-by-9 transverse hp mesh"
        ),
    )
    result.add_argument("--fit-x", type=float, default=0.008)
    result.add_argument("--fit-y", type=float, default=0.0025)
    result.add_argument("--longitudinal-maxh", type=float, default=0.005)
    result.add_argument(
        "--horizontal-subdivisions-per-macro-strip", type=int, default=1
    )
    result.add_argument("--vertical-layers", type=int, default=1)
    result.add_argument("--curve-order", type=int, default=2)
    result.add_argument("--hcurl-order", type=int, default=5)
    result.add_argument("--hdiv-order", type=int, default=4)
    result.add_argument(
        "--a-construction",
        choices=("exact", "quadrature"),
        default="exact",
        help=(
            "A source representation: 'exact' uses the analytic "
            "equivalent-current tetrahedron/triangle kernels that match the B "
            "route; 'quadrature' reproduces the recorded q14/q20/q24 "
            "point-dipole studies"
        ),
    )
    result.add_argument("--a-integration-order", type=int, default=6)
    result.add_argument("--a-maximum-integration-order", type=int, default=14)
    result.add_argument("--gauge-cycles", type=int, default=3)
    result.add_argument(
        "--gauge-constraint-subdivisions", type=int, default=12
    )
    result.add_argument("--gauge-audit-subdivisions", type=int, default=31)
    result.add_argument("--gauge-schur-rcond", type=float, default=1.0e-12)
    result.add_argument("--gauge-tolerance", type=float, default=3.0e-6)
    result.add_argument("--axial-gauge-tolerance", type=float, default=1.0e-6)
    result.add_argument("--curvature-sign", type=float, default=1.0)
    result.add_argument("--maximum-track-step", type=float, default=5.0e-4)
    result.add_argument("--exit-plane-tolerance", type=float, default=1.0e-10)
    result.add_argument(
        "--transverse-derivative-steps",
        type=float,
        nargs="+",
        default=[2.0e-6, 1.0e-6, 5.0e-7],
    )
    result.add_argument("--basis-derivative-step", type=float, default=1.0e-6)
    result.add_argument(
        "--derivative-initial",
        type=float,
        nargs=6,
        default=[2.0e-4, 1.0e-4, -1.0e-4, -5.0e-5, 0.0, 1.0e-4],
    )
    result.add_argument(
        "--basis-amplitudes",
        type=float,
        nargs=6,
        default=[2.0e-4, 2.0e-4, 1.0e-4, 1.0e-4, 0.0, 1.0e-4],
    )
    result.add_argument(
        "--aperture-scales", type=float, nargs="+", default=[0.5, 1.0, 2.0]
    )
    result.add_argument("--skip-centered-basis", action="store_true")
    result.add_argument("--skip-a-b-track", action="store_true")
    result.add_argument("--skip-three-way-b", action="store_true")
    result.add_argument(
        "--three-way-longitudinal-count", type=int, default=5
    )
    result.add_argument("--three-way-x-count", type=int, default=3)
    result.add_argument("--three-way-y-count", type=int, default=3)
    result.add_argument("--include-a-lie", action="store_true")
    result.add_argument(
        "--a-lie-scales", type=float, nargs="+", default=[1.0]
    )
    result.add_argument("--lie-map-step", type=float, default=5.0e-4)
    result.add_argument("--lie-generator-substeps", type=int, default=4)
    result.add_argument(
        "--lie-reference-orbit-tolerance", type=float, default=3.0e-6
    )
    result.add_argument(
        "--lie-factorization-tolerance", type=float, default=1.0e-6
    )
    return result


def configure_options(options):
    """Apply the named physical-aperture preset and validate its clearances."""
    if options.wide_20mm:
        options.half_width = 0.020
        options.half_height = 0.020
        options.fit_x = 0.018
        options.fit_y = 0.018
        options.pole_half_gap = 0.025
        options.pole_half_width = 0.030
        options.horizontal_subdivisions_per_macro_strip = 2
        options.vertical_layers = 9
        options.hcurl_order = max(int(options.hcurl_order), 6)
    if not 0.0 < float(options.half_height) < float(options.pole_half_gap):
        raise ValueError(
            "half-height must be positive and smaller than the pole half-gap"
        )
    if not 0.0 < float(options.half_width) < float(options.pole_half_width):
        raise ValueError(
            "half-width must be positive and smaller than the pole half-width"
        )
    if not 0.0 < float(options.fit_x) <= float(options.half_width):
        raise ValueError("fit-x must lie inside the beam half-width")
    if not 0.0 < float(options.fit_y) <= float(options.half_height):
        raise ValueError("fit-y must lie inside the beam half-height")
    return options


def main():
    options = configure_options(parser().parse_args())
    if options.output_json.exists() and not options.overwrite:
        raise FileExistsError(f"refusing to overwrite {options.output_json}")
    started = time.perf_counter()
    document = {
        "schema": "radia.validation.earlytimes-ctype-ab-study.v1",
        "status": "running",
        "hostname": platform.node(),
        "python": platform.python_version(),
        "ngsolve": getattr(ng, "__version__", "unknown"),
        "source_revision": options.source_revision,
        "cases": [],
    }
    write_result(options.output_json, document)
    try:
        for iron_maxh in options.iron_maxh_values:
            print(f"starting iron_maxh={iron_maxh:.6g} m", flush=True)
            case = run_case(options, iron_maxh)
            document["cases"].append(case)
            write_result(options.output_json, document)
            try:
                print(
                    f"completed iron_maxh={iron_maxh:.6g} m "
                    f"in {case['runtime_s']:.3f} s",
                    flush=True,
                )
            except OSError:
                pass
    except Exception as error:
        document["status"] = "failed"
        document["runtime_s"] = time.perf_counter() - started
        document["error"] = f"{type(error).__name__}: {error}"
        document["traceback"] = traceback.format_exc()
        write_result(options.output_json, document)
        raise
    document["status"] = "pass"
    document["runtime_s"] = time.perf_counter() - started
    write_result(options.output_json, document)
    try:
        print(
            json.dumps(
                {
                    "status": document["status"],
                    "hostname": document["hostname"],
                    "runtime_s": document["runtime_s"],
                    "case_count": len(document["cases"]),
                    "output_json": str(options.output_json),
                },
                indent=2,
            )
        )
    except OSError:
        pass


if __name__ == "__main__":
    main()
