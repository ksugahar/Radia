"""Air-volume-mesh-free HDiv-MMM proof of concept for FFAG cells.

The module turns a momentum-indexed family of periodic design orbits and
first-order cell maps into the multi-orbit contract in
``accelerator_magnet_topopt``.  The Bell--Abell non-scaling FFAG parameters
provide a reproducible soft-edge fixture, not a claim to reproduce their PTC
ring: their complete placement and closed-orbit files are not published in the
paper.

The optimization path uses exact combined-function matrix-exponential Frechet
derivatives.  Enge ``I1``/``I2`` values are diagnostics for a reduced edge map;
they are never applied on top of a field response which already samples the
soft fringe.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .accelerator_magnet_topopt import (
    CoilBuilderHDivSource,
    CoilHDivTotalField,
    MultiMomentumAcceleratorMagnetTopologyResult,
    MultiMomentumTransferMatrixObjective,
    PlanarDesignOrbit,
    build_multi_orbit_field_response_matrix,
    optimize_hdiv_mmm_magnet_from_transfer_matrices,
    planar_orbit_field_observations,
    solve_transfer_matrix_field_correction,
)
from .isochronous_topopt import (
    CombinedFunctionTransferMap,
    combined_function_transfer_map_from_field_response,
)


PROTON_REST_ENERGY_MEV = 938.27208816
GEV_C_PER_TESLA_METRE = 0.299792458


def magnetic_rigidity_from_kinetic_energy(
        kinetic_energy_mev, *, rest_energy_mev=PROTON_REST_ENERGY_MEV,
        charge_number=1.0):
    """Return relativistic magnetic rigidity ``B rho`` in tesla-metre."""
    kinetic = np.asarray(kinetic_energy_mev, dtype=float)
    rest = float(rest_energy_mev)
    charge = abs(float(charge_number))
    if (not np.all(np.isfinite(kinetic)) or np.any(kinetic < 0.0)
            or not np.isfinite(rest) or rest <= 0.0
            or not np.isfinite(charge) or charge <= 0.0):
        raise ValueError(
            "kinetic energy must be nonnegative and rest energy/charge "
            "must be positive")
    momentum_gev_c = np.sqrt(
        kinetic * (kinetic + 2.0 * rest)) / 1000.0
    result = momentum_gev_c / (GEV_C_PER_TESLA_METRE * charge)
    return float(result) if result.ndim == 0 else result


@dataclass(frozen=True)
class EngeFringeIntegrals:
    """Equal-integral boundary and dimensionless Enge form factors."""

    effective_boundary_m: float
    i1: float
    i2: float
    equal_integral_residual: float
    full_gap_m: float


def enge_fringe_integrals(
        coordinate_m, bending_field_t, *, body_field_t, full_gap_m
        ) -> EngeFringeIntegrals:
    """Evaluate exit-fringe ``I1`` and ``I2`` from a sampled field.

    Samples must run from the constant-field side to the zero-field side.
    ``sigma=(z-z_eff)/g`` uses the full gap ``g``.  The effective boundary
    ``z_eff`` equates the soft-field integral to a sharp unit step.  With
    ``f=B/B0`` and ``q=H(-sigma)-f``, the numerically stable first-moment form

    ``I1 = - integral sigma*q d sigma``

    is equivalent to the conventional nested integral when
    ``integral q d sigma=0``.  ``I2=integral f*(1-f) d sigma``.
    """
    coordinate = np.asarray(coordinate_m, dtype=float).reshape(-1)
    field = np.asarray(bending_field_t, dtype=float).reshape(-1)
    body = float(body_field_t)
    gap = float(full_gap_m)
    if (coordinate.size < 8 or field.shape != coordinate.shape
            or not np.all(np.isfinite(np.r_[coordinate, field, body, gap]))
            or np.any(np.diff(coordinate) <= 0.0) or body == 0.0
            or gap <= 0.0):
        raise ValueError(
            "Enge integration needs at least eight ordered finite samples, "
            "nonzero body field, and positive full gap")
    normalized = field / body
    end_tolerance = 5.0e-3
    if (abs(normalized[0] - 1.0) > end_tolerance
            or abs(normalized[-1]) > end_tolerance):
        raise ValueError(
            "samples must start in the body field and end in the zero field")
    if np.min(normalized) < -0.25 or np.max(normalized) > 1.25:
        raise ValueError(
            "normalized fringe field is too far outside the body/zero range")
    integral = float(np.trapezoid(normalized, coordinate))
    effective = float(coordinate[0] + integral)
    sigma = (coordinate - effective) / gap
    # Integrate the discontinuous sharp step analytically.  Sampling it with
    # a trapezoid would leave a grid-dependent half-cell in both I1 and the
    # equal-integral residual when z_eff happens to be a sample location.
    hard_moment = -0.5 * float(sigma[0] * sigma[0])
    soft_moment = float(np.trapezoid(sigma * normalized, sigma))
    i1 = soft_moment - hard_moment
    i2 = float(np.trapezoid(normalized * (1.0 - normalized), sigma))
    residual = float(
        ((effective - coordinate[0]) - integral) / gap)
    return EngeFringeIntegrals(
        effective, i1, i2, residual, gap)


def _tanh_window(coordinate, start, stop, epsilon):
    return 0.5 * (
        np.tanh((coordinate - start) / epsilon)
        - np.tanh((coordinate - stop) / epsilon))


@dataclass(frozen=True)
class FFAGSoftEdgeCellSpec:
    """Periodic doublet-cell parameters for a soft-edge FFAG PoC.

    The default constructor is generic.  :meth:`bell_abell` supplies Table 1
    of Bell and Abell, arXiv:1202.0805.  ``full_gap_m`` is deliberately a PoC
    input because that paper states that the 5 cm Enge scale is similar to the
    aperture but does not publish a full pole gap.
    """

    cell_count: int
    long_drift_m: float
    defocusing_length_m: float
    short_drift_m: float
    focusing_length_m: float
    defocusing_b0_t: float
    defocusing_gradient_t_per_m: float
    focusing_b0_t: float
    focusing_gradient_t_per_m: float
    fringe_epsilon_m: float
    full_gap_m: float

    def __post_init__(self):
        count = int(self.cell_count)
        lengths = np.asarray([
            self.long_drift_m, self.defocusing_length_m,
            self.short_drift_m, self.focusing_length_m,
            self.fringe_epsilon_m, self.full_gap_m], dtype=float)
        fields = np.asarray([
            self.defocusing_b0_t, self.defocusing_gradient_t_per_m,
            self.focusing_b0_t, self.focusing_gradient_t_per_m], dtype=float)
        if (count < 2 or not np.all(np.isfinite(lengths))
                or np.any(lengths <= 0.0) or not np.all(np.isfinite(fields))
                or (self.defocusing_gradient_t_per_m
                    + self.focusing_gradient_t_per_m) == 0.0):
            raise ValueError("invalid FFAG soft-edge cell specification")
        object.__setattr__(self, "cell_count", count)

    @classmethod
    def bell_abell(cls, *, full_gap_m=0.10):
        """Return the published 24-cell, 31--250 MeV proton fixture."""
        return cls(
            cell_count=24,
            long_drift_m=0.40,
            defocusing_length_m=0.22,
            short_drift_m=0.075,
            focusing_length_m=0.44,
            defocusing_b0_t=0.803952,
            defocusing_gradient_t_per_m=-12.8,
            focusing_b0_t=0.555057,
            focusing_gradient_t_per_m=8.0,
            fringe_epsilon_m=0.05,
            full_gap_m=full_gap_m)

    @property
    def cell_length_m(self) -> float:
        return float(
            self.long_drift_m + self.defocusing_length_m
            + self.short_drift_m + self.focusing_length_m)

    @property
    def cell_bend_angle_rad(self) -> float:
        return float(2.0 * np.pi / self.cell_count)

    @property
    def magnet_intervals_m(self):
        bd_start = 0.5 * self.long_drift_m
        bd_stop = bd_start + self.defocusing_length_m
        bf_start = bd_stop + self.short_drift_m
        bf_stop = bf_start + self.focusing_length_m
        return ((bd_start, bd_stop), (bf_start, bf_stop))

    def sampled_profiles(self, *, n_segments=256, periodic_images=2):
        """Return midpoint ``s,ds,B0(s),G(s)`` with overlapping fringes."""
        count = int(n_segments)
        images = int(periodic_images)
        if count < 16 or images < 1:
            raise ValueError(
                "FFAG profile needs at least 16 segments and one image")
        edges = np.linspace(0.0, self.cell_length_m, count + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        lengths = np.diff(edges)
        bd = np.zeros(count)
        bf = np.zeros(count)
        for image in range(-images, images + 1):
            shift = image * self.cell_length_m
            (bd_start, bd_stop), (bf_start, bf_stop) = (
                self.magnet_intervals_m)
            bd += _tanh_window(
                centers, bd_start + shift, bd_stop + shift,
                self.fringe_epsilon_m)
            bf += _tanh_window(
                centers, bf_start + shift, bf_stop + shift,
                self.fringe_epsilon_m)
        b0 = (self.defocusing_b0_t * bd
              + self.focusing_b0_t * bf)
        gradient = (self.defocusing_gradient_t_per_m * bd
                    + self.focusing_gradient_t_per_m * bf)
        return (np.ascontiguousarray(centers),
                np.ascontiguousarray(lengths),
                np.ascontiguousarray(b0),
                np.ascontiguousarray(gradient))

    def symmetric_tanh_fringe_integrals(self, *, sample_count=4001):
        """Return the isolated tanh edge diagnostic in the declared gap."""
        count = int(sample_count)
        if count < 101:
            raise ValueError("sample_count must be at least 101")
        extent = 12.0 * self.fringe_epsilon_m
        coordinate = np.linspace(-extent, extent, count)
        normalized = 0.5 * (
            1.0 - np.tanh(coordinate / self.fringe_epsilon_m))
        return enge_fringe_integrals(
            coordinate, normalized, body_field_t=1.0,
            full_gap_m=self.full_gap_m)


def _periodic_planar_orbit(curvature, segment_lengths, rigidity, bend_axis):
    curvature = np.asarray(curvature, dtype=float).reshape(-1)
    lengths = np.asarray(segment_lengths, dtype=float).reshape(-1)
    if curvature.shape != lengths.shape:
        raise ValueError("curvature and segment lengths must match")
    angle = np.r_[0.0, np.cumsum(curvature * lengths)]
    turning = np.diff(angle)
    midpoint_angle = angle[:-1] + 0.5 * turning
    steps = lengths[:, None] * np.column_stack((
        np.cos(midpoint_angle), np.sin(midpoint_angle)))
    relative = np.vstack((np.zeros(2), np.cumsum(steps, axis=0)))
    total_angle = float(angle[-1])
    rotation = np.array([
        [np.cos(total_angle), -np.sin(total_angle)],
        [np.sin(total_angle), np.cos(total_angle)],
    ])
    start = np.linalg.solve(rotation - np.eye(2), relative[-1])
    positions_2d = relative + start
    # All momenta must cross the same radial cell-boundary plane.  Rotate the
    # otherwise arbitrary local solution so its entrance lies on -y.  The
    # entrance tangent is allowed to vary with momentum, as it does in an FFAG.
    entrance_angle = float(np.arctan2(
        positions_2d[0, 1], positions_2d[0, 0]))
    alignment = -0.5 * np.pi - entrance_angle
    align_rotation = np.array([
        [np.cos(alignment), -np.sin(alignment)],
        [np.sin(alignment), np.cos(alignment)],
    ])
    positions_2d = positions_2d @ align_rotation.T
    tangent_2d = np.column_stack((np.cos(angle), np.sin(angle)))
    tangent_2d = tangent_2d @ align_rotation.T
    positions = np.column_stack((positions_2d, np.zeros(len(positions_2d))))
    tangents = np.column_stack((tangent_2d, np.zeros(len(angle))))
    return PlanarDesignOrbit(
        positions, tangents, magnetic_rigidity=float(rigidity),
        bend_axis=np.asarray(bend_axis, dtype=float),
        path_length_stations=np.r_[0.0, np.cumsum(lengths)])


@dataclass(frozen=True)
class FFAGCellReference:
    """One periodic reference orbit and its soft-edge first-order map."""

    kinetic_energy_mev: float
    magnetic_rigidity_tm: float
    transverse_offset_m: float
    orbit: PlanarDesignOrbit
    field_response: np.ndarray
    transfer: CombinedFunctionTransferMap
    bend_angle_rad: float
    periodic_position_residual_m: float
    periodic_tangent_residual: float


def build_ffag_cell_reference(
        spec: FFAGSoftEdgeCellSpec, kinetic_energy_mev, *,
        n_segments=256, response_entries=None) -> FFAGCellReference:
    """Construct the periodic reduced closed orbit for one kinetic energy.

    The Bell--Abell field law is linear in transverse displacement.  The
    single cell-wide displacement is therefore eliminated analytically from
    the total-bend condition.  This is a deterministic soft-edge target
    fixture.  A realized 3-D HDiv field must subsequently recover its own
    closed orbit; the reduced construction is not used as an acceptance
    substitute.
    """
    if not isinstance(spec, FFAGSoftEdgeCellSpec):
        raise TypeError("spec must be an FFAGSoftEdgeCellSpec")
    energy = float(kinetic_energy_mev)
    rigidity = magnetic_rigidity_from_kinetic_energy(energy)
    _, lengths, b0, gradient = spec.sampled_profiles(
        n_segments=n_segments)
    b0_integral = float(b0 @ lengths)
    gradient_integral = float(gradient @ lengths)
    scale = max(1.0, abs(b0_integral))
    if abs(gradient_integral) <= 1.0e-12 * scale:
        raise RuntimeError(
            "cell gradient integral cannot set the reference-orbit offset")
    offset = (
        rigidity * spec.cell_bend_angle_rad - b0_integral
    ) / gradient_integral
    field = b0 + gradient * offset
    curvature = field / rigidity
    orbit = _periodic_planar_orbit(
        curvature, lengths, rigidity, np.array([0.0, 0.0, 1.0]))
    raw = np.r_[field, gradient]
    transfer = combined_function_transfer_map_from_field_response(
        raw, lengths, rigidity, response_entries=response_entries)
    theta = spec.cell_bend_angle_rad
    rotation = np.array([
        [np.cos(theta), -np.sin(theta), 0.0],
        [np.sin(theta), np.cos(theta), 0.0],
        [0.0, 0.0, 1.0],
    ])
    position_residual = float(np.linalg.norm(
        orbit.positions[-1] - rotation @ orbit.positions[0]))
    tangent_residual = float(np.linalg.norm(
        orbit.tangents[-1] - rotation @ orbit.tangents[0]))
    bend_angle = float(np.sum(curvature * lengths))
    return FFAGCellReference(
        energy, rigidity, float(offset), orbit,
        np.ascontiguousarray(raw), transfer, bend_angle,
        position_residual, tangent_residual)


@dataclass(frozen=True)
class FFAGCellTargetFamily:
    """Momentum-indexed reduced cell targets ready for HDiv-MMM fusion."""

    spec: FFAGSoftEdgeCellSpec
    references: tuple[FFAGCellReference, ...]
    objective: MultiMomentumTransferMatrixObjective
    fringe_integrals: EngeFringeIntegrals

    @property
    def kinetic_energies_mev(self) -> np.ndarray:
        return np.asarray([
            reference.kinetic_energy_mev for reference in self.references])


def build_ffag_cell_target_family(
        kinetic_energies_mev, *, spec=None, n_segments=256,
        transfer_matrix_band=1.0e-3, bend_field_band=1.0e-3,
        response_entries=None) -> FFAGCellTargetFamily:
    """Create a multi-momentum cell objective from soft-edge FFAG targets."""
    if spec is None:
        spec = FFAGSoftEdgeCellSpec.bell_abell()
    energies = np.asarray(kinetic_energies_mev, dtype=float).reshape(-1)
    if (energies.size < 2 or not np.all(np.isfinite(energies))
            or np.any(energies <= 0.0)
            or np.any(np.diff(energies) <= 0.0)):
        raise ValueError(
            "FFAG target family needs at least two increasing positive "
            "kinetic energies")
    references = tuple(build_ffag_cell_reference(
        spec, energy, n_segments=n_segments,
        response_entries=response_entries) for energy in energies)
    entries = (references[0].transfer.response_entries
               if response_entries is None else tuple(response_entries))
    objective = MultiMomentumTransferMatrixObjective(
        tuple(reference.orbit for reference in references),
        np.asarray([reference.transfer.matrix for reference in references]),
        transfer_matrix_band, bend_field_band, entries)
    return FFAGCellTargetFamily(
        spec, references, objective,
        spec.symmetric_tanh_fringe_integrals())


def _evaluate_b_field(field, points) -> np.ndarray:
    """Evaluate a callable or ``b_field`` provider in tesla."""
    values = np.asarray(points, dtype=float)
    single = values.shape == (3,)
    points_2d = np.ascontiguousarray(values.reshape(-1, 3))
    if not np.all(np.isfinite(points_2d)):
        raise ValueError("field-evaluation points must be finite")
    evaluator = getattr(field, "b_field", None)
    if evaluator is not None:
        result = np.asarray(evaluator(points_2d), dtype=float)
    else:
        try:
            result = np.asarray(field(points_2d), dtype=float)
        except (TypeError, ValueError):
            result = np.asarray([
                field(float(point[0]), float(point[1]), float(point[2]))
                for point in points_2d], dtype=float)
    if result.shape == (3,) and len(points_2d) == 1:
        result = result[None, :]
    if result.shape != points_2d.shape or not np.all(np.isfinite(result)):
        raise ValueError(
            "magnetic field provider must return one finite 3-vector per "
            "point")
    return result[0] if single else result


def sample_planar_orbit_field_response(
        field, orbit: PlanarDesignOrbit, *, gradient_offset) -> np.ndarray:
    """Sample ``[B_binormal, dB_binormal/dnormal]`` on an orbit.

    The centered normal stencil defines the physical field observable.  It is
    not a design finite difference: whole-element HDiv-MMM sensitivities still
    use the analytic Schur/adjoint contraction of the corresponding rows.
    """
    points, weights = planar_orbit_field_observations(
        orbit, gradient_offset=gradient_offset)
    response = np.einsum(
        "rpc,pc->r", weights, _evaluate_b_field(field, points))
    return np.ascontiguousarray(response, dtype=float)


@dataclass(frozen=True)
class FullFieldClosedOrbit:
    """Periodic planar reference orbit recovered from a realized 3-D field."""

    magnetic_rigidity_tm: float
    orbit: PlanarDesignOrbit
    path_length_m: float
    entrance_radius_m: float
    entrance_incidence_angle_rad: float
    periodic_position_residual_m: float
    periodic_tangent_residual: float
    vertical_position_residual_m: float
    vertical_tangent_residual: float
    root_evaluations: int
    field_response: np.ndarray
    transfer: CombinedFunctionTransferMap

    @property
    def closure_residual(self) -> float:
        return float(max(
            self.periodic_position_residual_m,
            self.periodic_tangent_residual,
            abs(self.vertical_position_residual_m),
            abs(self.vertical_tangent_residual)))


def recover_periodic_planar_closed_orbit(
        field, *, magnetic_rigidity, cell_angle_rad, initial_radius_m,
        initial_incidence_angle_rad=0.0, n_segments=128,
        gradient_offset=1.0e-3, max_path_length_m=None,
        curvature_sign=1.0, position_tolerance=1.0e-9,
        tangent_tolerance=1.0e-9, root_max_evaluations=80,
        response_entries=None) -> FullFieldClosedOrbit:
    """Recover one-cell periodic orbit and its local transfer map.

    Two geometric unknowns are solved: the entrance radius and incidence
    angle on the radial cell boundary.  A trajectory is integrated through
    the supplied total magnetic field until it reaches the next radial cell
    boundary.  The exit point and tangent, rotated back by one cell angle,
    must equal the entrance data.  Numerical differencing used internally by
    the two-variable root finder concerns orbit recovery only; no design
    derivative or topology sensitivity is approximated.

    ``curvature_sign=+1`` matches :class:`PlanarTransferMatrixObjective`: a
    positive binormal field produces positive signed curvature.  It therefore
    fixes the otherwise conventional charge/field orientation in the Lorentz
    equation used here.
    """
    from scipy.integrate import solve_ivp
    from scipy.optimize import least_squares

    rigidity = float(magnetic_rigidity)
    theta = float(cell_angle_rad)
    radius_guess = float(initial_radius_m)
    alpha_guess = float(initial_incidence_angle_rad)
    segments = int(n_segments)
    gradient_offset = float(gradient_offset)
    curvature_sign = float(curvature_sign)
    position_tolerance = float(position_tolerance)
    tangent_tolerance = float(tangent_tolerance)
    max_evaluations = int(root_max_evaluations)
    if (not np.all(np.isfinite([
            rigidity, theta, radius_guess, alpha_guess, gradient_offset,
            curvature_sign, position_tolerance, tangent_tolerance]))
            or rigidity <= 0.0 or radius_guess <= 0.0
            or theta <= 0.0 or theta >= np.pi
            or segments < 8 or gradient_offset <= 0.0
            or curvature_sign == 0.0 or position_tolerance <= 0.0
            or tangent_tolerance <= 0.0 or max_evaluations < 4):
        raise ValueError("invalid periodic closed-orbit recovery settings")
    if max_path_length_m is None:
        max_path = max(
            4.0 * radius_guess * theta, 0.25 * radius_guess)
    else:
        max_path = float(max_path_length_m)
    if not np.isfinite(max_path) or max_path <= 0.0:
        raise ValueError("max_path_length_m must be positive")

    axis = np.array([0.0, 0.0, 1.0])
    radial_0 = np.array([0.0, -1.0, 0.0])
    tangent_0 = np.array([1.0, 0.0, 0.0])
    cosine = float(np.cos(theta))
    sine = float(np.sin(theta))
    rotation = np.array([
        [cosine, -sine, 0.0],
        [sine, cosine, 0.0],
        [0.0, 0.0, 1.0],
    ])
    radial_1 = rotation @ radial_0
    tangent_1 = rotation @ tangent_0
    radius_scale = max(radius_guess, 1.0e-6)

    def integrate(radius, alpha, *, dense_output=False):
        start_position = radius * radial_0
        start_tangent = (
            np.cos(alpha) * tangent_0 + np.sin(alpha) * radial_0)
        state_0 = np.r_[start_position, start_tangent]

        def ode(_path_length, state):
            tangent = state[3:]
            tangent_norm = float(np.linalg.norm(tangent))
            if tangent_norm <= 0.0:
                raise RuntimeError("particle tangent vanished")
            tangent = tangent / tangent_norm
            magnetic_field = _evaluate_b_field(field, state[:3])
            curvature = (-curvature_sign
                         * np.cross(tangent, magnetic_field) / rigidity)
            return np.r_[tangent, curvature]

        def next_boundary(_path_length, state):
            return float(state[:3] @ tangent_1)

        next_boundary.direction = 1.0
        next_boundary.terminal = True
        return solve_ivp(
            ode, (0.0, max_path), state_0, method="DOP853",
            events=next_boundary, dense_output=dense_output,
            rtol=min(1.0e-10, 0.1 * tangent_tolerance),
            atol=min(1.0e-12, 0.1 * position_tolerance),
            max_step=max_path / max(segments, 32))

    def terminal_state(radius, alpha, *, dense_output=False):
        solution = integrate(radius, alpha, dense_output=dense_output)
        if (not solution.success or len(solution.t_events) != 1
                or len(solution.t_events[0]) != 1):
            return None, solution
        terminal = np.asarray(solution.y_events[0][0], dtype=float)
        if terminal[:3] @ radial_1 <= 0.0:
            return None, solution
        terminal[3:] /= np.linalg.norm(terminal[3:])
        return terminal, solution

    def residual(parameters):
        radius, alpha = parameters
        terminal, _ = terminal_state(radius, alpha)
        if terminal is None:
            return np.array([1.0e3, 1.0e3])
        position_back = rotation.T @ terminal[:3]
        tangent_back = rotation.T @ terminal[3:]
        entrance_tangent = (
            np.cos(alpha) * tangent_0 + np.sin(alpha) * radial_0)
        tangent_cross = float(axis @ np.cross(
            entrance_tangent, tangent_back))
        tangent_dot = float(entrance_tangent @ tangent_back)
        return np.array([
            (position_back @ radial_0 - radius) / radius_scale,
            np.arctan2(tangent_cross, tangent_dot),
        ])

    root = least_squares(
        residual, np.array([radius_guess, alpha_guess]),
        bounds=(np.array([0.2 * radius_guess, -0.45 * np.pi]),
                np.array([5.0 * radius_guess, 0.45 * np.pi])),
        xtol=min(1.0e-12, position_tolerance / radius_scale),
        ftol=min(1.0e-12, tangent_tolerance),
        gtol=min(1.0e-12, tangent_tolerance),
        max_nfev=max_evaluations)
    radius, alpha = (float(value) for value in root.x)
    terminal, solution = terminal_state(radius, alpha, dense_output=True)
    if terminal is None:
        raise RuntimeError(
            "periodic orbit did not reach the next radial boundary")
    final_residual = residual((radius, alpha))
    position_residual = abs(float(final_residual[0])) * radius_scale
    tangent_residual = abs(float(final_residual[1]))
    position_back = rotation.T @ terminal[:3]
    tangent_back = rotation.T @ terminal[3:]
    if (not root.success or position_residual > position_tolerance
            or tangent_residual > tangent_tolerance):
        raise RuntimeError(
            "periodic orbit closure failed: position residual "
            f"{position_residual:.6e} m, tangent residual "
            f"{tangent_residual:.6e}")
    path_length = float(solution.t_events[0][0])
    path_stations = np.linspace(0.0, path_length, segments + 1)
    states = np.asarray(solution.sol(path_stations), dtype=float).T
    tangents = states[:, 3:]
    tangents /= np.linalg.norm(tangents, axis=1)[:, None]
    orbit = PlanarDesignOrbit(
        states[:, :3], tangents, magnetic_rigidity=rigidity,
        bend_axis=axis, path_length_stations=path_stations)
    response = sample_planar_orbit_field_response(
        field, orbit, gradient_offset=gradient_offset)
    transfer = combined_function_transfer_map_from_field_response(
        response, orbit.segment_lengths, rigidity,
        curvature_sign=curvature_sign,
        response_entries=response_entries)
    return FullFieldClosedOrbit(
        rigidity, orbit, path_length, radius, alpha,
        position_residual, tangent_residual,
        float(position_back[2]), float(tangent_back[2]), int(root.nfev),
        response, transfer)


def recover_ffag_closed_orbit(
        field, spec: FFAGSoftEdgeCellSpec, kinetic_energy_mev, *,
        initial_reference=None, n_segments=128, gradient_offset=1.0e-3,
        **recovery_options) -> FullFieldClosedOrbit:
    """Recover one realized periodic FFAG orbit using a reduced target seed."""
    if not isinstance(spec, FFAGSoftEdgeCellSpec):
        raise TypeError("spec must be an FFAGSoftEdgeCellSpec")
    energy = float(kinetic_energy_mev)
    reference = (build_ffag_cell_reference(
        spec, energy, n_segments=max(64, int(n_segments)))
        if initial_reference is None else initial_reference)
    if isinstance(reference, FFAGCellReference):
        if (abs(reference.kinetic_energy_mev - energy)
                > 1.0e-12 * max(1.0, energy)):
            raise ValueError("initial_reference kinetic energy does not match")
        seed_orbit = reference.orbit
        rigidity = reference.magnetic_rigidity_tm
        radius = None
        alpha = None
    elif isinstance(reference, FullFieldClosedOrbit):
        rigidity = magnetic_rigidity_from_kinetic_energy(energy)
        if (abs(reference.magnetic_rigidity_tm - rigidity)
                > 1.0e-10 * max(1.0, rigidity)):
            raise ValueError("initial full-field orbit rigidity does not match")
        seed_orbit = reference.orbit
        radius = reference.entrance_radius_m
        alpha = reference.entrance_incidence_angle_rad
    else:
        raise TypeError(
            "initial_reference must be an FFAGCellReference or "
            "FullFieldClosedOrbit")
    if radius is None:
        radial_0 = np.array([0.0, -1.0, 0.0])
        tangent_0 = np.array([1.0, 0.0, 0.0])
        radius = float(seed_orbit.positions[0] @ radial_0)
        entrance_tangent = seed_orbit.tangents[0]
        alpha = float(np.arctan2(
            entrance_tangent @ radial_0, entrance_tangent @ tangent_0))
    recovery_options.setdefault(
        "max_path_length_m", 2.5 * spec.cell_length_m)
    return recover_periodic_planar_closed_orbit(
        field, magnetic_rigidity=rigidity,
        cell_angle_rad=spec.cell_bend_angle_rad,
        initial_radius_m=radius,
        initial_incidence_angle_rad=alpha,
        n_segments=n_segments, gradient_offset=gradient_offset,
        **recovery_options)


def recover_ffag_closed_orbit_family(
        field, target_family: FFAGCellTargetFamily, *, n_segments=128,
        gradient_offset=1.0e-3, initial_references=None,
        **recovery_options
        ) -> tuple[FullFieldClosedOrbit, ...]:
    """Recover all momentum orbits from one shared realized magnet field."""
    if not isinstance(target_family, FFAGCellTargetFamily):
        raise TypeError("target_family must be an FFAGCellTargetFamily")
    references = (target_family.references if initial_references is None
                  else tuple(initial_references))
    if len(references) != len(target_family.references):
        raise ValueError(
            "initial_references must have one seed per target momentum")
    return tuple(recover_ffag_closed_orbit(
        field, target_family.spec, target.kinetic_energy_mev,
        initial_reference=seed, n_segments=n_segments,
        gradient_offset=gradient_offset, **recovery_options)
        for target, seed in zip(target_family.references, references))


def _realized_ffag_band_ratios(
        recovered_orbits, objective: MultiMomentumTransferMatrixObjective
        ) -> tuple[np.ndarray, np.ndarray]:
    field_ratios = []
    matrix_ratios = []
    for recovered, target in zip(recovered_orbits, objective.objectives):
        count = len(recovered.orbit.segment_lengths)
        if count != len(target.bend_field_band):
            raise ValueError(
                "recovered and target orbit segment counts must match")
        required_field = (
            recovered.orbit.magnetic_rigidity
            * recovered.orbit.signed_curvature / target.curvature_sign)
        field_ratios.append(float(np.max(np.abs(
            (recovered.field_response[:count] - required_field)
            / target.bend_field_band))))
        realized_values = np.asarray([
            recovered.transfer.matrix[row, column]
            for row, column in target.response_entries])
        target_values = np.asarray([
            target.target_matrix[row, column]
            for row, column in target.response_entries])
        bands = np.asarray([
            target.transfer_matrix_band[row, column]
            for row, column in target.response_entries])
        matrix_ratios.append(float(np.max(np.abs(
            (realized_values - target_values) / bands))))
    return np.asarray(field_ratios), np.asarray(matrix_ratios)


@dataclass(frozen=True)
class FFAGHDivMMMOuterIteration:
    """One Lego update followed by exact full-field orbit reconstruction."""

    index: int
    material_move_fraction: float
    active_count_before: int
    active_count_after: int
    source_scale_before: float
    source_scale_after: float
    max_band_ratio_before: float
    max_band_ratio_after: float
    max_position_closure_residual_m: float
    max_tangent_closure_residual: float
    accepted: bool
    reason: str
    topology_result: MultiMomentumAcceleratorMagnetTopologyResult


@dataclass(frozen=True)
class FFAGHDivMMMTopologyResult:
    """Full-field closed-orbit outer loop around binary HDiv-MMM growth."""

    target_family: FFAGCellTargetFamily
    source_scale: float
    active_elements: np.ndarray
    state: np.ndarray
    recovered_orbits: tuple[FullFieldClosedOrbit, ...]
    orbit_field_max_band_ratios: np.ndarray
    transfer_matrix_max_band_ratios: np.ndarray
    history: tuple[FFAGHDivMMMOuterIteration, ...]
    converged: bool
    stop_reason: str
    topology: object

    @property
    def realized_transfer_matrices(self) -> np.ndarray:
        return np.asarray([
            recovered.transfer.matrix for recovered in self.recovered_orbits])

    @property
    def max_band_ratio(self) -> float:
        return float(max(
            np.max(self.orbit_field_max_band_ratios),
            np.max(self.transfer_matrix_max_band_ratios)))


@dataclass(frozen=True)
class FFAGFixedOrbitMapTrial:
    """One exactly resolved material proposal checked against the one-pass map."""

    optics_iteration: int
    trial_index: int
    material_move_fraction: float | None
    active_count_before: int
    active_count_after: int
    max_band_ratio_before: float
    max_band_ratio_after: float
    accepted: bool
    reason: str
    proposal_model: str = "field-target"
    exact_search_trace: tuple = ()


@dataclass(frozen=True)
class FFAGFixedOrbitHDivMMMTopologyResult:
    """One-pass FFAG result about caller-supplied design orbits.

    Unlike :class:`FFAGHDivMMMTopologyResult`, this contract performs no
    periodic closed-orbit search.  The momentum-indexed design orbits in the
    target family are frozen observation paths from entrance to exit.  Their
    transfer matrices are therefore the one-pass maps that the material
    inverse must reproduce.
    """

    target_family: FFAGCellTargetFamily
    source_scale: float
    topology_result: MultiMomentumAcceleratorMagnetTopologyResult
    optics_history: tuple[MultiMomentumAcceleratorMagnetTopologyResult, ...]
    termination_reason: str
    initial_max_band_ratio: float
    map_trust_history: tuple[FFAGFixedOrbitMapTrial, ...]

    @property
    def active_elements(self) -> np.ndarray:
        return self.topology_result.active_elements

    @property
    def state(self) -> np.ndarray:
        return self.topology_result.generation.state

    @property
    def realized_transfer_matrices(self) -> np.ndarray:
        return self.topology_result.realized_transfer_matrices

    @property
    def orbit_field_max_band_ratios(self) -> np.ndarray:
        return self.topology_result.orbit_field_max_band_ratios

    @property
    def transfer_matrix_max_band_ratios(self) -> np.ndarray:
        return self.topology_result.transfer_matrix_max_band_ratios

    @property
    def topology(self):
        return self.topology_result.topology

    @property
    def field_correction(self):
        return self.topology_result.field_correction

    @property
    def history(self):
        return tuple(
            item
            for result in self.optics_history
            for item in result.generation.history)

    @property
    def converged(self) -> bool:
        return self.max_band_ratio <= 1.0

    @property
    def stop_reason(self) -> str:
        return self.termination_reason

    @property
    def max_band_ratio(self) -> float:
        return float(max(
            np.max(self.orbit_field_max_band_ratios),
            np.max(self.transfer_matrix_max_band_ratios)))


def optimize_ffag_hdiv_mmm_from_fixed_design_orbits(
        target_family: FFAGCellTargetFamily, *, source,
        charge_gram, fes, inv_chi, active_elements, element_volumes,
        volume_max, gradient_offset=1.0e-3, source_scale=1.0,
        optimize_source_scale=True,
        max_optics_iterations=1, material_iterations_per_optics=1,
        field_inverse_relative_tolerance=1.0e-3,
        field_inverse_basis=None,
        field_inverse_maximum_step_scale=1.0,
        field_inverse_line_search_steps=8,
        map_trust_region_trials=3, map_ratio_tolerance=1.0e-8,
        direct_map_oracle_fallback=False,
        direct_map_oracle_exact_beam_width=0,
        direct_map_oracle_exact_beam_depth=0,
        direct_map_oracle_graph_front_proposal_limit=0,
        **generation_options) -> FFAGFixedOrbitHDivMMMTopologyResult:
    """Optimize one magnet about fixed entrance-to-exit design orbits.

    This is the production proof-of-concept path when the design orbit and
    desired one-pass transfer matrix are inputs.  It deliberately omits ring
    closure and periodic-orbit recovery.  The optics TSVD first maps the
    transfer-matrix error to a sampled field target on those frozen paths;
    the independent Abe--Murata ACA--QR--TSVD material inverse then selects
    binary HDiv-MMM element additions/removals.  The field-to-map Jacobian is
    propagated by forward-mode AD with an exact Frechet matrix-exponential
    primitive.

    The caller owns ``ngsolve.TaskManager``.  No design finite difference,
    density interpolation, air volume mesh, or tracking root solve is used.
    """
    from .topology_optimization import solve_hdiv_mmm_active_elements

    if not isinstance(target_family, FFAGCellTargetFamily):
        raise TypeError("target_family must be an FFAGCellTargetFamily")
    if not isinstance(source, CoilBuilderHDivSource):
        raise TypeError("source must be a CoilBuilderHDivSource")
    active = np.asarray(active_elements, dtype=bool).reshape(-1).copy()
    volumes = np.asarray(element_volumes, dtype=float).reshape(-1)
    scale = float(source_scale)
    optics_count = int(max_optics_iterations)
    material_count = int(material_iterations_per_optics)
    map_trial_count = int(map_trust_region_trials)
    map_ratio_tolerance = float(map_ratio_tolerance)
    oracle_beam_width = int(direct_map_oracle_exact_beam_width)
    oracle_beam_depth = int(direct_map_oracle_exact_beam_depth)
    oracle_graph_limit = int(
        direct_map_oracle_graph_front_proposal_limit)
    if (active.shape != volumes.shape or not np.any(active)
            or not np.all(np.isfinite(volumes)) or np.any(volumes <= 0.0)
            or not np.isfinite(scale) or scale <= 0.0):
        raise ValueError("invalid fixed-design-orbit topology settings")
    if optics_count < 1 or material_count < 1 or map_trial_count < 1:
        raise ValueError(
            "fixed-design-orbit iteration counts must be positive")
    if not np.isfinite(map_ratio_tolerance) or map_ratio_tolerance < 0.0:
        raise ValueError(
            "map_ratio_tolerance must be nonnegative and finite")
    if (oracle_beam_width < 0 or oracle_beam_depth < 0
            or ((oracle_beam_width == 0)
                != (oracle_beam_depth == 0))):
        raise ValueError(
            "direct map oracle beam width and depth must both be zero or "
            "both be positive")
    if oracle_graph_limit < 0:
        raise ValueError(
            "direct map oracle graph-front proposal limit must be "
            "nonnegative")
    if "max_iterations" in generation_options:
        raise TypeError(
            "use material_iterations_per_optics instead of max_iterations")

    objective = target_family.objective
    source_rhs = source.assemble_hdiv_rhs(fes)
    rhs = scale * source_rhs
    state = solve_hdiv_mmm_active_elements(
        charge_gram=charge_gram, fes=fes, inv_chi=inv_chi, rhs=rhs,
        response_matrix=np.zeros((1, int(fes.ndof))),
        active_elements=active,
        solve_tolerance=float(generation_options.get(
            "solve_tolerance", 1.0e-9)),
        solve_max_iterations=int(generation_options.get(
            "solve_max_iterations", 5000)),
        mass_riesz=bool(generation_options.get("mass_riesz", True)),
        cluster_coarse_size=int(generation_options.get(
            "cluster_coarse_size", 0)),
        cluster_deflation_size=int(generation_options.get(
            "cluster_deflation_size", 0)),
        recycle_size=int(generation_options.get("recycle_size", 0)))[0]

    if optimize_source_scale:
        total_field = CoilHDivTotalField(
            source, charge_gram, state, source_scale=scale)
        sampled = []
        required = []
        bands = []
        for reference, target in zip(
                target_family.references, objective.objectives):
            response = sample_planar_orbit_field_response(
                total_field, reference.orbit,
                gradient_offset=gradient_offset)
            count = len(reference.orbit.segment_lengths)
            sampled.extend(response[:count])
            required.extend(target.required_bend_field)
            bands.extend(target.bend_field_band)
        sampled = np.asarray(sampled, dtype=float)
        required = np.asarray(required, dtype=float)
        weights = 1.0 / np.asarray(bands, dtype=float)
        denominator = float((weights * sampled) @ (weights * sampled))
        calibration = float(
            (weights * sampled) @ (weights * required) / denominator)
        if (not np.isfinite(calibration) or calibration <= 0.0
                or denominator <= np.finfo(float).tiny):
            raise RuntimeError(
                "CoilBuilder source cannot be positively calibrated to the "
                "fixed design-orbit bend fields")
        scale *= calibration
        rhs = scale * source_rhs
        state *= calibration

    response_matrix = build_multi_orbit_field_response_matrix(
        charge_gram, objective, gradient_offset=gradient_offset)
    incident = scale * source.incident_orbit_field_response(
        objective, gradient_offset=gradient_offset)
    optics_history = []
    map_trust_history = []
    accepted_result = None
    initial_max_band_ratio = None
    termination_reason = "maximum fixed-orbit optics iterations reached"
    for optics_iteration in range(optics_count):
        current_raw_field = np.asarray(
            response_matrix @ state + incident, dtype=float)
        current_ratio = float(np.max(np.abs(
            (objective.transform(current_raw_field)
             - objective.response_target) / objective.response_band)))
        if initial_max_band_ratio is None:
            initial_max_band_ratio = current_ratio
        field_correction = solve_transfer_matrix_field_correction(
            objective, current_raw_field,
            field_basis=field_inverse_basis,
            relative_tolerance=field_inverse_relative_tolerance,
            maximum_step_scale=field_inverse_maximum_step_scale,
            line_search_steps=field_inverse_line_search_steps)
        requested_initial = generation_options.get(
            "initial_material_move_fraction")
        requested_maximum = generation_options.get(
            "maximum_material_move_fraction")
        trial_fraction = (None if requested_initial is None else
                          float(requested_initial))
        accepted = False
        last_attempt = None
        for trial_index in range(map_trial_count):
            trial_options = dict(generation_options)
            if trial_fraction is not None:
                trial_options["initial_material_move_fraction"] = trial_fraction
                trial_options["maximum_material_move_fraction"] = min(
                    trial_fraction,
                    trial_fraction if requested_maximum is None else
                    float(requested_maximum))
            last_attempt = optimize_hdiv_mmm_magnet_from_transfer_matrices(
                objective.orbits, objective.target_matrices,
                transfer_matrix_band=objective.transfer_matrix_band,
                bend_field_band=objective.bend_field_band,
                charge_gram=charge_gram, fes=fes, inv_chi=inv_chi,
                rhs=rhs, field_response_matrix=response_matrix,
                incident_field_response=incident,
                field_correction=field_correction,
                active_elements=active, element_volumes=volumes,
                volume_max=volume_max,
                response_entries=objective.response_entries,
                curvature_sign=objective.curvature_sign,
                gradient_sign=objective.gradient_sign,
                max_iterations=material_count,
                **trial_options)
            next_active = np.asarray(
                last_attempt.active_elements, dtype=bool).copy()
            candidate_ratio = float(max(
                np.max(last_attempt.orbit_field_max_band_ratios),
                np.max(last_attempt.transfer_matrix_max_band_ratios)))
            changed = not np.array_equal(next_active, active)
            accepted = bool(
                candidate_ratio <= 1.0 + map_ratio_tolerance or
                (changed and candidate_ratio
                 < current_ratio - map_ratio_tolerance))
            if accepted:
                reason = "accepted by exact fixed one-pass map gate"
            elif not changed:
                reason = "material inverse proposed no active-set change"
            else:
                reason = "rejected by exact fixed one-pass map gate"
            map_trust_history.append(FFAGFixedOrbitMapTrial(
                optics_iteration, trial_index, trial_fraction,
                int(np.count_nonzero(active)),
                int(np.count_nonzero(next_active)), current_ratio,
                candidate_ratio, accepted, reason, "field-target",
                tuple(last_attempt.generation.exact_search_trace)))
            if accepted:
                accepted_result = last_attempt
                optics_history.append(last_attempt)
                active = next_active
                state = last_attempt.generation.state.copy()
                break
            if trial_fraction is None:
                break
            trial_fraction *= 0.5
        if not accepted and direct_map_oracle_fallback:
            # The two-stage transfer->field->material inverse can stall when
            # its reachable field target is a poor local surrogate for the
            # original map norm.  Reuse the same forward-mode AD field-to-map
            # Jacobian directly in the all-candidate material contraction as
            # a bounded fallback.  This is still an exact chain rule, not a
            # design finite difference or a density relaxation.
            oracle_options = dict(generation_options)
            oracle_fraction = (None if requested_initial is None else
                               float(requested_initial))
            if oracle_fraction is not None:
                oracle_options["initial_material_move_fraction"] = (
                    oracle_fraction)
                oracle_options["maximum_material_move_fraction"] = min(
                    oracle_fraction,
                    oracle_fraction if requested_maximum is None else
                    float(requested_maximum))
            # The default remains one bounded global all-candidate proposal.
            # A caller may explicitly enable shallow nonmonotone look-ahead
            # for the direct oracle after the primary field-target lane stalls.
            oracle_options["exact_beam_width"] = oracle_beam_width
            oracle_options["exact_beam_depth"] = oracle_beam_depth
            oracle_options["graph_front_proposal_limit"] = (
                oracle_graph_limit)
            last_attempt = optimize_hdiv_mmm_magnet_from_transfer_matrices(
                objective.orbits, objective.target_matrices,
                transfer_matrix_band=objective.transfer_matrix_band,
                bend_field_band=objective.bend_field_band,
                charge_gram=charge_gram, fes=fes, inv_chi=inv_chi,
                rhs=rhs, field_response_matrix=response_matrix,
                incident_field_response=incident,
                active_elements=active, element_volumes=volumes,
                volume_max=volume_max,
                response_entries=objective.response_entries,
                curvature_sign=objective.curvature_sign,
                gradient_sign=objective.gradient_sign,
                max_iterations=material_count,
                **oracle_options)
            next_active = np.asarray(
                last_attempt.active_elements, dtype=bool).copy()
            candidate_ratio = float(max(
                np.max(last_attempt.orbit_field_max_band_ratios),
                np.max(last_attempt.transfer_matrix_max_band_ratios)))
            changed = not np.array_equal(next_active, active)
            accepted = bool(
                candidate_ratio <= 1.0 + map_ratio_tolerance or
                (changed and candidate_ratio
                 < current_ratio - map_ratio_tolerance))
            reason = (
                "accepted by direct analytic map-Jacobian oracle"
                if accepted else
                ("direct map oracle proposed no active-set change"
                 if not changed else
                 "rejected by exact fixed one-pass map gate"))
            map_trust_history.append(FFAGFixedOrbitMapTrial(
                optics_iteration, map_trial_count, oracle_fraction,
                int(np.count_nonzero(active)),
                int(np.count_nonzero(next_active)), current_ratio,
                candidate_ratio, accepted, reason,
                "direct-map-jacobian",
                tuple(last_attempt.generation.exact_search_trace)))
            if accepted:
                accepted_result = last_attempt
                optics_history.append(last_attempt)
                active = next_active
                state = last_attempt.generation.state.copy()
        if not accepted:
            termination_reason = "map-level trust-region proposals rejected"
            break
        if candidate_ratio <= 1.0 + map_ratio_tolerance:
            termination_reason = "fixed one-pass transfer bands reached"
            break

    if accepted_result is None:
        # Preserve the incumbent when every material proposal is rejected.
        # This fallback performs one exact active solve only on that failure
        # path; it never returns the last rejected topology as the design.
        baseline_options = dict(generation_options)
        accepted_result = optimize_hdiv_mmm_magnet_from_transfer_matrices(
            objective.orbits, objective.target_matrices,
            transfer_matrix_band=objective.transfer_matrix_band,
            bend_field_band=objective.bend_field_band,
            charge_gram=charge_gram, fes=fes, inv_chi=inv_chi,
            rhs=rhs, field_response_matrix=response_matrix,
            incident_field_response=incident,
            field_correction=field_correction,
            active_elements=active, element_volumes=volumes,
            volume_max=volume_max,
            response_entries=objective.response_entries,
            curvature_sign=objective.curvature_sign,
            gradient_sign=objective.gradient_sign,
            max_iterations=0, **baseline_options)
    return FFAGFixedOrbitHDivMMMTopologyResult(
        target_family, scale, accepted_result, tuple(optics_history),
        termination_reason, float(initial_max_band_ratio),
        tuple(map_trust_history))


def optimize_ffag_hdiv_mmm_from_transfer_matrices(
        target_family: FFAGCellTargetFamily, *, source,
        charge_gram, fes, inv_chi, active_elements, element_volumes,
        volume_max, gradient_offset=1.0e-3, source_scale=1.0,
        optimize_source_scale=True,
        hdiv_order=1, orbit_segments=None, max_outer_iterations=8,
        inner_iterations=1, outer_initial_material_move_fraction=0.10,
        outer_trust_region_trials=3, outer_ratio_tolerance=1.0e-8,
        field_inverse_relative_tolerance=1.0e-3,
        field_inverse_basis=None,
        field_inverse_maximum_step_scale=1.0,
        field_inverse_line_search_steps=8,
        recovery_options=None, **generation_options
        ) -> FFAGHDivMMMTopologyResult:
    """Optimize a binary FFAG magnet with orbit recovery after every move.

    The fixed CoilBuilder source is assembled once into the HDiv RHS.  At each
    outer iteration the realized coil-plus-magnet field is tracked to recover
    every momentum-dependent periodic orbit.  A small dense optics TSVD first
    converts transfer-matrix error to a target correction of the sampled orbit
    field using the forward-mode AD field-to-map Jacobian.  Native HDiv rows
    are then rebuilt on those orbits and the separate Abe--Murata DUCAS
    ACA--QR--TSVD material inverse proposes exactly one (by default) Lego
    update.  HDiv candidates never enter the optics inverse.
    The topology is accepted only if a complete
    active-set solve, full-field orbit recovery, and transfer-map rebuild
    improve the actual band-normalized objective.  A rejected update shrinks
    the whole-element material trust region; it is never passed to Trafo.

    The caller owns ``ngsolve.TaskManager`` for the whole call.  No design
    finite difference, density interpolation, or air volume mesh is used.
    """
    from .topology_optimization import (
        ngsolve_growth_topology,
        solve_hdiv_mmm_active_elements,
    )

    if not isinstance(target_family, FFAGCellTargetFamily):
        raise TypeError("target_family must be an FFAGCellTargetFamily")
    if not isinstance(source, CoilBuilderHDivSource):
        raise TypeError("source must be a CoilBuilderHDivSource")
    active = np.asarray(active_elements, dtype=bool).reshape(-1).copy()
    volumes = np.asarray(element_volumes, dtype=float).reshape(-1)
    source_scale = float(source_scale)
    optimize_source_scale = bool(optimize_source_scale)
    outer_count = int(max_outer_iterations)
    inner_count = int(inner_iterations)
    trial_count = int(outer_trust_region_trials)
    move_fraction = float(outer_initial_material_move_fraction)
    ratio_tolerance = float(outer_ratio_tolerance)
    if (active.shape != volumes.shape or not np.any(active)
            or not np.all(np.isfinite(volumes)) or np.any(volumes <= 0.0)
            or not np.isfinite(source_scale) or source_scale <= 0.0
            or outer_count < 1 or inner_count < 1 or trial_count < 1
            or not np.isfinite(move_fraction)
            or move_fraction <= 0.0 or move_fraction > 1.0
            or not np.isfinite(ratio_tolerance) or ratio_tolerance < 0.0):
        raise ValueError("invalid FFAG HDiv-MMM outer-loop settings")
    reserved = {
        "max_iterations", "initial_material_move_fraction",
        "maximum_material_move_fraction", "source_calibration_rows",
        "source_calibration_target"}
    overlap = reserved.intersection(generation_options)
    if overlap:
        raise TypeError(
            "outer loop owns " + ", ".join(sorted(overlap)))
    segment_count = (len(target_family.references[0].orbit.segment_lengths)
                     if orbit_segments is None else int(orbit_segments))
    if (segment_count < 8 or any(
            len(reference.orbit.segment_lengths) != segment_count
            for reference in target_family.references)):
        raise ValueError(
            "orbit_segments must match every target-family orbit")
    recovery = {} if recovery_options is None else dict(recovery_options)
    recovery.setdefault("response_entries", target_family.objective.response_entries)

    source_rhs = source.assemble_hdiv_rhs(fes)
    rhs = source_scale * source_rhs
    zero_response = np.zeros((1, int(fes.ndof)))
    state = solve_hdiv_mmm_active_elements(
        charge_gram=charge_gram, fes=fes, inv_chi=inv_chi, rhs=rhs,
        response_matrix=zero_response, active_elements=active,
        solve_tolerance=float(generation_options.get(
            "solve_tolerance", 1.0e-9)),
        solve_max_iterations=int(generation_options.get(
            "solve_max_iterations", 5000)),
        mass_riesz=bool(generation_options.get("mass_riesz", True)),
        cluster_coarse_size=int(generation_options.get(
            "cluster_coarse_size", 0)),
        cluster_deflation_size=int(generation_options.get(
            "cluster_deflation_size", 0)),
        recycle_size=int(generation_options.get("recycle_size", 0)))[0]

    total_field = CoilHDivTotalField(
        source, charge_gram, state, source_scale=source_scale,
        hdiv_order=hdiv_order)
    if optimize_source_scale:
        sampled = []
        required = []
        bands = []
        for reference, target in zip(
                target_family.references,
                target_family.objective.objectives):
            response = sample_planar_orbit_field_response(
                total_field, reference.orbit,
                gradient_offset=gradient_offset)
            count = len(reference.orbit.segment_lengths)
            sampled.extend(response[:count])
            required.extend(target.required_bend_field)
            bands.extend(target.bend_field_band)
        sampled = np.asarray(sampled)
        required = np.asarray(required)
        weights = 1.0 / np.asarray(bands)
        denominator = float((weights * sampled) @ (weights * sampled))
        calibration = float(
            (weights * sampled) @ (weights * required) / denominator)
        if (not np.isfinite(calibration) or calibration <= 0.0
                or denominator <= np.finfo(float).tiny):
            raise RuntimeError(
                "CoilBuilder source cannot be positively calibrated to the "
                "target bend fields")
        source_scale *= calibration
        rhs = source_scale * source_rhs
        state *= calibration
        total_field = CoilHDivTotalField(
            source, charge_gram, state, source_scale=source_scale,
            hdiv_order=hdiv_order)
    recovered = recover_ffag_closed_orbit_family(
        total_field, target_family, n_segments=segment_count,
        gradient_offset=gradient_offset, **recovery)
    field_ratios, matrix_ratios = _realized_ffag_band_ratios(
        recovered, target_family.objective)
    current_ratio = float(max(np.max(field_ratios), np.max(matrix_ratios)))
    history = []
    stop_reason = "maximum outer iterations reached"

    for outer_index in range(outer_count):
        if current_ratio <= 1.0 + ratio_tolerance:
            stop_reason = "full-field orbit and transfer bands reached"
            break
        dynamic_objective = MultiMomentumTransferMatrixObjective(
            tuple(value.orbit for value in recovered),
            target_family.objective.target_matrices,
            target_family.objective.transfer_matrix_band,
            target_family.objective.bend_field_band,
            target_family.objective.response_entries,
            target_family.objective.curvature_sign,
            target_family.objective.gradient_sign)
        response_matrix = build_multi_orbit_field_response_matrix(
            charge_gram, dynamic_objective,
            gradient_offset=gradient_offset)
        incident = source_scale * source.incident_orbit_field_response(
            dynamic_objective, gradient_offset=gradient_offset)
        current_raw_field = np.asarray(
            response_matrix @ state + incident, dtype=float)
        field_correction = solve_transfer_matrix_field_correction(
            dynamic_objective, current_raw_field,
            field_basis=field_inverse_basis,
            relative_tolerance=field_inverse_relative_tolerance,
            maximum_step_scale=field_inverse_maximum_step_scale,
            line_search_steps=field_inverse_line_search_steps)
        calibration_rows = []
        calibration_target = []
        if optimize_source_scale:
            offsets = dynamic_objective.raw_offsets
            for index, objective in enumerate(dynamic_objective.objectives):
                count = len(objective.orbit.segment_lengths)
                calibration_rows.extend(range(
                    int(offsets[index]), int(offsets[index]) + count))
                calibration_target.extend(objective.required_bend_field)
            calibration_rows = np.asarray(calibration_rows, dtype=np.int64)
            calibration_target = np.asarray(calibration_target, dtype=float)
        accepted = False
        attempted_result = None
        trial_fraction = move_fraction
        trial_reason = "no improving full-field update"
        for _ in range(trial_count):
            attempted_result = optimize_hdiv_mmm_magnet_from_transfer_matrices(
                dynamic_objective.orbits,
                dynamic_objective.target_matrices,
                transfer_matrix_band=dynamic_objective.transfer_matrix_band,
                bend_field_band=dynamic_objective.bend_field_band,
                charge_gram=charge_gram, fes=fes, inv_chi=inv_chi,
                rhs=rhs, field_response_matrix=response_matrix,
                incident_field_response=incident,
                field_correction=field_correction,
                active_elements=active, element_volumes=volumes,
                volume_max=volume_max,
                response_entries=dynamic_objective.response_entries,
                curvature_sign=dynamic_objective.curvature_sign,
                gradient_sign=dynamic_objective.gradient_sign,
                max_iterations=inner_count,
                initial_material_move_fraction=trial_fraction,
                maximum_material_move_fraction=trial_fraction,
                source_calibration_rows=(
                    calibration_rows if optimize_source_scale else None),
                source_calibration_target=(
                    calibration_target if optimize_source_scale else None),
                **generation_options)
            candidate_active = attempted_result.active_elements
            if np.array_equal(candidate_active, active):
                trial_reason = attempted_result.generation.stop_reason
                break
            candidate_source_scale = (
                source_scale * attempted_result.generation.source_scale)
            candidate_field = CoilHDivTotalField(
                source, charge_gram, attempted_result.generation.state,
                source_scale=candidate_source_scale, hdiv_order=hdiv_order)
            try:
                candidate_recovered = recover_ffag_closed_orbit_family(
                    candidate_field, target_family, n_segments=segment_count,
                    gradient_offset=gradient_offset,
                    initial_references=recovered, **recovery)
            except RuntimeError as exc:
                trial_reason = "orbit recovery rejected: " + str(exc)
                trial_fraction *= 0.5
                continue
            candidate_field_ratios, candidate_matrix_ratios = (
                _realized_ffag_band_ratios(
                    candidate_recovered, target_family.objective))
            candidate_ratio = float(max(
                np.max(candidate_field_ratios),
                np.max(candidate_matrix_ratios)))
            if (candidate_ratio < current_ratio - ratio_tolerance
                    or candidate_ratio <= 1.0 + ratio_tolerance):
                accepted = True
                trial_reason = "accepted after full-field orbit rebuild"
                break
            trial_reason = (
                "full-field objective did not improve "
                f"({current_ratio:.6e} -> {candidate_ratio:.6e})")
            trial_fraction *= 0.5

        if attempted_result is None:
            raise RuntimeError("FFAG outer loop did not evaluate a proposal")
        if accepted:
            before_count = int(np.count_nonzero(active))
            source_scale_before = source_scale
            active = np.asarray(candidate_active, dtype=bool).copy()
            state = attempted_result.generation.state.copy()
            source_scale = candidate_source_scale
            rhs = source_scale * source_rhs
            recovered = candidate_recovered
            previous_ratio = current_ratio
            field_ratios = candidate_field_ratios
            matrix_ratios = candidate_matrix_ratios
            current_ratio = candidate_ratio
            move_fraction = min(1.0, 1.5 * trial_fraction)
            after_count = int(np.count_nonzero(active))
        else:
            before_count = int(np.count_nonzero(active))
            after_count = before_count
            source_scale_before = source_scale
            previous_ratio = current_ratio
        history.append(FFAGHDivMMMOuterIteration(
            outer_index, trial_fraction, before_count, after_count,
            source_scale_before, source_scale,
            previous_ratio, current_ratio,
            max(value.periodic_position_residual_m for value in recovered),
            max(value.periodic_tangent_residual for value in recovered),
            accepted, trial_reason, attempted_result))
        if not accepted:
            stop_reason = trial_reason
            break
    else:
        stop_reason = "maximum outer iterations reached"

    converged = bool(current_ratio <= 1.0 + ratio_tolerance)
    if converged:
        stop_reason = "full-field orbit and transfer bands reached"
    return FFAGHDivMMMTopologyResult(
        target_family, source_scale, active, state,
        tuple(recovered), field_ratios, matrix_ratios, tuple(history),
        converged, stop_reason,
        ngsolve_growth_topology(fes.mesh, active))


__all__ = [
    "EngeFringeIntegrals",
    "FFAGCellReference",
    "FFAGCellTargetFamily",
    "FFAGHDivMMMOuterIteration",
    "FFAGHDivMMMTopologyResult",
    "FFAGFixedOrbitHDivMMMTopologyResult",
    "FFAGFixedOrbitMapTrial",
    "FFAGSoftEdgeCellSpec",
    "FullFieldClosedOrbit",
    "GEV_C_PER_TESLA_METRE",
    "PROTON_REST_ENERGY_MEV",
    "build_ffag_cell_reference",
    "build_ffag_cell_target_family",
    "enge_fringe_integrals",
    "magnetic_rigidity_from_kinetic_energy",
    "optimize_ffag_hdiv_mmm_from_fixed_design_orbits",
    "optimize_ffag_hdiv_mmm_from_transfer_matrices",
    "recover_ffag_closed_orbit",
    "recover_ffag_closed_orbit_family",
    "recover_periodic_planar_closed_orbit",
    "sample_planar_orbit_field_response",
]
