"""Design-orbit-first HDiv-MMM topology optimization for accelerator magnets.

This proof-of-concept connects a prescribed planar reference orbit and a
target first-order transfer matrix to the existing whole-element HDiv-MMM
optimizer.  The electromagnetic problem supplies row-major response rows

``[B_binormal(segment 0..n-1), dB_binormal/dnormal(segment 0..n-1)]``.

The orbit fixes the required dipole field through ``B rho * curvature``.  The
same field response is converted to a 6-by-6 combined-function transfer map,
including its analytic Frechet Jacobian.  No design finite difference, density
interpolation, or gray material is used.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .isochronous_topopt import (
    MU0,
    CombinedFunctionTransferMap,
    combined_function_transfer_map_from_field_response,
)
from .topology_optimization import (
    HDivMMMGenerationResult,
    GrowthTopologyReport,
    grow_hdiv_mmm_by_superposition,
    ngsolve_growth_topology,
)


_ALL_TRANSFER_ENTRIES = tuple(
    (row, column) for row in range(6) for column in range(6))


def _finite_array(value, *, shape=None, name):
    array = np.asarray(value, dtype=float)
    if shape is not None and array.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


@dataclass(frozen=True)
class PlanarDesignOrbit:
    """Piecewise-smooth planar design orbit used by the magnet objective.

    ``positions`` and ``tangents`` contain the same ``n+1`` orbit stations.
    ``bend_axis`` is the constant plane normal and defines the positive signed
    turning angle.  The tangent rotation divided by the chord length supplies
    one signed curvature per electromagnetic response segment.

    ``magnetic_rigidity`` is the positive reference ``B rho`` in tesla-metre.
    The charge/bend orientation is represented by the signed curvature rather
    than by a signed rigidity.
    """

    positions: np.ndarray
    tangents: np.ndarray
    magnetic_rigidity: float
    bend_axis: np.ndarray

    def __post_init__(self):
        positions = _finite_array(self.positions, name="orbit positions")
        tangents = _finite_array(self.tangents, name="orbit tangents")
        axis = _finite_array(self.bend_axis, shape=(3,), name="bend_axis")
        rigidity = float(self.magnetic_rigidity)
        if (positions.ndim != 2 or positions.shape[1] != 3
                or positions.shape[0] < 2 or tangents.shape != positions.shape):
            raise ValueError(
                "positions and tangents need matching shape (n_station,3) "
                "with at least two stations")
        tangent_norm = np.linalg.norm(tangents, axis=1)
        axis_norm = float(np.linalg.norm(axis))
        if (np.any(tangent_norm <= 0.0) or axis_norm <= 0.0
                or not np.isfinite(rigidity) or rigidity <= 0.0):
            raise ValueError(
                "orbit tangents and bend_axis must be nonzero and "
                "magnetic_rigidity must be positive")
        tangents = tangents / tangent_norm[:, None]
        axis = axis / axis_norm
        segment_lengths = np.linalg.norm(np.diff(positions, axis=0), axis=1)
        scale = max(1.0, float(np.max(np.linalg.norm(
            positions - positions[0], axis=1))))
        if np.any(segment_lengths <= 1.0e-14 * scale):
            raise ValueError("consecutive design-orbit stations must differ")
        planar_error = np.max(np.abs((positions - positions[0]) @ axis))
        tangent_error = np.max(np.abs(tangents @ axis))
        if planar_error > 1.0e-9 * scale or tangent_error > 1.0e-9:
            raise ValueError(
                "design orbit and tangents must lie in the plane normal to "
                "bend_axis")
        chord = np.diff(positions, axis=0) / segment_lengths[:, None]
        alignment = np.einsum(
            "ij,ij->i", chord, tangents[:-1] + tangents[1:])
        if np.any(alignment <= 0.0):
            raise ValueError(
                "design-orbit tangents must point from entrance to exit")
        object.__setattr__(self, "positions", positions.copy())
        object.__setattr__(self, "tangents", tangents.copy())
        object.__setattr__(self, "bend_axis", axis.copy())
        object.__setattr__(self, "magnetic_rigidity", rigidity)

    @property
    def segment_lengths(self) -> np.ndarray:
        return np.linalg.norm(np.diff(self.positions, axis=0), axis=1)

    @property
    def sample_positions(self) -> np.ndarray:
        return 0.5 * (self.positions[:-1] + self.positions[1:])

    @property
    def signed_curvature(self) -> np.ndarray:
        left = self.tangents[:-1]
        right = self.tangents[1:]
        sine = np.einsum(
            "j,ij->i", self.bend_axis, np.cross(left, right))
        cosine = np.einsum("ij,ij->i", left, right)
        turning = np.arctan2(sine, cosine)
        return turning / self.segment_lengths


def planar_orbit_field_observations(
        orbit: PlanarDesignOrbit, *, gradient_offset) -> tuple[np.ndarray,
                                                               np.ndarray]:
    """Return points and vector weights for orbit ``B``/normal-gradient rows.

    Each dipole row samples the bend-axis component at a segment midpoint.
    Each gradient row applies a centered physical-space stencil along the local
    in-plane normal.  This stencil defines the field observable only; topology
    derivatives still come from analytic HDiv-MMM Schur contractions.
    """
    if not isinstance(orbit, PlanarDesignOrbit):
        raise TypeError("orbit must be a PlanarDesignOrbit")
    offset = float(gradient_offset)
    scale = max(1.0, float(np.max(orbit.segment_lengths)))
    if not np.isfinite(offset) or offset <= 1.0e-12 * scale:
        raise ValueError("gradient_offset must be a positive physical length")
    center = orbit.sample_positions
    tangent = orbit.tangents[:-1] + orbit.tangents[1:]
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    normal = np.cross(orbit.bend_axis[None, :], tangent)
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    count = len(center)
    points = np.vstack((center, center + offset*normal,
                        center - offset*normal))
    weights = np.zeros((2*count, 3*count, 3), dtype=float)
    weights[np.arange(count), np.arange(count), :] = orbit.bend_axis
    weights[count+np.arange(count), count+np.arange(count), :] = (
        orbit.bend_axis/(2.0*offset))
    weights[count+np.arange(count), 2*count+np.arange(count), :] = (
        -orbit.bend_axis/(2.0*offset))
    return np.ascontiguousarray(points), np.ascontiguousarray(weights)


def build_planar_orbit_field_response_matrix(
        charge_gram, orbit: PlanarDesignOrbit, *, gradient_offset,
        field_scale=MU0) -> np.ndarray:
    """Build exact native HDiv observation rows at the prescribed orbit.

    The caller owns the surrounding ``ngsolve.TaskManager``.  ``field_scale``
    defaults to ``mu0`` because the configured Laplace functional returns the
    magnetic-field kernel while accelerator optics consumes tesla.
    """
    native = getattr(charge_gram, "configured_field_functional_rows", None)
    if native is None:
        raise TypeError(
            "charge_gram must expose configured_field_functional_rows")
    scale = float(field_scale)
    if not np.isfinite(scale) or scale == 0.0:
        raise ValueError("field_scale must be finite and nonzero")
    points, weights = planar_orbit_field_observations(
        orbit, gradient_offset=gradient_offset)
    rows = scale*np.asarray(native(points, weights), dtype=float)
    if (rows.ndim != 2 or rows.shape[0] != 2*len(orbit.segment_lengths)
            or not np.all(np.isfinite(rows))):
        raise RuntimeError(
            "native configured-field API returned invalid orbit rows")
    return np.ascontiguousarray(rows)


def static_magnet_symplectic_residual(matrix) -> float:
    """Infinity-norm residual for the static-magnet coordinate convention."""
    value = _finite_array(
        matrix, shape=(6, 6), name="transfer matrix")
    form = np.zeros((6, 6))
    # ell has the opposite canonical sign to x and y for the state ordering
    # (x,x',y,y',ell,delta) used by isochronous_topopt.
    for left, right, sign in ((0, 1, 1.0), (2, 3, 1.0),
                              (4, 5, -1.0)):
        form[left, right] = sign
        form[right, left] = -sign
    return float(np.linalg.norm(value.T@form@value-form, ord=np.inf))


@dataclass(frozen=True)
class PlanarTransferMatrixObjective:
    """Joint design-orbit and transfer-matrix objective.

    Matrix coordinates are ordered ``(x,x',y,y',ell,delta)``.  By default all
    36 entries are checked, so entries that are invariant in the underlying
    combined-function model still fail loudly when the requested matrix is
    incompatible.  A caller may supply a physically relevant subset through
    ``response_entries``.
    """

    orbit: PlanarDesignOrbit
    target_matrix: np.ndarray
    transfer_matrix_band: np.ndarray | float
    bend_field_band: np.ndarray | float
    response_entries: tuple[tuple[int, int], ...] = _ALL_TRANSFER_ENTRIES
    curvature_sign: float = 1.0
    gradient_sign: float = 1.0

    def __post_init__(self):
        if not isinstance(self.orbit, PlanarDesignOrbit):
            raise TypeError("orbit must be a PlanarDesignOrbit")
        matrix = _finite_array(
            self.target_matrix, shape=(6, 6), name="target transfer matrix")
        matrix_band = np.asarray(self.transfer_matrix_band, dtype=float)
        try:
            matrix_band = np.broadcast_to(matrix_band, (6, 6)).copy()
        except ValueError as exc:
            raise ValueError(
                "transfer_matrix_band must be scalar or broadcast to (6,6)"
            ) from exc
        count = len(self.orbit.segment_lengths)
        bend_band = np.asarray(self.bend_field_band, dtype=float)
        try:
            bend_band = np.broadcast_to(bend_band, (count,)).copy()
        except ValueError as exc:
            raise ValueError(
                "bend_field_band must be scalar or match the orbit segments"
            ) from exc
        entries = tuple(tuple(int(value) for value in pair)
                        for pair in self.response_entries)
        curvature_sign = float(self.curvature_sign)
        gradient_sign = float(self.gradient_sign)
        if (not entries or len(set(entries)) != len(entries)
                or any(len(pair) != 2 for pair in entries)
                or any(row < 0 or row >= 6 or column < 0 or column >= 6
                       for row, column in entries)):
            raise ValueError(
                "response_entries must contain unique zero-based 6x6 indices")
        if (not np.all(np.isfinite(matrix_band)) or np.any(matrix_band <= 0.0)
                or not np.all(np.isfinite(bend_band))
                or np.any(bend_band <= 0.0)
                or not np.isfinite(curvature_sign) or curvature_sign == 0.0
                or not np.isfinite(gradient_sign) or gradient_sign == 0.0):
            raise ValueError(
                "objective bands and field-to-optics signs must be finite; "
                "bands must be positive and signs nonzero")
        object.__setattr__(self, "target_matrix", matrix.copy())
        object.__setattr__(self, "transfer_matrix_band", matrix_band)
        object.__setattr__(self, "bend_field_band", bend_band)
        object.__setattr__(self, "response_entries", entries)
        object.__setattr__(self, "curvature_sign", curvature_sign)
        object.__setattr__(self, "gradient_sign", gradient_sign)

    @property
    def raw_field_response_size(self) -> int:
        return 2 * len(self.orbit.segment_lengths)

    @property
    def required_bend_field(self) -> np.ndarray:
        return (self.orbit.magnetic_rigidity
                * self.orbit.signed_curvature / self.curvature_sign)

    @property
    def response_target(self) -> np.ndarray:
        matrix_values = np.asarray([
            self.target_matrix[row, column]
            for row, column in self.response_entries])
        return np.r_[self.required_bend_field, matrix_values]

    @property
    def response_band(self) -> np.ndarray:
        matrix_values = np.asarray([
            self.transfer_matrix_band[row, column]
            for row, column in self.response_entries])
        return np.r_[self.bend_field_band, matrix_values]

    def evaluate_transfer_map(
            self, field_response, *, field_response_jacobian=None
            ) -> CombinedFunctionTransferMap:
        values = _finite_array(field_response, name="field response").reshape(-1)
        if values.shape != (self.raw_field_response_size,):
            raise ValueError(
                "field response must contain B and normal-gradient rows for "
                "every design-orbit segment")
        return combined_function_transfer_map_from_field_response(
            values, self.orbit.segment_lengths,
            self.orbit.magnetic_rigidity,
            field_response_jacobian=field_response_jacobian,
            curvature_sign=self.curvature_sign,
            gradient_sign=self.gradient_sign,
            response_entries=self.response_entries)

    def transform(self, field_response) -> np.ndarray:
        values = np.asarray(field_response, dtype=float).reshape(-1)
        count = len(self.orbit.segment_lengths)
        transfer = self.evaluate_transfer_map(values)
        return np.r_[values[:count], transfer.response]

    def transform_jacobian(self, field_response) -> np.ndarray:
        values = np.asarray(field_response, dtype=float).reshape(-1)
        identity = np.eye(self.raw_field_response_size)
        transfer = self.evaluate_transfer_map(
            values, field_response_jacobian=identity)
        count = len(self.orbit.segment_lengths)
        return np.vstack((identity[:count], transfer.response_jacobian))


@dataclass(frozen=True)
class TransferMatrixFieldCorrection:
    """Optics inverse result used as the target of the material inverse.

    This is the explicit boundary between accelerator optics and the
    Abe--Murata DUCAS material step.  The small dense TSVD acts only on the
    design-orbit field coordinates, and a Chebyshev LP in its retained
    subspace aligns the update with the maximum engineering-band metric.
    HDiv-MMM candidate columns do not enter this solve.
    """

    current_field_response: np.ndarray
    target_field_response: np.ndarray
    field_correction: np.ndarray
    field_response_band: np.ndarray
    current_design_response: np.ndarray
    target_design_response: np.ndarray
    linearized_design_response: np.ndarray
    nonlinear_design_response: np.ndarray
    numerical_rank: int
    singular_values: np.ndarray
    normalized_mode_target_strengths: np.ndarray
    mode_field_amplitudes: np.ndarray
    basis_coefficients: np.ndarray
    relative_linearized_residual: float
    step_scale: float
    current_max_band_ratio: float
    linearized_max_band_ratio: float
    nonlinear_max_band_ratio: float
    status: str


@dataclass(frozen=True)
class AcceleratorMagnetTopologyResult:
    """Whole-element magnet and its exact post-optimization optics report."""

    objective: PlanarTransferMatrixObjective
    generation: HDivMMMGenerationResult
    realized_field_response: np.ndarray
    realized_transfer_matrix: np.ndarray
    orbit_field_max_band_ratio: float
    transfer_matrix_max_band_ratio: float
    target_symplectic_residual: float
    realized_symplectic_residual: float
    topology: GrowthTopologyReport
    field_correction: TransferMatrixFieldCorrection | None = None

    @property
    def active_elements(self) -> np.ndarray:
        return self.generation.active_elements

    @property
    def converged(self) -> bool:
        return bool(
            self.orbit_field_max_band_ratio <= 1.0
            and self.transfer_matrix_max_band_ratio <= 1.0)

    @property
    def field_target_converged(self) -> bool:
        return self.generation.converged


@dataclass(frozen=True)
class MultiMomentumTransferMatrixObjective:
    """Joint orbit/map contract at several magnetic rigidities.

    Every orbit owns its physical observation points and may contain a
    different number of longitudinal segments.  Raw rows are concatenated in
    momentum order, with each orbit retaining the row order
    ``[B_binormal..., dB_binormal/dnormal...]``.  The transformed response is
    the corresponding concatenation of
    :class:`PlanarTransferMatrixObjective` responses.  Its Jacobian is block
    diagonal and uses the same analytic matrix-exponential Frechet chain; no
    momentum or design finite difference is introduced.
    """

    orbits: tuple[PlanarDesignOrbit, ...]
    target_matrices: np.ndarray
    transfer_matrix_band: np.ndarray | float
    bend_field_band: object
    response_entries: tuple[tuple[int, int], ...] = _ALL_TRANSFER_ENTRIES
    curvature_sign: float = 1.0
    gradient_sign: float = 1.0

    def __post_init__(self):
        orbits = tuple(self.orbits)
        if not orbits or any(
                not isinstance(orbit, PlanarDesignOrbit) for orbit in orbits):
            raise TypeError(
                "orbits must be a non-empty sequence of PlanarDesignOrbit")
        matrices = _finite_array(
            self.target_matrices, name="target_matrices")
        if matrices.shape != (len(orbits), 6, 6):
            raise ValueError(
                "target_matrices must have shape (n_momentum,6,6)")
        matrix_band = np.asarray(self.transfer_matrix_band, dtype=float)
        try:
            matrix_band = np.broadcast_to(
                matrix_band, matrices.shape).copy()
        except ValueError as exc:
            raise ValueError(
                "transfer_matrix_band must broadcast to "
                "(n_momentum,6,6)") from exc

        raw_bend_band = self.bend_field_band
        if np.isscalar(raw_bend_band):
            bend_bands = tuple(float(raw_bend_band) for _ in orbits)
        else:
            try:
                bend_bands = tuple(raw_bend_band)
            except TypeError as exc:
                raise ValueError(
                    "bend_field_band must be scalar or have one entry per "
                    "momentum") from exc
            if len(bend_bands) != len(orbits):
                raise ValueError(
                    "bend_field_band must be scalar or have one entry per "
                    "momentum")
        objectives = tuple(
            PlanarTransferMatrixObjective(
                orbit, matrix, band, bend_band, self.response_entries,
                self.curvature_sign, self.gradient_sign)
            for orbit, matrix, band, bend_band in zip(
                orbits, matrices, matrix_band, bend_bands))
        object.__setattr__(self, "orbits", orbits)
        object.__setattr__(self, "target_matrices", matrices.copy())
        object.__setattr__(self, "transfer_matrix_band", matrix_band)
        object.__setattr__(self, "bend_field_band", tuple(
            objective.bend_field_band.copy() for objective in objectives))
        object.__setattr__(self, "response_entries",
                           objectives[0].response_entries)
        object.__setattr__(self, "_objectives", objectives)

    @property
    def objectives(self) -> tuple[PlanarTransferMatrixObjective, ...]:
        return self._objectives

    @property
    def raw_field_response_size(self) -> int:
        return sum(
            objective.raw_field_response_size
            for objective in self.objectives)

    @property
    def raw_offsets(self) -> np.ndarray:
        sizes = [objective.raw_field_response_size
                 for objective in self.objectives]
        return np.r_[0, np.cumsum(sizes)].astype(np.int64)

    @property
    def response_target(self) -> np.ndarray:
        return np.concatenate([
            objective.response_target for objective in self.objectives])

    @property
    def response_band(self) -> np.ndarray:
        return np.concatenate([
            objective.response_band for objective in self.objectives])

    def split_raw_response(self, field_response) -> tuple[np.ndarray, ...]:
        values = _finite_array(
            field_response, name="multi-momentum field response").reshape(-1)
        if values.shape != (self.raw_field_response_size,):
            raise ValueError(
                "field response does not match the multi-momentum raw-row "
                "contract")
        offsets = self.raw_offsets
        return tuple(values[left:right]
                     for left, right in zip(offsets[:-1], offsets[1:]))

    def transform(self, field_response) -> np.ndarray:
        return np.concatenate([
            objective.transform(values)
            for objective, values in zip(
                self.objectives, self.split_raw_response(field_response))])

    def transform_jacobian(self, field_response) -> np.ndarray:
        values = self.split_raw_response(field_response)
        blocks = [objective.transform_jacobian(value)
                  for objective, value in zip(self.objectives, values)]
        rows = sum(block.shape[0] for block in blocks)
        jacobian = np.zeros((rows, self.raw_field_response_size))
        raw_offsets = self.raw_offsets
        response_offset = 0
        for block, left, right in zip(
                blocks, raw_offsets[:-1], raw_offsets[1:]):
            next_response = response_offset + block.shape[0]
            jacobian[response_offset:next_response, left:right] = block
            response_offset = next_response
        return jacobian


def solve_transfer_matrix_field_correction(
        objective, current_field_response, *, field_basis=None,
        relative_tolerance=1.0e-3, maximum_step_scale=1.0,
        line_search_steps=8) -> TransferMatrixFieldCorrection:
    """Invert transfer-matrix error to a target orbit-field correction.

    The current recovered design orbit is fixed during this small optics
    inverse.  If ``F(B)`` is the bend-field/transfer-map response, the method
    solves

    ``diag(1/band) dF/dB field_basis dq = diag(1/band) (target-F(B))``

    by dense TSVD followed by a Chebyshev minimax solve in the retained modal
    subspace.  A cheap nonlinear line search evaluates ``F`` itself; it
    performs no particle tracking and no HDiv-MMM candidate solve.  The
    returned field target is subsequently passed to the separate
    ACA--QR--TSVD DUCAS material inverse.

    ``field_basis`` may restrict the correction to smooth body, gradient, or
    entrance/exit fringe modes.  Identity is the explicit default and keeps
    one coordinate for every sampled binormal-field/normal-gradient value.
    """
    required = (
        "raw_field_response_size", "response_target", "response_band",
        "transform", "transform_jacobian")
    if any(not hasattr(objective, name) for name in required):
        raise TypeError(
            "objective must expose the planar transfer-matrix field contract")
    raw_size = int(objective.raw_field_response_size)
    current = _finite_array(
        current_field_response, name="current_field_response").reshape(-1)
    target = _finite_array(
        objective.response_target, name="target_design_response").reshape(-1)
    band = _finite_array(
        objective.response_band, name="design_response_band").reshape(-1)
    if (current.shape != (raw_size,) or target.size == 0
            or band.shape != target.shape or np.any(band <= 0.0)):
        raise ValueError(
            "field response and positive design-response bands must match "
            "the transfer-matrix objective")
    if field_basis is None:
        basis = np.eye(raw_size)
    else:
        basis = _finite_array(field_basis, name="field_basis")
        if basis.ndim != 2 or basis.shape[0] != raw_size or basis.shape[1] == 0:
            raise ValueError(
                "field_basis must have shape (raw_field_response_size,n_mode)")
    tolerance = float(relative_tolerance)
    maximum_scale = float(maximum_step_scale)
    line_search_steps = int(line_search_steps)
    if (not np.isfinite(tolerance) or not 0.0 < tolerance < 1.0
            or not np.isfinite(maximum_scale)
            or not 0.0 < maximum_scale <= 1.0
            or line_search_steps < 1):
        raise ValueError(
            "field inverse requires 0<tolerance<1, 0<step<=1, and a "
            "positive line-search count")

    current_design = _finite_array(
        objective.transform(current), name="current_design_response"
    ).reshape(-1)
    jacobian = _finite_array(
        objective.transform_jacobian(current),
        name="field_to_transfer_jacobian")
    if (current_design.shape != target.shape
            or jacobian.shape != (target.size, raw_size)):
        raise RuntimeError(
            "transfer objective returned an incompatible response or Jacobian")
    normalized_residual = (target - current_design) / band
    normalized_operator = (jacobian @ basis) / band[:, None]
    U, singular, Vh = np.linalg.svd(
        normalized_operator, full_matrices=False)
    if singular.size and singular[0] > 0.0:
        rank = int(np.count_nonzero(
            singular >= tolerance * singular[0]))
    else:
        rank = 0
    mode_projection = U.T @ normalized_residual
    normalized_mode_target_strengths = (
        mode_projection / np.sqrt(float(target.size)))
    mode_field_amplitudes = np.zeros_like(mode_projection)
    coefficients = np.zeros(basis.shape[1], dtype=float)
    if rank:
        # Least squares can increase the largest engineering-band error even
        # while reducing its aggregate two-norm.  That failure is common when
        # several transfer entries compete for one field mode: the subsequent
        # max-band line search then rejects every step.  Keep TSVD as the
        # reachable optics subspace, but choose its modal response by a small
        # Chebyshev LP.  The zero correction is feasible, so this stage never
        # predicts a worse linearized infinity norm.
        from scipy.optimize import linprog

        retained_left = U[:, :rank]
        lp_objective = np.r_[np.zeros(rank), 1.0]
        lp_matrix = np.vstack((
            np.c_[-retained_left, -np.ones(target.size)],
            np.c_[retained_left, -np.ones(target.size)],
        ))
        lp_rhs = np.r_[-normalized_residual, normalized_residual]
        lp = linprog(
            lp_objective, A_ub=lp_matrix, b_ub=lp_rhs,
            bounds=[(None, None)] * rank + [(0.0, None)],
            method="highs")
        if not lp.success:
            raise RuntimeError(
                "transfer-matrix TSVD Chebyshev solve failed: "
                + str(lp.message))
        retained_response_amplitudes = np.asarray(
            lp.x[:rank], dtype=float)
        mode_field_amplitudes[:rank] = (
            retained_response_amplitudes / singular[:rank])
        coefficients = (
            Vh[:rank].T @ mode_field_amplitudes[:rank])
    unconstrained_correction = basis @ coefficients
    if not np.all(np.isfinite(unconstrained_correction)):
        raise RuntimeError(
            "transfer-matrix field inverse produced a non-finite correction")

    def max_ratio(design_response):
        return float(np.max(np.abs(
            (np.asarray(design_response, dtype=float) - target) / band)))

    current_ratio = max_ratio(current_design)
    selected_scale = 0.0
    nonlinear_design = current_design.copy()
    nonlinear_ratio = current_ratio
    if current_ratio > 1.0:
        for index in range(line_search_steps):
            scale = maximum_scale * (0.5 ** index)
            candidate = _finite_array(
                objective.transform(current + scale * unconstrained_correction),
                name="line_search_design_response").reshape(-1)
            candidate_ratio = max_ratio(candidate)
            if candidate_ratio < nonlinear_ratio:
                selected_scale = scale
                nonlinear_design = candidate
                nonlinear_ratio = candidate_ratio
    correction = selected_scale * unconstrained_correction
    coefficients = selected_scale * coefficients
    mode_field_amplitudes = selected_scale * mode_field_amplitudes
    linearized_design = current_design + jacobian @ correction
    linearized_ratio = max_ratio(linearized_design)
    normalized_linearized_residual = (
        (target - linearized_design) / band)
    residual_norm = float(np.linalg.norm(normalized_residual))
    relative_residual = float(
        np.linalg.norm(normalized_linearized_residual)
        / max(np.finfo(float).tiny, residual_norm))

    # One raw-field change equal to field_response_band produces at most one
    # normalized design-band change in the local Jacobian.  This supplies the
    # separate DUCAS material inverse with physically interpretable B/G bands
    # without mixing tesla and tesla/metre by their numerical magnitudes.
    sensitivity = np.max(np.abs(jacobian / band[:, None]), axis=0)
    fallback_band = np.maximum(1.0, np.abs(current))
    field_response_band = np.divide(
        1.0, sensitivity, out=fallback_band.copy(), where=sensitivity > 0.0)
    status = (
        "target already within transfer bands" if current_ratio <= 1.0 else
        "dense optics TSVD-Chebyshev solve produced no improving field "
        "correction"
        if selected_scale == 0.0 else
        "dense optics TSVD-Chebyshev solve to orbit-field target")
    return TransferMatrixFieldCorrection(
        current_field_response=current.copy(),
        target_field_response=np.asarray(current + correction, dtype=float),
        field_correction=np.asarray(correction, dtype=float),
        field_response_band=np.asarray(field_response_band, dtype=float),
        current_design_response=current_design.copy(),
        target_design_response=target.copy(),
        linearized_design_response=np.asarray(
            linearized_design, dtype=float),
        nonlinear_design_response=np.asarray(nonlinear_design, dtype=float),
        numerical_rank=rank,
        singular_values=np.asarray(singular, dtype=float),
        normalized_mode_target_strengths=np.asarray(
            normalized_mode_target_strengths, dtype=float),
        mode_field_amplitudes=np.asarray(
            mode_field_amplitudes, dtype=float),
        basis_coefficients=np.asarray(coefficients, dtype=float),
        relative_linearized_residual=relative_residual,
        step_scale=float(selected_scale),
        current_max_band_ratio=current_ratio,
        linearized_max_band_ratio=linearized_ratio,
        nonlinear_max_band_ratio=nonlinear_ratio,
        status=status)


def build_multi_orbit_field_response_matrix(
        charge_gram, objective: MultiMomentumTransferMatrixObjective, *,
        gradient_offset, field_scale=MU0) -> np.ndarray:
    """Build native HDiv rows for every momentum-indexed orbit.

    The caller owns ``ngsolve.TaskManager``.  Each orbit is submitted as one
    native row batch, avoiding a dense block-diagonal point-weight tensor
    whose zero blocks would grow quadratically with momentum count.
    """
    if not isinstance(objective, MultiMomentumTransferMatrixObjective):
        raise TypeError(
            "objective must be a MultiMomentumTransferMatrixObjective")
    offsets = gradient_offset
    if np.isscalar(offsets):
        offsets = (float(offsets),) * len(objective.orbits)
    else:
        offsets = tuple(offsets)
        if len(offsets) != len(objective.orbits):
            raise ValueError(
                "gradient_offset must be scalar or have one value per orbit")
    rows = [build_planar_orbit_field_response_matrix(
                charge_gram, orbit, gradient_offset=offset,
                field_scale=field_scale)
            for orbit, offset in zip(objective.orbits, offsets)]
    return np.ascontiguousarray(np.vstack(rows))


@dataclass(frozen=True)
class CoilBuilderHDivSource:
    """One or more closed CoilBuilder paths as an HDiv incident source.

    The same finite-filament representation owns both sides of the coupling:
    its differentiable NGSolve ``CoefficientFunction`` assembles the iron RHS,
    while the vectorized analytic Biot--Savart evaluator supplies the incident
    vacuum field rows and full-field tracker.  This prevents a solid-current
    model and a different filament model from silently driving the solve and
    acceptance calculation.
    """

    segment_groups: tuple[tuple[np.ndarray, float], ...]

    def __post_init__(self):
        groups = []
        for index, group in enumerate(tuple(self.segment_groups)):
            if len(group) != 2:
                raise ValueError(
                    "each coil source group must contain segments and current")
            segments = np.asarray(group[0], dtype=float)
            current = float(group[1])
            if (segments.ndim != 3 or segments.shape[1:] != (2, 3)
                    or segments.shape[0] == 0
                    or not np.all(np.isfinite(segments))
                    or not np.isfinite(current) or current == 0.0):
                raise ValueError(
                    f"invalid finite-filament coil group {index}")
            if np.any(np.linalg.norm(
                    segments[:, 1] - segments[:, 0], axis=1) <= 0.0):
                raise ValueError("coil filaments must have positive length")
            groups.append((np.ascontiguousarray(segments), current))
        if not groups:
            raise ValueError("at least one coil source group is required")
        object.__setattr__(self, "segment_groups", tuple(groups))

    @classmethod
    def from_coilbuilders(cls, coils, *, n_arc=80,
                          closure_tolerance=1.0e-9):
        """Create the shared source from closed CoilBuilder objects."""
        from .coil_builder import CoilBuilder

        builders = (coils,) if isinstance(coils, CoilBuilder) else tuple(coils)
        n_arc = int(n_arc)
        tolerance = float(closure_tolerance)
        if (not builders or n_arc < 1 or not np.isfinite(tolerance)
                or tolerance < 0.0
                or any(not isinstance(value, CoilBuilder)
                       for value in builders)):
            raise ValueError(
                "coils must contain CoilBuilder objects, n_arc must be "
                "positive, and closure_tolerance nonnegative")
        groups = []
        for index, builder in enumerate(builders):
            if float(builder.gap) > tolerance:
                raise ValueError(
                    f"CoilBuilder path {index} is open by {builder.gap:.6e} m")
            segments, current = builder.to_wire_segments(n_arc=n_arc)
            groups.append((np.asarray(segments, dtype=float), float(current)))
        return cls(tuple(groups))

    @property
    def segment_count(self) -> int:
        return sum(len(segments) for segments, _ in self.segment_groups)

    def h_field(self, points) -> np.ndarray:
        """Evaluate the source ``H`` field in A/m at one or many points."""
        from .biot_savart import h_segments_batch

        raw = np.asarray(points, dtype=float)
        single = raw.shape == (3,)
        values = np.ascontiguousarray(raw.reshape(-1, 3))
        if not np.all(np.isfinite(values)):
            raise ValueError("coil field points must be finite")
        field = np.zeros((len(values), 3), dtype=float)
        for segments, current in self.segment_groups:
            field += h_segments_batch(segments, values, current=current)
        return field[0] if single else field

    def b_field(self, points) -> np.ndarray:
        """Evaluate the source magnetic flux density in tesla."""
        return MU0 * self.h_field(points)

    def coefficient_function(self):
        """Return the differentiable NGSolve source ``H`` field."""
        from .biot_savart import h_segments_cf

        fields = [h_segments_cf(segments, current=current)
                  for segments, current in self.segment_groups]
        while len(fields) > 1:
            fields = [
                fields[index] + fields[index + 1]
                if index + 1 < len(fields) else fields[index]
                for index in range(0, len(fields), 2)]
        return fields[0]

    def assemble_hdiv_rhs(self, fes, *, bonus_intorder=10) -> np.ndarray:
        """Assemble ``integral H_coil . v``; caller owns TaskManager."""
        import ngsolve as ng

        form = ng.LinearForm(fes)
        form += ng.InnerProduct(
            self.coefficient_function(), fes.TestFunction()) * ng.dx(
                bonus_intorder=int(bonus_intorder))
        form.Assemble()
        values = np.asarray(form.vec.FV().NumPy(), dtype=float).copy()
        if values.shape != (int(fes.ndof),) or not np.all(np.isfinite(values)):
            raise RuntimeError("CoilBuilder HDiv source assembly failed")
        return values

    def incident_orbit_field_response(
            self, objective: MultiMomentumTransferMatrixObjective, *,
            gradient_offset) -> np.ndarray:
        """Sample the same CoilBuilder source on all objective orbits."""
        if not isinstance(objective, MultiMomentumTransferMatrixObjective):
            raise TypeError(
                "objective must be a MultiMomentumTransferMatrixObjective")
        offsets = gradient_offset
        if np.isscalar(offsets):
            offsets = (float(offsets),) * len(objective.orbits)
        else:
            offsets = tuple(offsets)
            if len(offsets) != len(objective.orbits):
                raise ValueError(
                    "gradient_offset must be scalar or have one value per "
                    "orbit")
        responses = []
        for orbit, offset in zip(objective.orbits, offsets):
            points, weights = planar_orbit_field_observations(
                orbit, gradient_offset=offset)
            responses.append(np.einsum(
                "rpc,pc->r", weights, self.b_field(points)))
        result = np.ascontiguousarray(np.concatenate(responses))
        if result.shape != (objective.raw_field_response_size,):
            raise RuntimeError("CoilBuilder orbit response size mismatch")
        return result


class CoilHDivTotalField:
    """Persistent total ``B`` evaluator for one solved active iron state."""

    def __init__(self, source: CoilBuilderHDivSource, charge_gram, state, *,
                 source_scale=1.0, hdiv_order=1, algorithm="auto"):
        from .vim._field_batch import _create_field_evaluator

        if not isinstance(source, CoilBuilderHDivSource):
            raise TypeError("source must be a CoilBuilderHDivSource")
        coefficients = np.ascontiguousarray(state, dtype=float).reshape(-1)
        scale = float(source_scale)
        order = int(hdiv_order)
        if (coefficients.size == 0 or not np.all(np.isfinite(coefficients))
                or not np.isfinite(scale) or scale <= 0.0
                or order not in (1, 2)):
            raise ValueError("invalid solved HDiv total-field state")
        evaluator, stats = _create_field_evaluator(
            charge_gram, coefficients, order)
        self.source = source
        self.charge_gram = charge_gram
        self.state = coefficients.copy()
        self.source_scale = scale
        self.hdiv_order = order
        self.algorithm = str(algorithm)
        self._evaluator = evaluator
        self.stats = stats

    def b_field(self, points) -> np.ndarray:
        raw = np.asarray(points, dtype=float)
        single = raw.shape == (3,)
        values = np.ascontiguousarray(raw.reshape(-1, 3))
        demag_h = np.asarray(
            self._evaluator.field(values, self.algorithm), dtype=float
        ) / (4.0 * np.pi)
        total = MU0 * (
            self.source_scale * self.source.h_field(values) + demag_h)
        return total[0] if single else total

    def __call__(self, x, y=None, z=None):
        point = (np.asarray(x, dtype=float) if y is None and z is None else
                 np.asarray([x, y, z], dtype=float))
        return self.b_field(point)


@dataclass(frozen=True)
class MultiMomentumAcceleratorMagnetTopologyResult:
    """Exact whole-element result checked at every requested rigidity."""

    objective: MultiMomentumTransferMatrixObjective
    generation: HDivMMMGenerationResult
    realized_field_responses: tuple[np.ndarray, ...]
    realized_transfer_matrices: np.ndarray
    orbit_field_max_band_ratios: np.ndarray
    transfer_matrix_max_band_ratios: np.ndarray
    target_symplectic_residuals: np.ndarray
    realized_symplectic_residuals: np.ndarray
    topology: GrowthTopologyReport
    field_correction: TransferMatrixFieldCorrection | None = None

    @property
    def active_elements(self) -> np.ndarray:
        return self.generation.active_elements

    @property
    def converged(self) -> bool:
        return bool(
            np.max(self.orbit_field_max_band_ratios) <= 1.0
            and np.max(self.transfer_matrix_max_band_ratios) <= 1.0)

    @property
    def field_target_converged(self) -> bool:
        return self.generation.converged


def optimize_hdiv_mmm_magnet_from_transfer_matrix(
        design_orbit: PlanarDesignOrbit, target_transfer_matrix, *,
        transfer_matrix_band, bend_field_band,
        charge_gram, fes, inv_chi, rhs, field_response_matrix,
        active_elements, element_volumes, volume_max,
        incident_field_response=None, field_correction=None,
        response_entries=None,
        curvature_sign=1.0, gradient_sign=1.0,
        **generation_options) -> AcceleratorMagnetTopologyResult:
    """Create a whole-element HDiv-MMM magnet from orbit and map inputs.

    CAD/mesh construction and coil excitation remain caller-owned.  The field
    response matrix must already sample the total binormal field and its normal
    gradient at ``design_orbit.sample_positions`` using NGSolve-owned FE
    evaluation.  This function owns the physics-to-optics chain and the binary
    material optimization, returning the active HEX/TET/WEDGE elements that
    constitute the PoC magnet.
    """
    entries = (_ALL_TRANSFER_ENTRIES if response_entries is None else
               tuple(response_entries))
    objective = PlanarTransferMatrixObjective(
        design_orbit, target_transfer_matrix, transfer_matrix_band,
        bend_field_band, entries, curvature_sign, gradient_sign)
    response_matrix = _finite_array(
        field_response_matrix, name="field_response_matrix")
    if (response_matrix.ndim != 2
            or response_matrix.shape != (
                objective.raw_field_response_size, int(fes.ndof))):
        raise ValueError(
            "field_response_matrix must have shape "
            "(2*n_orbit_segment,fes.ndof)")
    incident = (np.zeros(objective.raw_field_response_size)
                if incident_field_response is None else
                _finite_array(
                    incident_field_response,
                    name="incident_field_response").reshape(-1))
    if incident.shape != (objective.raw_field_response_size,):
        raise ValueError(
            "incident_field_response must match the raw field response size")
    reserved = {
        "response_matrix", "response_target", "response_band",
        "response_transform", "response_transform_jacobian",
        "incident_response",
    }
    overlap = reserved.intersection(generation_options)
    if overlap:
        raise TypeError(
            "generation_options cannot override the orbit/map contract: "
            + ", ".join(sorted(overlap)))
    if field_correction is None:
        generation_target = objective.response_target
        generation_band = objective.response_band
        response_transform = objective.transform
        response_transform_jacobian = objective.transform_jacobian
    else:
        if not isinstance(field_correction, TransferMatrixFieldCorrection):
            raise TypeError(
                "field_correction must be a TransferMatrixFieldCorrection")
        generation_target = _finite_array(
            field_correction.target_field_response,
            name="field_correction.target_field_response").reshape(-1)
        generation_band = _finite_array(
            field_correction.field_response_band,
            name="field_correction.field_response_band").reshape(-1)
        if (generation_target.shape != (objective.raw_field_response_size,)
                or generation_band.shape != generation_target.shape
                or np.any(generation_band <= 0.0)):
            raise ValueError(
                "field correction must match the raw orbit field rows")
        response_transform = None
        response_transform_jacobian = None
    generation = grow_hdiv_mmm_by_superposition(
        charge_gram=charge_gram, fes=fes, inv_chi=inv_chi, rhs=rhs,
        response_matrix=response_matrix, active_elements=active_elements,
        element_volumes=element_volumes,
        response_target=generation_target,
        response_band=generation_band, volume_max=volume_max,
        incident_response=incident,
        response_transform=response_transform,
        response_transform_jacobian=response_transform_jacobian,
        **generation_options)
    raw = response_matrix @ generation.state + incident
    transfer = objective.evaluate_transfer_map(raw)
    count = len(design_orbit.segment_lengths)
    orbit_ratio = float(np.max(np.abs(
        (raw[:count] - objective.required_bend_field)
        / objective.bend_field_band)))
    target_values = np.asarray([
        objective.target_matrix[row, column]
        for row, column in objective.response_entries])
    target_band = np.asarray([
        objective.transfer_matrix_band[row, column]
        for row, column in objective.response_entries])
    transfer_ratio = float(np.max(np.abs(
        (transfer.response - target_values) / target_band)))
    topology = ngsolve_growth_topology(fes.mesh, generation.active_elements)
    return AcceleratorMagnetTopologyResult(
        objective, generation, np.asarray(raw, dtype=float),
        np.asarray(transfer.matrix, dtype=float), orbit_ratio,
        transfer_ratio,
        static_magnet_symplectic_residual(objective.target_matrix),
        static_magnet_symplectic_residual(transfer.matrix), topology,
        field_correction)


def optimize_hdiv_mmm_magnet_from_transfer_matrices(
        design_orbits, target_transfer_matrices, *,
        transfer_matrix_band, bend_field_band,
        charge_gram, fes, inv_chi, rhs, field_response_matrix,
        active_elements, element_volumes, volume_max,
        incident_field_response=None, field_correction=None,
        response_entries=None,
        curvature_sign=1.0, gradient_sign=1.0,
        **generation_options) -> MultiMomentumAcceleratorMagnetTopologyResult:
    """Create one binary HDiv-MMM magnet for several orbit/map targets.

    This is the FFAG operating-point fusion API.  The field-response matrix
    contains all momentum-indexed vacuum observation rows, while the physical
    solve and whole-element add/remove proposal remain shared.  The nonlinear
    optics transform and its exact Frechet Jacobian are block-assembled before
    the existing batched adjoint contraction and ACA--QR--TSVD master step.
    """
    entries = (_ALL_TRANSFER_ENTRIES if response_entries is None else
               tuple(response_entries))
    objective = MultiMomentumTransferMatrixObjective(
        tuple(design_orbits), target_transfer_matrices,
        transfer_matrix_band, bend_field_band, entries,
        curvature_sign, gradient_sign)
    response_matrix = _finite_array(
        field_response_matrix, name="field_response_matrix")
    if (response_matrix.ndim != 2
            or response_matrix.shape != (
                objective.raw_field_response_size, int(fes.ndof))):
        raise ValueError(
            "field_response_matrix must have shape "
            "(sum(2*n_orbit_segment),fes.ndof)")
    incident = (np.zeros(objective.raw_field_response_size)
                if incident_field_response is None else
                _finite_array(
                    incident_field_response,
                    name="incident_field_response").reshape(-1))
    if incident.shape != (objective.raw_field_response_size,):
        raise ValueError(
            "incident_field_response must match all raw orbit rows")
    reserved = {
        "response_matrix", "response_target", "response_band",
        "response_transform", "response_transform_jacobian",
        "incident_response",
    }
    overlap = reserved.intersection(generation_options)
    if overlap:
        raise TypeError(
            "generation_options cannot override the multi-momentum "
            "orbit/map contract: " + ", ".join(sorted(overlap)))
    if field_correction is None:
        generation_target = objective.response_target
        generation_band = objective.response_band
        response_transform = objective.transform
        response_transform_jacobian = objective.transform_jacobian
    else:
        if not isinstance(field_correction, TransferMatrixFieldCorrection):
            raise TypeError(
                "field_correction must be a TransferMatrixFieldCorrection")
        generation_target = _finite_array(
            field_correction.target_field_response,
            name="field_correction.target_field_response").reshape(-1)
        generation_band = _finite_array(
            field_correction.field_response_band,
            name="field_correction.field_response_band").reshape(-1)
        if (generation_target.shape != (objective.raw_field_response_size,)
                or generation_band.shape != generation_target.shape
                or np.any(generation_band <= 0.0)):
            raise ValueError(
                "field correction must match all raw orbit field rows")
        response_transform = None
        response_transform_jacobian = None
    generation = grow_hdiv_mmm_by_superposition(
        charge_gram=charge_gram, fes=fes, inv_chi=inv_chi, rhs=rhs,
        response_matrix=response_matrix, active_elements=active_elements,
        element_volumes=element_volumes,
        response_target=generation_target,
        response_band=generation_band, volume_max=volume_max,
        incident_response=incident,
        response_transform=response_transform,
        response_transform_jacobian=response_transform_jacobian,
        **generation_options)
    raw = response_matrix @ generation.state + incident
    split = objective.split_raw_response(raw)
    transfers = tuple(
        item.evaluate_transfer_map(values)
        for item, values in zip(objective.objectives, split))
    orbit_ratios = []
    transfer_ratios = []
    for item, values, transfer in zip(
            objective.objectives, split, transfers):
        count = len(item.orbit.segment_lengths)
        orbit_ratios.append(float(np.max(np.abs(
            (values[:count] - item.required_bend_field)
            / item.bend_field_band))))
        target_values = np.asarray([
            item.target_matrix[row, column]
            for row, column in item.response_entries])
        target_band = np.asarray([
            item.transfer_matrix_band[row, column]
            for row, column in item.response_entries])
        transfer_ratios.append(float(np.max(np.abs(
            (transfer.response - target_values) / target_band))))
    target_symplectic = np.asarray([
        static_magnet_symplectic_residual(matrix)
        for matrix in objective.target_matrices])
    realized_matrices = np.asarray([
        transfer.matrix for transfer in transfers])
    realized_symplectic = np.asarray([
        static_magnet_symplectic_residual(matrix)
        for matrix in realized_matrices])
    topology = ngsolve_growth_topology(fes.mesh, generation.active_elements)
    return MultiMomentumAcceleratorMagnetTopologyResult(
        objective, generation, tuple(np.asarray(value, dtype=float)
                                     for value in split),
        realized_matrices, np.asarray(orbit_ratios),
        np.asarray(transfer_ratios), target_symplectic,
        realized_symplectic, topology, field_correction)


__all__ = [
    "AcceleratorMagnetTopologyResult",
    "CoilBuilderHDivSource",
    "CoilHDivTotalField",
    "MultiMomentumAcceleratorMagnetTopologyResult",
    "MultiMomentumTransferMatrixObjective",
    "PlanarDesignOrbit",
    "PlanarTransferMatrixObjective",
    "TransferMatrixFieldCorrection",
    "build_multi_orbit_field_response_matrix",
    "build_planar_orbit_field_response_matrix",
    "optimize_hdiv_mmm_magnet_from_transfer_matrix",
    "optimize_hdiv_mmm_magnet_from_transfer_matrices",
    "planar_orbit_field_observations",
    "solve_transfer_matrix_field_correction",
    "static_magnet_symplectic_residual",
]
