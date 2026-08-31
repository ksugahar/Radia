"""Design-orbit-first HDiv-MMM topology optimization for accelerator magnets.

This proof-of-concept connects a prescribed planar reference orbit and a
target first-order transfer matrix to the existing whole-element HDiv-MMM
optimizer.  The electromagnetic problem supplies row-major response rows

``[B_binormal(segment 0..n-1), dB_binormal/dnormal(segment 0..n-1)]``.

The orbit fixes the required dipole field through ``B rho * curvature``.  The
same field response is converted to a 6-by-6 combined-function transfer map,
including its forward-mode AD Jacobian.  The matrix exponential is an exact
Frechet-differentiated primitive.  No design finite difference, density
interpolation, or gray material is used.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np

from .isochronous_topopt import (
    MU0,
    CombinedFunctionTransferMap,
    combined_function_transfer_map_from_field_response,
)
from .topology_optimization import (
    GrowthTopologyReport,
    HDivMMMGenerationResult,
    TSVDElementCandidateSelection,
    grow_hdiv_mmm_by_superposition,
    ngsolve_growth_topology,
    select_tsvd_element_candidates,
)

_ALL_TRANSFER_ENTRIES = tuple(
    (row, column) for row in range(6) for column in range(6))

_STATIC_MAGNET_TRANSFER_COMPONENT_GROUPS = {
    # State ordering: (x,x',y,y',ell,delta).  The two focusing-strength
    # entries are the most direct first PoC for pole-face control.  Complete
    # transverse blocks remain available when phase advance and imaging must
    # be constrained together.
    "horizontal_focusing": ((1, 0),),
    "vertical_focusing": ((3, 2),),
    "horizontal_block": ((0, 0), (0, 1), (1, 0), (1, 1)),
    "vertical_block": ((2, 2), (2, 3), (3, 2), (3, 3)),
    "horizontal_dispersion": ((0, 5), (1, 5)),
    "path_length": ((4, 0), (4, 1), (4, 5)),
}


def static_magnet_transfer_component_entries(component_groups):
    """Resolve named, physically interpretable transfer-map objectives.

    The returned zero-based entries use the static-magnet state ordering
    ``(x,x',y,y',ell,delta)``.  This selects which entries are compared; it
    does not make the remaining entries independent.  Every realized map is
    still produced by the Hamiltonian field model and remains symplectic.
    """
    if isinstance(component_groups, str):
        component_groups = (component_groups,)
    try:
        names = tuple(
            str(value).strip().lower().replace("-", "_")
            for value in component_groups)
    except TypeError as exc:
        raise ValueError(
            "component_groups must be a non-empty sequence of names") from exc
    if not names or any(not name for name in names):
        raise ValueError(
            "component_groups must be a non-empty sequence of names")
    unknown = tuple(
        name for name in names
        if name not in _STATIC_MAGNET_TRANSFER_COMPONENT_GROUPS)
    if unknown:
        available = ", ".join(_STATIC_MAGNET_TRANSFER_COMPONENT_GROUPS)
        raise ValueError(
            f"unknown transfer component group {unknown[0]!r}; "
            f"available groups are {available}")
    entries = []
    for name in names:
        for entry in _STATIC_MAGNET_TRANSFER_COMPONENT_GROUPS[name]:
            if entry not in entries:
                entries.append(entry)
    return tuple(entries)


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
    ``path_length_stations`` optionally preserves the RK independent variable;
    when present it is the source of truth for each ``Delta s``.  Otherwise
    chord lengths are used for backward-compatible manually supplied orbits.
    ``bend_axis`` is the constant plane normal and defines the positive signed
    turning angle.  By default, tangent rotation divided by segment length
    supplies the signed curvature.  A tracked orbit may instead provide
    ``signed_curvature_per_m`` from the ODE/B-field collocation points; the
    tangent-turning value remains available as ``geometric_signed_curvature``
    for a finite-station consistency diagnostic.

    ``magnetic_rigidity`` is the positive reference ``B rho`` in tesla-metre.
    The charge/bend orientation is represented by the signed curvature rather
    than by a signed rigidity.
    """

    positions: np.ndarray
    tangents: np.ndarray
    magnetic_rigidity: float
    bend_axis: np.ndarray
    path_length_stations: np.ndarray | None = None
    signed_curvature_per_m: np.ndarray | None = None

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
        chord_lengths = np.linalg.norm(np.diff(positions, axis=0), axis=1)
        scale = max(1.0, float(np.max(np.linalg.norm(
            positions - positions[0], axis=1))))
        if np.any(chord_lengths <= 1.0e-14 * scale):
            raise ValueError("consecutive design-orbit stations must differ")
        path_stations = self.path_length_stations
        if path_stations is not None:
            path_stations = _finite_array(
                path_stations,
                shape=(len(positions),),
                name="path_length_stations",
            )
            path_stations = path_stations - path_stations[0]
            if np.any(np.diff(path_stations) <= 0.0):
                raise ValueError(
                    "path_length_stations must be strictly increasing"
                )
        prescribed_curvature = self.signed_curvature_per_m
        if prescribed_curvature is not None:
            prescribed_curvature = _finite_array(
                prescribed_curvature,
                shape=(len(positions) - 1,),
                name="signed_curvature_per_m",
            )
        planar_error = np.max(np.abs((positions - positions[0]) @ axis))
        tangent_error = np.max(np.abs(tangents @ axis))
        if planar_error > 1.0e-9 * scale or tangent_error > 1.0e-9:
            raise ValueError(
                "design orbit and tangents must lie in the plane normal to "
                "bend_axis")
        chord = np.diff(positions, axis=0) / chord_lengths[:, None]
        alignment = np.einsum(
            "ij,ij->i", chord, tangents[:-1] + tangents[1:])
        if np.any(alignment <= 0.0):
            raise ValueError(
                "design-orbit tangents must point from entrance to exit")
        object.__setattr__(self, "positions", positions.copy())
        object.__setattr__(self, "tangents", tangents.copy())
        object.__setattr__(self, "bend_axis", axis.copy())
        object.__setattr__(self, "magnetic_rigidity", rigidity)
        object.__setattr__(
            self,
            "path_length_stations",
            None if path_stations is None else path_stations.copy(),
        )
        object.__setattr__(
            self,
            "signed_curvature_per_m",
            (
                None
                if prescribed_curvature is None
                else prescribed_curvature.copy()
            ),
        )

    @property
    def segment_lengths(self) -> np.ndarray:
        if self.path_length_stations is not None:
            return np.diff(self.path_length_stations)
        return self.chord_lengths

    @property
    def arc_length_stations(self) -> np.ndarray:
        """Return the design-orbit station coordinates in metres."""
        if self.path_length_stations is not None:
            return self.path_length_stations.copy()
        return np.r_[0.0, np.cumsum(self.chord_lengths)]

    @property
    def length_m(self) -> float:
        """Return the represented design-orbit length."""
        return float(self.arc_length_stations[-1])

    @property
    def chord_lengths(self) -> np.ndarray:
        return np.linalg.norm(np.diff(self.positions, axis=0), axis=1)

    @property
    def sample_positions(self) -> np.ndarray:
        # Cubic-Hermite midpoint stays on the tracked curve much more closely
        # than a chord midpoint, while remaining exact on straight segments.
        return (
            0.5 * (self.positions[:-1] + self.positions[1:])
            + self.segment_lengths[:, None]
            * (self.tangents[:-1] - self.tangents[1:])
            / 8.0
        )

    @property
    def geometric_signed_curvature(self) -> np.ndarray:
        """Return tangent-turning curvature for each orbit segment."""
        left = self.tangents[:-1]
        right = self.tangents[1:]
        sine = np.einsum(
            "j,ij->i", self.bend_axis, np.cross(left, right))
        cosine = np.einsum("ij,ij->i", left, right)
        turning = np.arctan2(sine, cosine)
        return turning / self.segment_lengths

    @property
    def signed_curvature(self) -> np.ndarray:
        """Return the collocated design curvature ``h(s)`` per segment.

        A tracked orbit may supply the ODE/B-field curvature at each segment
        collocation point.  Otherwise the tangent-turning average remains the
        backward-compatible source.  Keeping both values exposes finite-
        station geometry/collocation error instead of silently folding it into
        the canonical Hamiltonian's linear term.
        """
        if self.signed_curvature_per_m is not None:
            return self.signed_curvature_per_m.copy()
        return self.geometric_signed_curvature

    def _segment_coordinates(self, s_m):
        raw = np.asarray(s_m)
        if np.iscomplexobj(raw):
            raise ValueError("s_m must contain finite real arc lengths")
        query = np.asarray(raw, dtype=float)
        if not np.all(np.isfinite(query)):
            raise ValueError("s_m must contain finite real arc lengths")
        stations = self.arc_length_stations
        tolerance = 64.0 * np.finfo(float).eps * max(1.0, self.length_m)
        if np.any(query < -tolerance) or np.any(query > self.length_m + tolerance):
            raise ValueError("s_m is outside the represented design orbit")
        clipped = np.clip(query, 0.0, self.length_m)
        segment = np.searchsorted(stations, clipped, side="right") - 1
        segment = np.clip(segment, 0, len(stations) - 2)
        local = (clipped - stations[segment]) / self.segment_lengths[segment]
        return query.ndim == 0, segment, local

    def position_at(self, s_m) -> np.ndarray:
        """Evaluate global ``(X(s), Y(s), Z(s))`` by cubic Hermite interpolation."""
        scalar, segment, local = self._segment_coordinates(s_m)
        u = np.asarray(local, dtype=float).reshape(-1)
        indices = np.asarray(segment, dtype=np.int64).reshape(-1)
        ds = self.segment_lengths[indices]
        p0 = self.positions[indices]
        p1 = self.positions[indices + 1]
        t0 = self.tangents[indices]
        t1 = self.tangents[indices + 1]
        h00 = 2.0 * u**3 - 3.0 * u**2 + 1.0
        h10 = u**3 - 2.0 * u**2 + u
        h01 = -2.0 * u**3 + 3.0 * u**2
        h11 = u**3 - u**2
        result = (
            h00[:, None] * p0
            + (h10 * ds)[:, None] * t0
            + h01[:, None] * p1
            + (h11 * ds)[:, None] * t1
        )
        return result[0] if scalar else result.reshape(np.shape(s_m) + (3,))

    def tangent_at(self, s_m) -> np.ndarray:
        """Evaluate the unit design-orbit tangent at arc length ``s_m``."""
        scalar, segment, local = self._segment_coordinates(s_m)
        u = np.asarray(local, dtype=float).reshape(-1)
        indices = np.asarray(segment, dtype=np.int64).reshape(-1)
        ds = self.segment_lengths[indices]
        p0 = self.positions[indices]
        p1 = self.positions[indices + 1]
        t0 = self.tangents[indices]
        t1 = self.tangents[indices + 1]
        derivative = (
            ((6.0 * u**2 - 6.0 * u) / ds)[:, None] * p0
            + (3.0 * u**2 - 4.0 * u + 1.0)[:, None] * t0
            + ((-6.0 * u**2 + 6.0 * u) / ds)[:, None] * p1
            + (3.0 * u**2 - 2.0 * u)[:, None] * t1
        )
        derivative /= np.linalg.norm(derivative, axis=1)[:, None]
        return derivative[0] if scalar else derivative.reshape(np.shape(s_m) + (3,))

    def signed_curvature_at(self, s_m) -> np.ndarray | float:
        """Evaluate the segment-averaged design curvature ``h(s)``."""
        scalar, segment, _ = self._segment_coordinates(s_m)
        result = self.signed_curvature[np.asarray(segment, dtype=np.int64)]
        return float(result) if scalar else np.asarray(result, dtype=float)

    def frame_at(self, s_m) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return local ``(horizontal, vertical, tangent)`` axes in global xyz."""
        tangent = self.tangent_at(s_m)
        horizontal = np.cross(self.bend_axis, tangent)
        horizontal /= np.linalg.norm(horizontal, axis=-1, keepdims=True)
        vertical = np.broadcast_to(self.bend_axis, np.shape(horizontal)).copy()
        return horizontal, vertical, tangent

    def local_to_global(self, s_m, x_m=0.0, y_m=0.0) -> np.ndarray:
        """Map planar moving-frame coordinates ``(s,x,y)`` to global xyz."""
        position = self.position_at(s_m)
        horizontal, vertical, _ = self.frame_at(s_m)
        x = np.asarray(x_m, dtype=float)
        y = np.asarray(y_m, dtype=float)
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
            raise ValueError("x_m and y_m must be finite")
        return position + x[..., None] * horizontal + y[..., None] * vertical


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


@dataclass(frozen=True)
class MeasuredMedianPlaneFieldTarget:
    """Measured local B samples used only as an HDiv-MMM inverse target.

    The data contain the local components named by ``components`` at physical
    ``(s,x,y=0)`` points.  A full ``(B_x,B_y,B_s)`` last axis is also accepted
    and reduced to the selected components.  The values are never interpolated
    into an off-plane field and never converted into a vector potential.
    Instead, HDiv-MMM observation rows compare a three-dimensional magnet solve
    with these samples while the optimizer changes the pole topology.
    """

    orbit: PlanarDesignOrbit
    s_m: np.ndarray
    x_m: np.ndarray
    measured_B_local_t: np.ndarray
    measurement_band_t: np.ndarray | float
    components: tuple[str, ...] = ("x", "y", "s")

    def __post_init__(self):
        if not isinstance(self.orbit, PlanarDesignOrbit):
            raise TypeError("orbit must be a PlanarDesignOrbit")
        s = _finite_array(self.s_m, name="s_m").reshape(-1)
        x = _finite_array(self.x_m, name="x_m").reshape(-1)
        if s.size == 0 or x.size == 0 or np.any(np.diff(s) <= 0.0):
            raise ValueError("s_m must increase strictly and x_m must not be empty")
        components = tuple(str(value).lower() for value in self.components)
        allowed = {"x", "y", "s"}
        if (
            not components
            or len(set(components)) != len(components)
            or any(value not in allowed for value in components)
        ):
            raise ValueError(
                "components must be a non-empty unique subset of ('x','y','s')"
            )
        indices = tuple({"x": 0, "y": 1, "s": 2}[value] for value in components)
        measured_raw = _finite_array(
            self.measured_B_local_t, name="measured_B_local_t"
        )
        selected_shape = (len(s), len(x), len(components))
        full_shape = (len(s), len(x), 3)
        if measured_raw.shape == full_shape:
            measured = measured_raw[:, :, indices]
        elif measured_raw.shape == selected_shape:
            measured = measured_raw
        else:
            raise ValueError(
                "measured_B_local_t must have shape (n_s,n_x,n_component) "
                "or (n_s,n_x,3)"
            )
        band_raw = np.asarray(self.measurement_band_t, dtype=float)
        try:
            band = np.broadcast_to(band_raw, measured.shape).copy()
        except ValueError:
            try:
                band = np.broadcast_to(band_raw, full_shape)[:, :, indices].copy()
            except ValueError as exc:
                raise ValueError(
                    "measurement_band_t must broadcast to the selected or "
                    "full local-component shape"
                ) from exc
        if not np.all(np.isfinite(band)) or np.any(band <= 0.0):
            raise ValueError("measurement_band_t must be finite and positive")
        # This also rejects stations outside the orbit interpolation domain.
        self.orbit.position_at(s)
        object.__setattr__(self, "s_m", np.ascontiguousarray(s))
        object.__setattr__(self, "x_m", np.ascontiguousarray(x))
        object.__setattr__(self, "measured_B_local_t", measured.copy())
        object.__setattr__(self, "measurement_band_t", band)
        object.__setattr__(self, "components", components)

    @property
    def component_indices(self) -> tuple[int, ...]:
        lookup = {"x": 0, "y": 1, "s": 2}
        return tuple(lookup[value] for value in self.components)

    @property
    def observation_points_m(self) -> np.ndarray:
        s, x = np.meshgrid(self.s_m, self.x_m, indexing="ij")
        return np.ascontiguousarray(
            self.orbit.local_to_global(s.reshape(-1), x.reshape(-1), 0.0)
        )

    @property
    def observation_weights(self) -> np.ndarray:
        horizontal, vertical, tangent = self.orbit.frame_at(self.s_m)
        local_basis = np.stack((horizontal, vertical, tangent), axis=1)
        point_basis = np.repeat(local_basis, len(self.x_m), axis=0)
        indices = self.component_indices
        point_count = len(point_basis)
        weights = np.zeros((point_count * len(indices), point_count, 3))
        for point in range(point_count):
            for local, component in enumerate(indices):
                weights[point * len(indices) + local, point] = point_basis[
                    point, component
                ]
        return np.ascontiguousarray(weights)

    @property
    def response_target(self) -> np.ndarray:
        return np.ascontiguousarray(self.measured_B_local_t.reshape(-1))

    @property
    def response_band(self) -> np.ndarray:
        return np.ascontiguousarray(self.measurement_band_t.reshape(-1))

    @property
    def raw_field_response_size(self) -> int:
        return int(self.response_target.size)


def build_measured_median_plane_field_response_matrix(
    charge_gram,
    target: MeasuredMedianPlaneFieldTarget,
    *,
    field_scale=MU0,
) -> np.ndarray:
    """Build native HDiv-MMM rows at the actual measurement locations."""
    if not isinstance(target, MeasuredMedianPlaneFieldTarget):
        raise TypeError("target must be a MeasuredMedianPlaneFieldTarget")
    native = getattr(charge_gram, "configured_field_functional_rows", None)
    if native is None:
        raise TypeError(
            "charge_gram must expose configured_field_functional_rows"
        )
    scale = float(field_scale)
    if not np.isfinite(scale) or scale == 0.0:
        raise ValueError("field_scale must be finite and nonzero")
    rows = scale * np.asarray(
        native(target.observation_points_m, target.observation_weights),
        dtype=float,
    )
    if (
        rows.ndim != 2
        or rows.shape[0] != target.raw_field_response_size
        or not np.all(np.isfinite(rows))
    ):
        raise RuntimeError(
            "native configured-field API returned invalid measurement rows"
        )
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
    derivative_backend: str = "forward-mode-expm-frechet-ad"
    field_to_design_jacobian: np.ndarray | None = None


@dataclass(frozen=True)
class TransferMatrixAutomaticDifferentiation:
    """Forward-mode derivative from sampled magnetic field to optics rows.

    ``field_basis`` is the seed matrix.  The full Jacobian differentiates the
    bend-field rows and selected transfer-matrix entries with respect to every
    raw ``[B..., dB/dx...]`` sample.  ``directional_jacobian`` is the same
    derivative after applying an optional smooth-field basis.
    """

    backend: str
    current_field_response: np.ndarray
    field_basis: np.ndarray
    design_response: np.ndarray
    full_jacobian: np.ndarray
    directional_jacobian: np.ndarray


@dataclass(frozen=True)
class TransferMatrixMaterialInversePipelineResult:
    """Auditable field -> map -> field error -> material proposal chain."""

    objective: object
    field_distribution: np.ndarray
    realized_transfer_matrices: np.ndarray
    target_transfer_matrices: np.ndarray
    transfer_matrix_difference: np.ndarray
    normalized_transfer_matrix_difference: np.ndarray
    automatic_differentiation: TransferMatrixAutomaticDifferentiation
    field_correction: TransferMatrixFieldCorrection
    material_selection: TSVDElementCandidateSelection
    stage_order: tuple[str, ...] = (
        "magnetic-field-distribution",
        "forward-ad-transfer-matrix",
        "target-transfer-matrix-difference",
        "tsvd-minimax-field-correction",
        "aca-thin-qr-tsvd-material-inverse",
    )

    @property
    def proposed_field_distribution(self) -> np.ndarray:
        return np.asarray(
            self.material_selection.predicted_response, dtype=float)

    @property
    def status(self) -> str:
        return self.material_selection.status


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
class MeasuredMedianPlaneTopologyResult:
    """Pole-topology result re-solved against measured median-plane B."""

    target: MeasuredMedianPlaneFieldTarget
    generation: HDivMMMGenerationResult
    realized_field_response_t: np.ndarray
    maximum_measurement_band_ratio: float
    topology: GrowthTopologyReport

    @property
    def active_elements(self) -> np.ndarray:
        return self.generation.active_elements

    @property
    def converged(self) -> bool:
        return bool(
            self.generation.converged
            and self.maximum_measurement_band_ratio <= 1.0
        )


@dataclass(frozen=True)
class MultiMomentumTransferMatrixObjective:
    """Joint orbit/map contract at several magnetic rigidities.

    Every orbit owns its physical observation points and may contain a
    different number of longitudinal segments.  Raw rows are concatenated in
    momentum order, with each orbit retaining the row order
    ``[B_binormal..., dB_binormal/dnormal...]``.  The transformed response is
    the corresponding concatenation of
    :class:`PlanarTransferMatrixObjective` responses.  Its Jacobian is block
    diagonal and uses the same forward-mode matrix-exponential AD chain; no
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
                     for left, right in itertools.pairwise(offsets))

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


def differentiate_transfer_matrix_field_response(
        objective, current_field_response, *, field_basis=None
        ) -> TransferMatrixAutomaticDifferentiation:
    """Differentiate sampled field -> selected optics response by AD.

    The objective seeds every raw field coordinate and propagates those
    tangents through its declared map.  An optional ``field_basis`` is
    composed with that full Jacobian afterward.  The first-order combined-
    function objective uses forward scalar rules and a Fréchet-differentiated
    matrix exponential; higher-order objectives may declare another exact AD
    backend through ``derivative_backend``.  This is algorithmic
    differentiation of the optics model, not a design finite difference.
    """
    required = (
        "raw_field_response_size", "transform", "transform_jacobian")
    if any(not hasattr(objective, name) for name in required):
        raise TypeError(
            "objective must expose the planar transfer-matrix field contract")
    raw_size = int(objective.raw_field_response_size)
    current = _finite_array(
        current_field_response, name="current_field_response").reshape(-1)
    if current.shape != (raw_size,):
        raise ValueError(
            "current_field_response must match the objective raw-field size")
    if field_basis is None:
        basis = np.eye(raw_size)
    else:
        basis = _finite_array(field_basis, name="field_basis")
        if basis.ndim != 2 or basis.shape[0] != raw_size or basis.shape[1] == 0:
            raise ValueError(
                "field_basis must have shape "
                "(raw_field_response_size,n_mode)")
    design = _finite_array(
        objective.transform(current), name="current_design_response"
    ).reshape(-1)
    jacobian = _finite_array(
        objective.transform_jacobian(current),
        name="field_to_transfer_jacobian")
    if jacobian.shape != (design.size, raw_size):
        raise RuntimeError(
            "transfer objective returned an incompatible AD Jacobian")
    backend = str(getattr(
        objective, "derivative_backend",
        "forward-mode-expm-frechet-ad"))
    return TransferMatrixAutomaticDifferentiation(
        backend=backend,
        current_field_response=current.copy(),
        field_basis=np.asarray(basis, dtype=float).copy(),
        design_response=design.copy(),
        full_jacobian=np.asarray(jacobian, dtype=float).copy(),
        directional_jacobian=np.asarray(
            jacobian @ basis, dtype=float))


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
    automatic = differentiate_transfer_matrix_field_response(
        objective, current, field_basis=field_basis)
    basis = automatic.field_basis
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

    current_design = automatic.design_response
    jacobian = automatic.full_jacobian
    if (current_design.shape != target.shape
            or jacobian.shape != (target.size, raw_size)):
        raise RuntimeError(
            "transfer objective returned an incompatible response or Jacobian")
    normalized_residual = (target - current_design) / band
    normalized_operator = automatic.directional_jacobian / band[:, None]
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
        status=status,
        derivative_backend=automatic.backend,
        field_to_design_jacobian=automatic.full_jacobian.copy())


def run_transfer_matrix_material_inverse_pipeline(
        objective, current_field_response, *, candidate_elements,
        candidate_field_response_delta, candidate_volumes, volume_budget,
        field_basis=None, field_inverse_relative_tolerance=1.0e-3,
        field_inverse_maximum_step_scale=1.0,
        field_inverse_line_search_steps=8,
        material_relative_tolerance=1.0e-3,
        material_improvement_capture=0.9,
        ratio_tolerance=1.0e-12, active_elements=None,
        predecessor_elements=None, candidate_volume_changes=None,
        candidate_material_active=None, candidate_exclusion_groups=None,
        maximum_changed_volume=None, maximum_changed_elements=None,
        candidate_secondary_cost=None
        ) -> TransferMatrixMaterialInversePipelineResult:
    """Run the complete field-to-binary-material inverse pipeline.

    The five explicit stages are:

    1. accept the sampled ``[B..., dB/dx...]`` magnetic-field distribution;
    2. build the realized transfer matrix and its forward-mode AD Jacobian;
    3. form the target-minus-realized transfer-matrix difference;
    4. solve a band-normalized optics TSVD/minimax inverse for ``delta B``;
    5. approximate all candidate ``delta B`` columns with ACA -> thin QR ->
       TSVD and solve the whole-element binary material proposal.

    The returned proposal is a screening result.  A production topology loop
    must still evaluate its exact Schur block, completely re-solve the active
    HDiv-MMM system, rebuild the transfer matrix, and accept only an improving
    exact result.
    """
    if not isinstance(
            objective,
            (PlanarTransferMatrixObjective,
             MultiMomentumTransferMatrixObjective)):
        raise TypeError(
            "objective must be PlanarTransferMatrixObjective or "
            "MultiMomentumTransferMatrixObjective")
    current = _finite_array(
        current_field_response, name="current_field_response").reshape(-1)
    automatic = differentiate_transfer_matrix_field_response(
        objective, current, field_basis=field_basis)
    correction = solve_transfer_matrix_field_correction(
        objective, current, field_basis=automatic.field_basis,
        relative_tolerance=field_inverse_relative_tolerance,
        maximum_step_scale=field_inverse_maximum_step_scale,
        line_search_steps=field_inverse_line_search_steps)

    if isinstance(objective, PlanarTransferMatrixObjective):
        raw_by_objective = (current,)
        objectives = (objective,)
        targets = objective.target_matrix[None, :, :]
        matrix_bands = objective.transfer_matrix_band[None, :, :]
    else:
        raw_by_objective = objective.split_raw_response(current)
        objectives = objective.objectives
        targets = objective.target_matrices
        matrix_bands = objective.transfer_matrix_band
    realized = np.asarray([
        item.evaluate_transfer_map(raw).matrix
        for item, raw in zip(objectives, raw_by_objective)], dtype=float)
    difference = np.asarray(targets - realized, dtype=float)
    normalized_difference = difference / matrix_bands

    selection = select_tsvd_element_candidates(
        current_response=current,
        response_target=correction.target_field_response,
        response_band=correction.field_response_band,
        candidate_elements=candidate_elements,
        candidate_response_delta=candidate_field_response_delta,
        candidate_volumes=candidate_volumes,
        volume_budget=volume_budget,
        active_elements=active_elements,
        predecessor_elements=predecessor_elements,
        relative_tolerance=material_relative_tolerance,
        improvement_capture=material_improvement_capture,
        ratio_tolerance=ratio_tolerance,
        candidate_volume_changes=candidate_volume_changes,
        candidate_material_active=candidate_material_active,
        candidate_exclusion_groups=candidate_exclusion_groups,
        maximum_changed_volume=maximum_changed_volume,
        maximum_changed_elements=maximum_changed_elements,
        candidate_secondary_cost=candidate_secondary_cost)
    return TransferMatrixMaterialInversePipelineResult(
        objective=objective,
        field_distribution=current.copy(),
        realized_transfer_matrices=realized,
        target_transfer_matrices=np.asarray(targets, dtype=float).copy(),
        transfer_matrix_difference=difference,
        normalized_transfer_matrix_difference=normalized_difference,
        automatic_differentiation=automatic,
        field_correction=correction,
        material_selection=selection)


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


def multi_orbit_field_observations(
        objective: MultiMomentumTransferMatrixObjective, *,
        gradient_offset, field_scale=MU0) -> tuple[np.ndarray, np.ndarray]:
    """Materialize the common native observation contract for all orbits.

    The returned points have shape ``(n_point,3)`` and the row-major vector
    weights have shape ``(objective.raw_field_response_size,n_point,3)``.
    Unlike :func:`build_multi_orbit_field_response_matrix`, this deliberately
    constructs the zero-padded block tensor required by the native analytic
    configured-field *shape derivative*.  Ordinary field-row assembly should
    keep using the batched builder above to avoid this quadratic workspace.
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
    scale = float(field_scale)
    if not np.isfinite(scale) or scale == 0.0:
        raise ValueError("field_scale must be finite and nonzero")
    batches = tuple(
        planar_orbit_field_observations(orbit, gradient_offset=offset)
        for orbit, offset in zip(objective.orbits, offsets))
    points = np.ascontiguousarray(np.vstack([batch[0] for batch in batches]))
    weights = np.zeros(
        (objective.raw_field_response_size, len(points), 3), dtype=float)
    raw_offsets = objective.raw_offsets
    point_offsets = np.r_[0, np.cumsum([len(batch[0]) for batch in batches])]
    for index, (_, block) in enumerate(batches):
        weights[
            raw_offsets[index]:raw_offsets[index + 1],
            point_offsets[index]:point_offsets[index + 1], :
        ] = scale * block
    return points, np.ascontiguousarray(weights)


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

    def to_radia_object(self, *, closure_tolerance=1.0e-9) -> int:
        """Materialize the same closed filaments as one native Radia object.

        ``CoilBuilderHDivSource`` owns the finite-filament representation used
        for both HDiv RHS assembly and incident-field sampling.  Native orbit
        tracking cannot call the Python Biot--Savart evaluator, so it receives
        an ``ObjFlmCur`` container built from those exact segments rather than
        a separately discretized solid-current coil.

        The conversion deliberately accepts only continuous closed paths.
        ``from_coilbuilders`` already establishes this contract; the explicit
        check keeps manually constructed sources from silently acquiring an
        unphysical closing segment at the pybind boundary.
        """
        import radia as rad

        tolerance = float(closure_tolerance)
        if not np.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError("closure_tolerance must be finite and nonnegative")
        objects = []
        for index, (segments, current) in enumerate(self.segment_groups):
            joins = np.linalg.norm(segments[1:, 0] - segments[:-1, 1], axis=1)
            maximum_join = float(np.max(joins, initial=0.0))
            closure_gap = float(np.linalg.norm(segments[-1, 1] - segments[0, 0]))
            if maximum_join > tolerance or closure_gap > tolerance:
                raise ValueError(
                    "native Radia coil materialization requires a continuous "
                    f"closed filament path (group {index}: max join "
                    f"{maximum_join:.6e} m, closure {closure_gap:.6e} m)"
                )
            points = np.vstack((segments[:, 0], segments[0, 0]))
            objects.append(rad.ObjFlmCur(points.tolist(), float(current)))
        return int(rad.ObjCnt(objects))

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


def optimize_hdiv_mmm_magnet_to_measured_median_plane(
    target: MeasuredMedianPlaneFieldTarget,
    *,
    charge_gram,
    fes,
    inv_chi,
    rhs,
    field_response_matrix,
    active_elements,
    element_volumes,
    volume_max,
    incident_field_response=None,
    **generation_options,
) -> MeasuredMedianPlaneTopologyResult:
    """Adjust whole pole elements so the 3-D solution matches measured B.

    The measured samples are used only in the band-normalized topology
    objective.  Every accepted candidate is a complete HDiv-MMM re-solve;
    this routine does not form an interpolated field, an off-plane B-spline,
    or an A-map from the measurements.  The accepted 3-D magnet must be solved
    separately for its HCurl A-map and independent direct B-map before Lie/RK
    validation.
    """
    if not isinstance(target, MeasuredMedianPlaneFieldTarget):
        raise TypeError("target must be a MeasuredMedianPlaneFieldTarget")
    response_matrix = _finite_array(
        field_response_matrix, name="field_response_matrix"
    )
    expected_shape = (target.raw_field_response_size, int(fes.ndof))
    if response_matrix.ndim != 2 or response_matrix.shape != expected_shape:
        raise ValueError(
            f"field_response_matrix must have shape {expected_shape}"
        )
    incident = (
        np.zeros(target.raw_field_response_size)
        if incident_field_response is None
        else _finite_array(
            incident_field_response, name="incident_field_response"
        ).reshape(-1)
    )
    if incident.shape != (target.raw_field_response_size,):
        raise ValueError(
            "incident_field_response must match the measurement response size"
        )
    reserved = {
        "response_matrix",
        "response_target",
        "response_band",
        "response_transform",
        "response_transform_jacobian",
        "incident_response",
    }
    overlap = reserved.intersection(generation_options)
    if overlap:
        raise TypeError(
            "generation_options cannot override the measured-field contract: "
            + ", ".join(sorted(overlap))
        )
    generation = grow_hdiv_mmm_by_superposition(
        charge_gram=charge_gram,
        fes=fes,
        inv_chi=inv_chi,
        rhs=rhs,
        response_matrix=response_matrix,
        active_elements=active_elements,
        element_volumes=element_volumes,
        response_target=target.response_target,
        response_band=target.response_band,
        volume_max=volume_max,
        incident_response=incident,
        **generation_options,
    )
    realized = response_matrix @ generation.state + incident
    ratio = float(
        np.max(
            np.abs((realized - target.response_target) / target.response_band),
            initial=0.0,
        )
    )
    return MeasuredMedianPlaneTopologyResult(
        target=target,
        generation=generation,
        realized_field_response_t=np.ascontiguousarray(realized),
        maximum_measurement_band_ratio=ratio,
        topology=ngsolve_growth_topology(fes.mesh, generation.active_elements),
    )


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
    optics transform and its forward-mode AD Jacobian are block-assembled
    before the existing batched adjoint contraction and ACA--QR--TSVD master
    step.
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
    "MeasuredMedianPlaneFieldTarget",
    "MeasuredMedianPlaneTopologyResult",
    "MultiMomentumAcceleratorMagnetTopologyResult",
    "MultiMomentumTransferMatrixObjective",
    "PlanarDesignOrbit",
    "PlanarTransferMatrixObjective",
    "TransferMatrixAutomaticDifferentiation",
    "TransferMatrixFieldCorrection",
    "TransferMatrixMaterialInversePipelineResult",
    "build_measured_median_plane_field_response_matrix",
    "build_multi_orbit_field_response_matrix",
    "build_planar_orbit_field_response_matrix",
    "differentiate_transfer_matrix_field_response",
    "optimize_hdiv_mmm_magnet_from_transfer_matrices",
    "optimize_hdiv_mmm_magnet_from_transfer_matrix",
    "optimize_hdiv_mmm_magnet_to_measured_median_plane",
    "multi_orbit_field_observations",
    "planar_orbit_field_observations",
    "run_transfer_matrix_material_inverse_pipeline",
    "solve_transfer_matrix_field_correction",
    "static_magnet_transfer_component_entries",
    "static_magnet_symplectic_residual",
]
