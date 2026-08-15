"""Curvilinear beam-tube meshes on a Bishop rotation-minimizing frame.

The geometric frame uses the Wang--Juttler--Zheng--Liu double-reflection
discretization used by EarlyTimes' native NGSolve field adapter.  It is a
Frenet--Serret-compatible frame on curved planar trajectories, but remains
defined at straight and zero-curvature stations where the Frenet normal does
not.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise
from typing import Any

import numpy as np


def _finite_array(value, name: str, *, shape=None) -> np.ndarray:
    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise TypeError(f"{name} must be real")
    array = np.ascontiguousarray(raw, dtype=np.float64)
    if shape is not None and array.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _normalize(vector: np.ndarray, name: str) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if not np.isfinite(norm) or norm <= 64.0 * np.finfo(float).eps:
        raise ValueError(f"{name} must have nonzero finite norm")
    return vector / norm


def _project_normal(
    vector: np.ndarray, tangent: np.ndarray, name: str
) -> np.ndarray:
    return _normalize(vector - np.dot(vector, tangent) * tangent, name)


def _seed_horizontal(initial: np.ndarray, tangent: np.ndarray) -> np.ndarray:
    projected = initial - np.dot(initial, tangent) * tangent
    if np.linalg.norm(projected) > 64.0 * np.finfo(float).eps:
        return _normalize(projected, "initial horizontal axis")
    axes = np.eye(3)
    least_aligned = axes[np.argmin(np.abs(axes @ tangent))]
    return _project_normal(
        least_aligned, tangent, "fallback horizontal axis"
    )


def _reflect(
    vector: np.ndarray, plane_normal: np.ndarray, normal_squared: float
) -> np.ndarray:
    return vector - 2.0 * np.dot(plane_normal, vector) / normal_squared * plane_normal


def _transport_double_reflection(
    previous_position: np.ndarray,
    position: np.ndarray,
    previous_tangent: np.ndarray,
    tangent: np.ndarray,
    previous_horizontal: np.ndarray,
) -> np.ndarray:
    chord = position - previous_position
    chord_squared = float(np.dot(chord, chord))
    position_scale = max(
        1.0, float(np.linalg.norm(previous_position)), float(np.linalg.norm(position))
    )
    tolerance = 64.0 * np.finfo(float).eps * position_scale
    if chord_squared <= tolerance * tolerance:
        raise ValueError(
            "consecutive reference positions must be distinct for the "
            "Bishop double-reflection frame"
        )

    reflected_horizontal = _reflect(previous_horizontal, chord, chord_squared)
    reflected_tangent = _reflect(previous_tangent, chord, chord_squared)
    tangent_difference = tangent - reflected_tangent
    difference_squared = float(np.dot(tangent_difference, tangent_difference))
    horizontal = reflected_horizontal
    tangent_tolerance = 64.0 * np.finfo(float).eps
    if difference_squared > tangent_tolerance * tangent_tolerance:
        horizontal = _reflect(
            reflected_horizontal, tangent_difference, difference_squared
        )
    return _project_normal(
        horizontal, tangent, "Bishop transported horizontal axis"
    )


def _rotate_about_tangent(
    normal: np.ndarray, tangent: np.ndarray, angle: float
) -> np.ndarray:
    rotated = (
        np.cos(angle) * normal
        + np.sin(angle) * np.cross(tangent, normal)
    )
    return _project_normal(rotated, tangent, "minimal-twist horizontal axis")


@dataclass(frozen=True)
class BishopRMFFrame:
    """Station-wise right-handed ``(x, y, s)`` moving frame.

    ``horizontal`` is the local ``x`` direction, ``vertical`` is local ``y``,
    and ``tangent`` is local ``s``.  Thus ``tangent x horizontal = vertical``.
    ``arc_length_m`` is the cumulative chord length from the first station.
    """

    positions_m: np.ndarray
    tangent: np.ndarray
    horizontal: np.ndarray
    vertical: np.ndarray
    arc_length_m: np.ndarray
    periodic: bool = False
    holonomy_correction_rad: float = 0.0

    def local_to_global(self, station: int, x_m=0.0, y_m=0.0) -> np.ndarray:
        """Map local transverse coordinates at one station to global xyz."""
        index = int(station)
        if index < 0 or index >= len(self.positions_m):
            raise IndexError("station is outside the moving frame")
        x = float(x_m)
        y = float(y_m)
        if not np.all(np.isfinite([x, y])):
            raise ValueError("x_m and y_m must be finite")
        return (
            self.positions_m[index]
            + x * self.horizontal[index]
            + y * self.vertical[index]
        )

    def vector_components(self, station: int, vector) -> np.ndarray:
        """Return global-vector components in ``(x, y, s)`` order."""
        index = int(station)
        if index < 0 or index >= len(self.positions_m):
            raise IndexError("station is outside the moving frame")
        value = _finite_array(vector, "vector", shape=(3,))
        return np.array(
            [
                np.dot(value, self.horizontal[index]),
                np.dot(value, self.vertical[index]),
                np.dot(value, self.tangent[index]),
            ]
        )


def bishop_rmf_frame(
    positions_m,
    tangents,
    *,
    initial_horizontal=(1.0, 0.0, 0.0),
    periodic: bool = False,
) -> BishopRMFFrame:
    """Build the same double-reflection Bishop/RMF used by EarlyTimes.

    The input stations can contain straight-to-curved transitions and
    inflection/zero-curvature points.  A periodic frame expects distinct first
    and last stations on a sampled closed loop; it distributes the one-turn
    holonomy uniformly in chord arc length.
    """
    positions = _finite_array(positions_m, "positions_m")
    tangent_values = _finite_array(tangents, "tangents")
    if (
        positions.ndim != 2
        or positions.shape[1] != 3
        or positions.shape[0] < 2
        or tangent_values.shape != positions.shape
    ):
        raise ValueError(
            "positions_m and tangents need matching shape (n_station, 3) "
            "with at least two stations"
        )
    tangent_norm = np.linalg.norm(tangent_values, axis=1)
    if np.any(tangent_norm <= 64.0 * np.finfo(float).eps):
        raise ValueError("tangents must have nonzero finite norm")
    tangent_values = tangent_values / tangent_norm[:, None]
    seed = _finite_array(initial_horizontal, "initial_horizontal", shape=(3,))

    horizontal = np.empty_like(tangent_values)
    horizontal[0] = _seed_horizontal(seed, tangent_values[0])
    for index in range(1, len(positions)):
        horizontal[index] = _transport_double_reflection(
            positions[index - 1],
            positions[index],
            tangent_values[index - 1],
            tangent_values[index],
            horizontal[index - 1],
        )

    chords = np.linalg.norm(np.diff(positions, axis=0), axis=1)
    arc_length = np.concatenate(([0.0], np.cumsum(chords)))
    correction = 0.0
    if periodic:
        if len(positions) < 3:
            raise ValueError(
                "periodic Bishop frame needs at least three reference stations"
            )
        closure = _transport_double_reflection(
            positions[-1],
            positions[0],
            tangent_values[-1],
            tangent_values[0],
            horizontal[-1],
        )
        correction = float(
            np.arctan2(
                np.dot(
                    tangent_values[0], np.cross(closure, horizontal[0])
                ),
                np.dot(closure, horizontal[0]),
            )
        )
        total_length = float(arc_length[-1] + np.linalg.norm(positions[0] - positions[-1]))
        if not np.isfinite(total_length) or total_length <= 0.0:
            raise ValueError("periodic Bishop frame needs a nondegenerate closed path")
        for index in range(1, len(positions)):
            horizontal[index] = _rotate_about_tangent(
                horizontal[index],
                tangent_values[index],
                correction * arc_length[index] / total_length,
            )

    vertical = np.cross(tangent_values, horizontal)
    vertical /= np.linalg.norm(vertical, axis=1)[:, None]
    return BishopRMFFrame(
        positions_m=positions.copy(),
        tangent=tangent_values.copy(),
        horizontal=horizontal,
        vertical=vertical,
        arc_length_m=arc_length,
        periodic=bool(periodic),
        holonomy_correction_rad=correction,
    )


@dataclass(frozen=True)
class CurvilinearBeamMesh:
    """Structured four-HEX-per-cell tube on Bishop/RMF sections.

    The transverse topology is fixed to four quadrilateral x strips with
    ``y=0`` inside every strip and ``x=0`` on the central internal face.  The
    stored vertices/connectivity make that finite-element contract auditable
    without reverse-engineering an NGSolve mesh.
    """

    mesh: Any
    shape: Any
    frame: BishopRMFFrame
    half_width_m: float
    half_height_m: float
    maxh_m: float
    material: str
    boundary: str
    curve_order: int
    topology: str
    x_nodes_m: np.ndarray
    y_nodes_m: np.ndarray
    longitudinal_stations_m: np.ndarray
    vertices_m: np.ndarray
    hex_connectivity: np.ndarray

    @property
    def hexes_per_longitudinal_cell(self) -> int:
        return 4


@dataclass(frozen=True)
class OrbitGaugeVectorPotential:
    """HCurl vector potential in the ``A_s=A_y=0`` design-orbit gauge."""

    vector_potential: Any
    ungauged_vector_potential: Any
    gauge_potential: Any
    order: int
    frame: BishopRMFFrame
    As_before_t_m: np.ndarray
    Ay_before_t_m: np.ndarray
    As_after_t_m: np.ndarray
    Ay_after_t_m: np.ndarray
    rbf_width_m: float
    rbf_condition: float
    curl_change_l2_t_m32: float

    @property
    def maximum_orbit_gauge_residual_t_m(self) -> float:
        return float(
            max(
                np.max(np.abs(self.As_after_t_m), initial=0.0),
                np.max(np.abs(self.Ay_after_t_m), initial=0.0),
            )
        )


@dataclass(frozen=True)
class TransverseVectorPotentialData:
    """HCurl A-map sampled on a full local ``(x,y)`` tensor grid."""

    s_m: np.ndarray
    x_m: np.ndarray
    y_m: np.ndarray
    global_points_m: np.ndarray
    horizontal: np.ndarray
    vertical: np.ndarray
    tangent: np.ndarray
    Ax_t_m: np.ndarray
    Ay_t_m: np.ndarray
    As_t_m: np.ndarray
    grid_function_space_class: str
    grid_function_space_order: int | None

    @property
    def design_orbit_indices(self) -> tuple[int, int]:
        """Return the transverse sample indexes nearest ``x=y=0``."""
        return (
            int(np.argmin(np.abs(self.x_m))),
            int(np.argmin(np.abs(self.y_m))),
        )

    @property
    def maximum_design_orbit_residual_t_m(self) -> float:
        """Return max ``|Ax|,|Ay|,|As|`` on the sampled design orbit."""
        x_index, y_index = self.design_orbit_indices
        return float(
            max(
                np.max(np.abs(self.Ax_t_m[:, x_index, y_index]), initial=0.0),
                np.max(np.abs(self.Ay_t_m[:, x_index, y_index]), initial=0.0),
                np.max(np.abs(self.As_t_m[:, x_index, y_index]), initial=0.0),
            )
        )


@dataclass(frozen=True)
class EarlyTimesHCurlFieldCertificate:
    """Acceptance certificate for the constrained EarlyTimes HCurl A-map.

    The certificate is deliberately about the finite-element field, not about
    a subsequently fitted polynomial.  It records the three conditions used
    by the HCurl-specialized Lie-map boundary: transverse axial gauge
    ``A_x=0`` throughout the sampled aperture, ``A_y=A_s=0`` on the design
    orbit, and median-plane parity ``A_x,A_s`` even / ``A_y`` odd.
    """

    grid_function_space_class: str
    grid_function_space_order: int | None
    required_order: int
    symmetry_class: str
    maximum_Ax_t_m: float
    maximum_orbit_Ay_As_t_m: float
    maximum_Ax_even_parity_defect_t_m: float
    maximum_Ay_odd_parity_defect_t_m: float
    maximum_As_even_parity_defect_t_m: float
    maximum_Ay_even_parity_defect_t_m: float
    maximum_As_odd_parity_defect_t_m: float
    axial_gauge_tolerance_t_m: float
    orbit_gauge_tolerance_t_m: float
    symmetry_tolerance_t_m: float

    @property
    def maximum_symmetry_defect_t_m(self) -> float:
        """Return the largest forbidden component for the declared symmetry."""
        if self.symmetry_class == "normal":
            return float(
                max(
                    self.maximum_Ax_even_parity_defect_t_m,
                    self.maximum_Ay_odd_parity_defect_t_m,
                    self.maximum_As_even_parity_defect_t_m,
                )
            )
        if self.symmetry_class == "skew":
            return float(
                max(
                    self.maximum_Ax_even_parity_defect_t_m,
                    self.maximum_Ay_even_parity_defect_t_m,
                    self.maximum_As_odd_parity_defect_t_m,
                )
            )
        return 0.0


@dataclass(frozen=True)
class TransverseVectorPotentialPolynomialFit:
    """Internal full triangular jet recovered from a certified HCurl field."""

    data: TransverseVectorPotentialData
    field_certificate: EarlyTimesHCurlFieldCertificate
    degree: int
    symmetry_class: str
    canonical_powers: tuple[tuple[int, int], ...]
    Ay_coefficients_t_m: np.ndarray
    As_coefficients_t_m: np.ndarray
    Ay_sample_to_coefficients: np.ndarray
    As_sample_to_coefficients: np.ndarray
    maximum_Ax_t_m: float
    maximum_orbit_Ay_As_t_m: float
    maximum_Ay_fit_residual_t_m: float
    maximum_As_fit_residual_t_m: float
    maximum_left_right_scaled_coefficient_discrepancy_t_m: float
    scaled_Ay_design_condition: float
    scaled_As_design_condition: float
    scaled_left_Ay_design_condition: float
    scaled_right_Ay_design_condition: float
    scaled_left_As_design_condition: float
    scaled_right_As_design_condition: float

    @property
    def derivative_backend(self) -> str:
        """Name the exact linear derivative used for internal jet recovery."""
        return "scaled-least-squares-exact-linear-jacobian"

    def lie_parameter_response(
        self,
        Ay_sample_response,
        As_sample_response,
    ) -> np.ndarray:
        """Map sampled-HCurl responses to the internal Lie-jet ordering."""
        Ay_response = _finite_array(Ay_sample_response, "Ay_sample_response")
        As_response = _finite_array(As_sample_response, "As_sample_response")
        expected = (
            len(self.data.s_m),
            len(self.data.x_m),
            len(self.data.y_m),
        )
        if (
            Ay_response.ndim != 4
            or Ay_response.shape != As_response.shape
            or Ay_response.shape[:3] != expected
            or Ay_response.shape[3] < 1
        ):
            raise ValueError(
                "sample responses need matching shape "
                "(n_segment,n_x,n_y,n_mode)"
            )
        station_count = expected[0]
        mode_count = Ay_response.shape[3]
        Ay_flat = Ay_response.reshape(station_count, -1, mode_count)
        As_flat = As_response.reshape(station_count, -1, mode_count)
        Ay_coefficients = np.einsum(
            "pq,sqm->psm",
            self.Ay_sample_to_coefficients,
            Ay_flat,
            optimize=True,
        )
        As_coefficients = np.einsum(
            "pq,sqm->psm",
            self.As_sample_to_coefficients,
            As_flat,
            optimize=True,
        )
        return np.ascontiguousarray(
            np.concatenate((Ay_coefficients, As_coefficients), axis=0).reshape(
                -1, mode_count
            )
        )


def certify_earlytimes_hcurl_vector_potential(
    data: TransverseVectorPotentialData,
    *,
    required_order: int = 5,
    symmetry_class: str = "normal",
    axial_gauge_tolerance_t_m=1.0e-9,
    orbit_gauge_tolerance_t_m=1.0e-9,
    symmetry_tolerance_t_m=1.0e-9,
) -> EarlyTimesHCurlFieldCertificate:
    """Certify the constrained HCurl field accepted by the Radia Lie map.

    Both sides of the median plane must be sampled on a symmetric ``y`` grid.
    This function does not repair or invent missing off-plane data.  It rejects
    a field whose gauge or parity is incompatible with the specialized
    EarlyTimes contract before any local Taylor jet is recovered.
    """
    if not isinstance(data, TransverseVectorPotentialData):
        raise TypeError("data must be TransverseVectorPotentialData")
    symmetry = str(symmetry_class).lower()
    if symmetry not in {"normal", "skew", "general"}:
        raise ValueError("symmetry_class must be 'normal', 'skew', or 'general'")
    order = int(required_order)
    if isinstance(required_order, bool) or order != required_order or order < 1:
        raise ValueError("required_order must be a positive integer")
    tolerances = np.asarray(
        [
            axial_gauge_tolerance_t_m,
            orbit_gauge_tolerance_t_m,
            symmetry_tolerance_t_m,
        ],
        dtype=float,
    )
    if np.any(~np.isfinite(tolerances)) or np.any(tolerances <= 0.0):
        raise ValueError("HCurl gauge/symmetry tolerances must be finite and positive")
    if (
        data.grid_function_space_order is not None
        and data.grid_function_space_order < order
    ):
        raise ValueError(
            f"EarlyTimes fourth-order Lie input requires HCurl order >= {order}"
        )

    coordinate_scale = max(
        1.0,
        float(np.max(np.abs(data.y_m), initial=0.0)),
    )
    coordinate_tolerance = 128.0 * np.finfo(float).eps * coordinate_scale
    if not np.allclose(
        data.y_m,
        -data.y_m[::-1],
        rtol=0.0,
        atol=coordinate_tolerance,
    ):
        raise ValueError(
            "median-plane certification requires symmetric upper/lower y samples"
        )

    maximum_Ax = float(np.max(np.abs(data.Ax_t_m), initial=0.0))
    orbit_x, orbit_y = data.design_orbit_indices
    maximum_orbit = float(
        max(
            np.max(np.abs(data.Ay_t_m[:, orbit_x, orbit_y]), initial=0.0),
            np.max(np.abs(data.As_t_m[:, orbit_x, orbit_y]), initial=0.0),
        )
    )
    Ax_even_defect = float(
        np.max(np.abs(data.Ax_t_m - data.Ax_t_m[:, :, ::-1]), initial=0.0)
    )
    Ay_odd_defect = float(
        np.max(np.abs(data.Ay_t_m + data.Ay_t_m[:, :, ::-1]), initial=0.0)
    )
    As_even_defect = float(
        np.max(np.abs(data.As_t_m - data.As_t_m[:, :, ::-1]), initial=0.0)
    )
    Ay_even_defect = float(
        np.max(np.abs(data.Ay_t_m - data.Ay_t_m[:, :, ::-1]), initial=0.0)
    )
    As_odd_defect = float(
        np.max(np.abs(data.As_t_m + data.As_t_m[:, :, ::-1]), initial=0.0)
    )
    axial_tolerance, orbit_tolerance, symmetry_tolerance = tolerances
    if maximum_Ax > axial_tolerance:
        raise ValueError(
            "EarlyTimes HCurl axial gauge failed: "
            f"max |Ax|={maximum_Ax:.6e} T*m"
        )
    if maximum_orbit > orbit_tolerance:
        raise ValueError(
            "EarlyTimes HCurl design-orbit gauge failed: max |Ay|/|As|="
            f"{maximum_orbit:.6e} T*m"
        )
    if symmetry == "normal":
        maximum_symmetry = max(Ax_even_defect, Ay_odd_defect, As_even_defect)
    elif symmetry == "skew":
        maximum_symmetry = max(Ax_even_defect, Ay_even_defect, As_odd_defect)
    else:
        maximum_symmetry = 0.0
    if maximum_symmetry > symmetry_tolerance:
        raise ValueError(
            "EarlyTimes HCurl median-plane symmetry failed: max parity defect="
            f"{maximum_symmetry:.6e} T*m"
        )
    return EarlyTimesHCurlFieldCertificate(
        grid_function_space_class=data.grid_function_space_class,
        grid_function_space_order=data.grid_function_space_order,
        required_order=order,
        symmetry_class=symmetry,
        maximum_Ax_t_m=maximum_Ax,
        maximum_orbit_Ay_As_t_m=maximum_orbit,
        maximum_Ax_even_parity_defect_t_m=Ax_even_defect,
        maximum_Ay_odd_parity_defect_t_m=Ay_odd_defect,
        maximum_As_even_parity_defect_t_m=As_even_defect,
        maximum_Ay_even_parity_defect_t_m=Ay_even_defect,
        maximum_As_odd_parity_defect_t_m=As_odd_defect,
        axial_gauge_tolerance_t_m=float(axial_tolerance),
        orbit_gauge_tolerance_t_m=float(orbit_tolerance),
        symmetry_tolerance_t_m=float(symmetry_tolerance),
    )


def _triangular_powers(degree: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        (x_power, total_degree - x_power)
        for total_degree in range(1, degree + 1)
        for x_power in range(total_degree, -1, -1)
    )


def _scaled_polynomial_fit(
    values: np.ndarray,
    x_values: np.ndarray,
    y_values: np.ndarray,
    powers: tuple[tuple[int, int], ...],
    x_radius: float,
    y_radius: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    design = np.column_stack(
        [
            (x_values / x_radius) ** x_power
            * (y_values / y_radius) ** y_power
            for x_power, y_power in powers
        ]
    )
    if np.linalg.matrix_rank(design) < len(powers):
        raise ValueError("transverse samples do not determine the requested jet")
    condition = float(np.linalg.cond(design))
    scaled_operator = np.linalg.pinv(design)
    scales = np.asarray(
        [x_radius**x_power * y_radius**y_power for x_power, y_power in powers]
    )
    operator = scaled_operator / scales[:, None]
    coefficients = values @ operator.T
    fitted = coefficients @ np.column_stack(
        [x_values**x_power * y_values**y_power for x_power, y_power in powers]
    ).T
    residual = float(np.max(np.abs(fitted - values), initial=0.0))
    return coefficients, operator, scales, condition, residual


def fit_transverse_vector_potential_polynomials(
    data: TransverseVectorPotentialData,
    *,
    degree: int,
    symmetry_class: str = "normal",
    axial_gauge_tolerance_t_m=1.0e-9,
    orbit_gauge_tolerance_t_m=1.0e-9,
    symmetry_tolerance_t_m=1.0e-9,
    fit_tolerance_t_m=None,
    left_right_tolerance_t_m=None,
) -> TransverseVectorPotentialPolynomialFit:
    """Recover the full transverse A jet and its exact sample Jacobian.

    This is an internal representation step after HCurl-field certification;
    it is not an alternative public field-map input.  The fit uses every
    upper/lower sample and all declared-symmetry monomials ``x**i*y**j`` with
    ``1 <= i+j <= degree``.  Normal symmetry uses ``A_y`` odd / ``A_s`` even;
    skew symmetry uses ``A_y`` even / ``A_s`` odd; general mode retains both.
    Independent
    ``x<=0`` and ``x>=0`` fits expose the derivative mismatch across the
    design-orbit face without averaging two HCurl normal traces.
    """
    if not isinstance(data, TransverseVectorPotentialData):
        raise TypeError("data must be TransverseVectorPotentialData")
    symmetry = str(symmetry_class).lower()
    if symmetry not in {"normal", "skew", "general"}:
        raise ValueError("symmetry_class must be 'normal', 'skew', or 'general'")
    polynomial_degree = int(degree)
    if (
        isinstance(degree, bool)
        or polynomial_degree != degree
        or polynomial_degree < 1
        or polynomial_degree > 5
    ):
        raise ValueError("degree must be an integer from one through five")
    axial_tolerance = float(axial_gauge_tolerance_t_m)
    orbit_tolerance = float(orbit_gauge_tolerance_t_m)
    fit_tolerance = None if fit_tolerance_t_m is None else float(fit_tolerance_t_m)
    side_tolerance = (
        None if left_right_tolerance_t_m is None else float(left_right_tolerance_t_m)
    )
    tolerances = [axial_tolerance, orbit_tolerance]
    tolerances.extend(
        value for value in (fit_tolerance, side_tolerance) if value is not None
    )
    if any(not np.isfinite(value) or value <= 0.0 for value in tolerances):
        raise ValueError("gauge/fit/left-right tolerances must be finite and positive")

    field_certificate = certify_earlytimes_hcurl_vector_potential(
        data,
        required_order=polynomial_degree,
        symmetry_class=symmetry,
        axial_gauge_tolerance_t_m=axial_tolerance,
        orbit_gauge_tolerance_t_m=orbit_tolerance,
        symmetry_tolerance_t_m=symmetry_tolerance_t_m,
    )
    maximum_Ax = field_certificate.maximum_Ax_t_m
    maximum_orbit = field_certificate.maximum_orbit_Ay_As_t_m

    x_radius = float(np.max(np.abs(data.x_m), initial=0.0))
    y_radius = float(np.max(np.abs(data.y_m), initial=0.0))
    if min(x_radius, y_radius) <= 0.0:
        raise ValueError("transverse samples must span nonzero x and y")
    x_grid, y_grid = np.meshgrid(data.x_m, data.y_m, indexing="ij")
    x_flat = x_grid.reshape(-1)
    y_flat = y_grid.reshape(-1)
    Ay_values = data.Ay_t_m.reshape(len(data.s_m), -1)
    As_values = data.As_t_m.reshape(len(data.s_m), -1)
    canonical_powers = _triangular_powers(polynomial_degree)
    Ay_powers = tuple(
        power
        for power in canonical_powers
        if symmetry == "general"
        or (symmetry == "normal" and power[1] % 2 == 1)
        or (symmetry == "skew" and power[1] % 2 == 0)
    )
    As_powers = tuple(
        power
        for power in canonical_powers
        if symmetry == "general"
        or (symmetry == "normal" and power[1] % 2 == 0)
        or (symmetry == "skew" and power[1] % 2 == 1)
    )

    Ay_fit, Ay_operator, Ay_scales, Ay_condition, Ay_residual = (
        _scaled_polynomial_fit(
            Ay_values,
            x_flat,
            y_flat,
            Ay_powers,
            x_radius,
            y_radius,
        )
    )
    As_fit, As_operator, As_scales, As_condition, As_residual = (
        _scaled_polynomial_fit(
            As_values,
            x_flat,
            y_flat,
            As_powers,
            x_radius,
            y_radius,
        )
    )
    if fit_tolerance is not None and max(Ay_residual, As_residual) > fit_tolerance:
        raise RuntimeError(
            "transverse A polynomial fit exceeds fit_tolerance_t_m: "
            f"maximum_Ay_fit_residual_t_m={Ay_residual:.17g}, "
            f"maximum_As_fit_residual_t_m={As_residual:.17g}, "
            f"fit_tolerance_t_m={fit_tolerance:.17g}"
        )

    coordinate_tolerance = 64.0 * np.finfo(float).eps * max(x_radius, 1.0)
    left = x_flat <= coordinate_tolerance
    right = x_flat >= -coordinate_tolerance
    left_Ay, _, _, left_Ay_condition, _ = _scaled_polynomial_fit(
        Ay_values[:, left],
        x_flat[left],
        y_flat[left],
        Ay_powers,
        x_radius,
        y_radius,
    )
    right_Ay, _, _, right_Ay_condition, _ = _scaled_polynomial_fit(
        Ay_values[:, right],
        x_flat[right],
        y_flat[right],
        Ay_powers,
        x_radius,
        y_radius,
    )
    left_As, _, _, left_As_condition, _ = _scaled_polynomial_fit(
        As_values[:, left],
        x_flat[left],
        y_flat[left],
        As_powers,
        x_radius,
        y_radius,
    )
    right_As, _, _, right_As_condition, _ = _scaled_polynomial_fit(
        As_values[:, right],
        x_flat[right],
        y_flat[right],
        As_powers,
        x_radius,
        y_radius,
    )
    side_discrepancy = float(
        max(
            np.max(np.abs((left_Ay - right_Ay) * Ay_scales), initial=0.0),
            np.max(np.abs((left_As - right_As) * As_scales), initial=0.0),
        )
    )
    if side_tolerance is not None and side_discrepancy > side_tolerance:
        raise RuntimeError(
            "left/right transverse A jets exceed left_right_tolerance_t_m: "
            "maximum_left_right_scaled_coefficient_discrepancy_t_m="
            f"{side_discrepancy:.17g}, "
            f"left_right_tolerance_t_m={side_tolerance:.17g}"
        )

    station_count = len(data.s_m)
    coefficient_shape = (
        station_count,
        polynomial_degree + 1,
        polynomial_degree + 1,
    )
    Ay_coefficients = np.zeros(coefficient_shape)
    As_coefficients = np.zeros(coefficient_shape)
    Ay_sample_operator = np.zeros((len(canonical_powers), x_flat.size))
    As_sample_operator = np.zeros_like(Ay_sample_operator)
    canonical_index = {power: index for index, power in enumerate(canonical_powers)}
    for local, (x_power, y_power) in enumerate(Ay_powers):
        Ay_coefficients[:, x_power, y_power] = Ay_fit[:, local]
        Ay_sample_operator[canonical_index[(x_power, y_power)]] = Ay_operator[local]
    for local, (x_power, y_power) in enumerate(As_powers):
        As_coefficients[:, x_power, y_power] = As_fit[:, local]
        As_sample_operator[canonical_index[(x_power, y_power)]] = As_operator[local]

    return TransverseVectorPotentialPolynomialFit(
        data=data,
        field_certificate=field_certificate,
        degree=polynomial_degree,
        symmetry_class=symmetry,
        canonical_powers=canonical_powers,
        Ay_coefficients_t_m=np.ascontiguousarray(Ay_coefficients),
        As_coefficients_t_m=np.ascontiguousarray(As_coefficients),
        Ay_sample_to_coefficients=np.ascontiguousarray(Ay_sample_operator),
        As_sample_to_coefficients=np.ascontiguousarray(As_sample_operator),
        maximum_Ax_t_m=maximum_Ax,
        maximum_orbit_Ay_As_t_m=maximum_orbit,
        maximum_Ay_fit_residual_t_m=Ay_residual,
        maximum_As_fit_residual_t_m=As_residual,
        maximum_left_right_scaled_coefficient_discrepancy_t_m=side_discrepancy,
        scaled_Ay_design_condition=Ay_condition,
        scaled_As_design_condition=As_condition,
        scaled_left_Ay_design_condition=left_Ay_condition,
        scaled_right_Ay_design_condition=right_Ay_condition,
        scaled_left_As_design_condition=left_As_condition,
        scaled_right_As_design_condition=right_As_condition,
    )


def sample_transverse_vector_potential(
    vector_potential,
    orbit,
    x_m,
    y_m,
    *,
    s_m=None,
) -> TransverseVectorPotentialData:
    """Sample HCurl ``A`` throughout a local transverse ``(x,y)`` patch.

    NGSolve retains ownership of the HCurl basis, Piola mapping, element
    orientation, curved geometry, and point evaluation.  Radia only maps the
    physical vectors into the design-orbit frame.  Both transverse arrays must
    contain zero so the gauge residual remains an explicit sampled quantity.

    The caller owns the surrounding ``ngsolve.TaskManager``.
    """
    from .accelerator_magnet_topopt import PlanarDesignOrbit

    if not isinstance(orbit, PlanarDesignOrbit):
        raise TypeError("orbit must be a PlanarDesignOrbit")
    space = getattr(vector_potential, "space", None)
    if space is None or "HCurl" not in type(space).__name__:
        raise TypeError(
            "vector_potential must be an NGSolve GridFunction on HCurl(order=p)"
        )
    x_offsets = _finite_array(x_m, "x_m").reshape(-1)
    y_offsets = _finite_array(y_m, "y_m").reshape(-1)
    for offsets, name in ((x_offsets, "x_m"), (y_offsets, "y_m")):
        if offsets.size < 2 or np.any(np.diff(offsets) <= 0.0):
            raise ValueError(
                f"{name} must be a strictly increasing array with at least two points"
            )
        scale = max(1.0, float(np.max(np.abs(offsets))))
        if np.min(np.abs(offsets)) > 64.0 * np.finfo(float).eps * scale:
            raise ValueError(f"{name} must contain the design-orbit coordinate zero")
    if s_m is None:
        stations = orbit.arc_length_stations
        longitudinal = 0.5 * (stations[:-1] + stations[1:])
    else:
        longitudinal = _finite_array(s_m, "s_m").reshape(-1)
        if longitudinal.size == 0:
            raise ValueError("s_m must not be empty")

    horizontal, vertical, tangent = orbit.frame_at(longitudinal)
    reference = orbit.position_at(longitudinal)
    points = (
        reference[:, None, None, :]
        + x_offsets[None, :, None, None] * horizontal[:, None, None, :]
        + y_offsets[None, None, :, None] * vertical[:, None, None, :]
    )
    mesh = space.mesh
    sampled = np.empty_like(points)
    for station in range(len(longitudinal)):
        for x_index in range(len(x_offsets)):
            for y_index in range(len(y_offsets)):
                try:
                    value = vector_potential(mesh(*points[station, x_index, y_index]))
                except Exception as error:
                    raise ValueError(
                        "cannot evaluate HCurl vector potential at transverse "
                        f"sample (s={station}, x={x_index}, y={y_index})"
                    ) from error
                sampled[station, x_index, y_index] = _finite_array(
                    value,
                    "transverse vector potential",
                    shape=(3,),
                )
    Ax = np.einsum("sxyi,si->sxy", sampled, horizontal)
    Ay = np.einsum("sxyi,si->sxy", sampled, vertical)
    As = np.einsum("sxyi,si->sxy", sampled, tangent)
    order_value = getattr(space, "globalorder", None)
    if order_value is None:
        order_value = getattr(space, "order", None)
    try:
        order_value = None if order_value is None else int(order_value)
    except (TypeError, ValueError):
        order_value = None
    return TransverseVectorPotentialData(
        s_m=np.ascontiguousarray(longitudinal),
        x_m=np.ascontiguousarray(x_offsets),
        y_m=np.ascontiguousarray(y_offsets),
        global_points_m=np.ascontiguousarray(points),
        horizontal=np.ascontiguousarray(horizontal),
        vertical=np.ascontiguousarray(vertical),
        tangent=np.ascontiguousarray(tangent),
        Ax_t_m=np.ascontiguousarray(Ax),
        Ay_t_m=np.ascontiguousarray(Ay),
        As_t_m=np.ascontiguousarray(As),
        grid_function_space_class=type(space).__name__,
        grid_function_space_order=order_value,
    )


def build_curvilinear_beam_mesh(
    orbit,
    *,
    half_width_m: float,
    half_height_m: float,
    maxh_m: float,
    inner_half_width_m: float | None = None,
    initial_horizontal=None,
    curve_order: int = 1,
    material: str = "beam_tube",
    boundary: str = "beam_tube_boundary",
) -> CurvilinearBeamMesh:
    """Build the EarlyTimes four-QUAD/four-HEX Bishop/RMF loft chain.

    ``orbit`` is a :class:`PlanarDesignOrbit`.  Its plane normal seeds the
    frame so that local ``y`` equals ``bend_axis`` and local ``x`` lies in the
    median plane.  Passing ``initial_horizontal`` overrides that convention.
    Each section has x nodes ``[-a,-c,0,c,a]`` and only the two y nodes
    ``[-b,b]``.  Four quadrilateral strips are matched between consecutive
    sections, producing exactly four HEX elements per longitudinal cell.
    Consequently the median plane lies inside the HEX volume, while the
    design orbit lies inside the central ``x=0`` face.  ``A_s`` and ``A_y``
    are tangential HCurl traces on that face; ``A_x`` is its normal component.

    ``maxh_m`` controls the maximum station spacing along the design orbit.
    The geometry between stations is the NGSolve HEX mapping, so longitudinal
    refinement remains the convergence control even when ``curve_order`` is
    greater than one.  No regular Cartesian field grid is made.

    This first mesh constructor is for an open beam line.  A closed ring needs
    a periodic mesh seam and is rejected explicitly for now.
    """
    from .accelerator_magnet_topopt import PlanarDesignOrbit

    if not isinstance(orbit, PlanarDesignOrbit):
        raise TypeError("orbit must be a PlanarDesignOrbit")
    width = float(half_width_m)
    height = float(half_height_m)
    maxh = float(maxh_m)
    inner_width = width / 2.0 if inner_half_width_m is None else float(
        inner_half_width_m
    )
    if not np.all(np.isfinite([width, height, maxh])) or min(width, height, maxh) <= 0.0:
        raise ValueError(
            "half_width_m, half_height_m, and maxh_m must be finite and positive"
        )
    if not np.isfinite(inner_width) or not 0.0 < inner_width < width:
        raise ValueError(
            "inner_half_width_m must be finite and satisfy "
            "0 < inner_half_width_m < half_width_m"
        )
    order = int(curve_order)
    if isinstance(curve_order, bool) or order != curve_order or order < 1:
        raise ValueError("curve_order must be a positive integer")
    if not material or not boundary:
        raise ValueError("material and boundary must be non-empty strings")
    scale = max(1.0, float(np.max(np.linalg.norm(
        orbit.positions - orbit.positions[0], axis=1
    ))))
    if np.linalg.norm(orbit.positions[-1] - orbit.positions[0]) <= 1.0e-12 * scale:
        raise NotImplementedError(
            "closed-ring beam meshes need a periodic mesh seam; provide an open "
            "orbit segment for this constructor"
        )

    original_stations = orbit.arc_length_stations
    longitudinal = [float(original_stations[0])]
    for left, right in pairwise(original_stations):
        subdivisions = max(1, int(np.ceil((right - left) / maxh)))
        longitudinal.extend(
            np.linspace(left, right, subdivisions + 1, dtype=float)[1:].tolist()
        )
    longitudinal = np.asarray(longitudinal, dtype=float)
    positions = np.asarray(orbit.position_at(longitudinal), dtype=float)
    tangents = np.asarray(orbit.tangent_at(longitudinal), dtype=float)

    seed = (
        np.cross(orbit.bend_axis, tangents[0])
        if initial_horizontal is None
        else initial_horizontal
    )
    frame = bishop_rmf_frame(
        positions,
        tangents,
        initial_horizontal=seed,
        periodic=False,
    )

    from netgen.meshing import (
        Element2D,
        Element3D,
        FaceDescriptor,
        MeshPoint,
        Pnt,
    )
    from netgen.meshing import Mesh as NetgenMesh
    from ngsolve import Mesh as NGSolveMesh

    x_nodes = np.asarray(
        [-width, -inner_width, 0.0, inner_width, width], dtype=float
    )
    y_nodes = np.asarray([-height, height], dtype=float)
    vertices = (
        frame.positions_m[:, None, None, :]
        + x_nodes[None, None, :, None] * frame.horizontal[:, None, None, :]
        + y_nodes[None, :, None, None] * frame.vertical[:, None, None, :]
    )

    ngmesh = NetgenMesh(dim=3)
    ngmesh.SetMaterial(1, str(material))
    point_ids = np.empty(vertices.shape[:-1], dtype=object)
    for station in range(len(longitudinal)):
        for y_index in range(2):
            for x_index in range(5):
                point_ids[station, y_index, x_index] = ngmesh.Add(
                    MeshPoint(Pnt(*vertices[station, y_index, x_index]))
                )

    connectivity = []
    for station in range(len(longitudinal) - 1):
        for strip in range(4):
            gmsh_order = [
                point_ids[station, 0, strip],
                point_ids[station, 0, strip + 1],
                point_ids[station, 1, strip + 1],
                point_ids[station, 1, strip],
                point_ids[station + 1, 0, strip],
                point_ids[station + 1, 0, strip + 1],
                point_ids[station + 1, 1, strip + 1],
                point_ids[station + 1, 1, strip],
            ]
            netgen_order = [gmsh_order[index] for index in (0, 1, 5, 4, 3, 2, 6, 7)]
            ngmesh.Add(Element3D(1, netgen_order))
            connectivity.append([
                (station * 2 + 0) * 5 + strip,
                (station * 2 + 0) * 5 + strip + 1,
                ((station + 1) * 2 + 0) * 5 + strip + 1,
                ((station + 1) * 2 + 0) * 5 + strip,
                (station * 2 + 1) * 5 + strip,
                (station * 2 + 1) * 5 + strip + 1,
                ((station + 1) * 2 + 1) * 5 + strip + 1,
                ((station + 1) * 2 + 1) * 5 + strip,
            ])

    face = ngmesh.Add(FaceDescriptor(bc=1, domin=1, domout=0))
    ngmesh.SetBCName(0, str(boundary))

    def add_surface(nodes):
        ngmesh.Add(Element2D(face, list(nodes)))

    last = len(longitudinal) - 1
    for strip in range(4):
        add_surface([
            point_ids[0, 0, strip], point_ids[0, 0, strip + 1],
            point_ids[0, 1, strip + 1], point_ids[0, 1, strip],
        ])
        add_surface([
            point_ids[last, 0, strip], point_ids[last, 1, strip],
            point_ids[last, 1, strip + 1], point_ids[last, 0, strip + 1],
        ])
    for station in range(last):
        for strip in range(4):
            add_surface([
                point_ids[station, 0, strip],
                point_ids[station + 1, 0, strip],
                point_ids[station + 1, 0, strip + 1],
                point_ids[station, 0, strip + 1],
            ])
            add_surface([
                point_ids[station, 1, strip],
                point_ids[station, 1, strip + 1],
                point_ids[station + 1, 1, strip + 1],
                point_ids[station + 1, 1, strip],
            ])
        add_surface([
            point_ids[station, 0, 0], point_ids[station, 1, 0],
            point_ids[station + 1, 1, 0], point_ids[station + 1, 0, 0],
        ])
        add_surface([
            point_ids[station, 0, 4], point_ids[station + 1, 0, 4],
            point_ids[station + 1, 1, 4], point_ids[station, 1, 4],
        ])

    mesh = NGSolveMesh(ngmesh)
    if order > 1:
        mesh.Curve(order)
    volumes = np.asarray(
        __import__("ngsolve").Integrate(1.0, mesh, element_wise=True),
        dtype=float,
    )
    if volumes.shape != (4 * last,) or np.any(volumes <= 0.0):
        raise RuntimeError("four-HEX loft chain contains an invalid element")
    return CurvilinearBeamMesh(
        mesh=mesh,
        shape=None,
        frame=frame,
        half_width_m=width,
        half_height_m=height,
        maxh_m=maxh,
        material=str(material),
        boundary=str(boundary),
        curve_order=order,
        topology="bishop-rmf-four-quad-loft-chain",
        x_nodes_m=x_nodes,
        y_nodes_m=y_nodes,
        longitudinal_stations_m=longitudinal,
        vertices_m=np.ascontiguousarray(vertices),
        hex_connectivity=np.ascontiguousarray(connectivity, dtype=np.int64),
    )


def project_design_orbit_gauge(
    vector_potential_coefficient,
    beam_mesh: CurvilinearBeamMesh,
    *,
    order: int = 5,
    constraint_subdivisions: int = 2,
    rbf_width_factor: float = 4.0,
    maximum_rbf_condition: float = 1.0e14,
    gauge_tolerance: float = 1.0e-6,
    verify_curl: bool = True,
    name: str = "A_orbit_gauge",
) -> OrbitGaugeVectorPotential:
    """Project ``A`` into HCurl with ``A_s=A_y=0`` on the design orbit.

    A smooth scalar gauge potential is built by Hermite radial-basis
    interpolation of the station-wise constraints
    ``grad(chi).tangent=-A_s`` and ``grad(chi).vertical=-A_y``.  It is then
    interpolated into ``H1(order+1)`` and its NGSolve-native gradient is added
    before canonical interpolation into ``HCurl(order)``.  Consequently the
    correction belongs to the discrete de Rham gradient subspace and changes
    ``curl(A)`` only by interpolation roundoff.

    Each basis uses distance only in the plane normal to its local horizontal
    axis.  It is therefore constant along local ``x`` lines and does not undo
    an already-applied ``A_x=0`` gauge on a straight loft cell; curvature
    leakage remains controlled by the final HCurl certificate.

    ``constraint_subdivisions=2`` includes every original orbit station and
    every segment midpoint used by EarlyTimes' field fit.  Larger values add
    longitudinal gauge constraints between those points.

    Odd and even positive orders are both supported.  For a decapole-level
    fourth-order LIE map, ``order=5`` is the minimum polynomial candidate;
    ``order=6`` generally gives a tighter orbit-gauge residual on curved
    meshes.  The caller owns :class:`ngsolve.TaskManager`.

    The broad default RBF needed by the projected p=5 gauge can have a
    condition number around ``1e13`` on short curved chains.  The condition
    number is returned for audit, while the accepted field must still satisfy
    both the orbit-gauge residual and the curl-invariance certificate below.
    """
    if not isinstance(beam_mesh, CurvilinearBeamMesh):
        raise TypeError("beam_mesh must be a CurvilinearBeamMesh")
    p = int(order)
    if isinstance(order, bool) or p != order or p < 1:
        raise ValueError("order must be a positive integer")
    subdivisions = int(constraint_subdivisions)
    if (
        isinstance(constraint_subdivisions, bool)
        or subdivisions != constraint_subdivisions
        or subdivisions < 1
    ):
        raise ValueError("constraint_subdivisions must be a positive integer")
    width_factor = float(rbf_width_factor)
    maximum_condition = float(maximum_rbf_condition)
    tolerance = float(gauge_tolerance)
    if (
        not np.all(np.isfinite([width_factor, maximum_condition, tolerance]))
        or width_factor <= 0.0
        or maximum_condition <= 1.0
        or tolerance <= 0.0
    ):
        raise ValueError(
            "RBF width/condition and gauge tolerance must be finite and positive"
        )

    import ngsolve as ng

    mesh = beam_mesh.mesh
    base_frame = beam_mesh.frame
    if subdivisions == 1:
        frame = base_frame
    else:
        positions_refined = []
        tangents_refined = []
        for segment in range(len(base_frame.positions_m) - 1):
            for local in range(subdivisions):
                fraction = local / subdivisions
                positions_refined.append(
                    (1.0 - fraction) * base_frame.positions_m[segment]
                    + fraction * base_frame.positions_m[segment + 1]
                )
                tangent = (
                    (1.0 - fraction) * base_frame.tangent[segment]
                    + fraction * base_frame.tangent[segment + 1]
                )
                tangents_refined.append(_normalize(tangent, "refined tangent"))
        positions_refined.append(base_frame.positions_m[-1])
        tangents_refined.append(base_frame.tangent[-1])
        frame = bishop_rmf_frame(
            np.asarray(positions_refined),
            np.asarray(tangents_refined),
            initial_horizontal=base_frame.horizontal[0],
        )
    positions = frame.positions_m
    station_count = len(positions)
    if station_count < 2:
        raise ValueError("the beam frame needs at least two stations")

    sampled = np.empty((station_count, 3))
    for index, point in enumerate(positions):
        try:
            value = vector_potential_coefficient(mesh(*point))
        except Exception as error:
            raise ValueError(
                f"cannot evaluate vector potential at orbit station {index}"
            ) from error
        sampled[index] = _finite_array(
            value, f"vector potential at station {index}", shape=(3,)
        )
    As = np.einsum("ij,ij->i", sampled, frame.tangent)
    Ay = np.einsum("ij,ij->i", sampled, frame.vertical)

    spacing = np.linalg.norm(np.diff(positions, axis=0), axis=1)
    width = width_factor * float(np.mean(spacing))
    if not np.isfinite(width) or width <= 0.0:
        raise ValueError("orbit station spacing must be finite and positive")
    directions = np.concatenate((frame.tangent, frame.vertical), axis=0)
    centers = np.concatenate((positions, positions), axis=0)
    horizontal_directions = np.concatenate(
        (frame.horizontal, frame.horizontal), axis=0
    )
    targets = -np.concatenate((As, Ay))
    count = 2 * station_count
    collocation = np.empty((count, count))
    width_squared = width * width
    for row, (evaluation_point, evaluation_direction) in enumerate(
        zip(centers, directions)
    ):
        for column, (center, basis_direction, basis_horizontal) in enumerate(
            zip(centers, directions, horizontal_directions)
        ):
            displacement = evaluation_point - center
            horizontal_displacement = float(
                np.dot(basis_horizontal, displacement)
            )
            transverse_displacement = (
                displacement
                - horizontal_displacement * basis_horizontal
            )
            radial = np.exp(
                -float(
                    np.dot(
                        transverse_displacement,
                        transverse_displacement,
                    )
                )
                / (2.0 * width_squared)
            )
            gradient = (
                basis_direction
                - float(np.dot(basis_direction, displacement))
                * transverse_displacement
                / width_squared
            ) * radial
            collocation[row, column] = float(
                np.dot(evaluation_direction, gradient)
            )
    condition = float(np.linalg.cond(collocation))
    if not np.isfinite(condition) or condition > maximum_condition:
        raise RuntimeError(
            "design-orbit gauge RBF is ill-conditioned "
            f"(condition={condition:.6g}, limit={maximum_condition:.6g}); "
            "adjust station spacing or rbf_width_factor"
        )
    coefficients = np.linalg.solve(collocation, targets)

    coordinates = (ng.x, ng.y, ng.z)
    chi_cf = ng.CoefficientFunction(0.0)
    for coefficient, center, direction, horizontal in zip(
        coefficients, centers, directions, horizontal_directions
    ):
        displacement = tuple(
            coordinates[axis] - float(center[axis]) for axis in range(3)
        )
        linear = sum(
            float(direction[axis]) * displacement[axis] for axis in range(3)
        )
        axial_displacement = sum(
            float(horizontal[axis]) * displacement[axis]
            for axis in range(3)
        )
        radius_squared = (
            sum(value * value for value in displacement)
            - axial_displacement * axial_displacement
        )
        chi_cf = chi_cf + float(coefficient) * linear * ng.exp(
            -radius_squared / (2.0 * width_squared)
        )

    gauge_potential = ng.GridFunction(
        ng.H1(mesh, order=p + 1), name=f"{name}_chi"
    )
    space = ng.HCurl(mesh, order=p)
    ungauged = ng.GridFunction(space, name=f"{name}_ungauged")
    gauged = ng.GridFunction(space, name=name)
    # NGSolve's SIMD local L2 projection of p=5 HCurl HEX fields can exceed
    # the TaskManager's fixed 10 MB per-task LocalHeap.  The public Set option
    # selects the equivalent scalar NGSolve projector without duplicating FE
    # basis/orientation logic in Radia.
    gauge_potential.Set(chi_cf, use_simd=False)
    ungauged.Set(vector_potential_coefficient, use_simd=False)
    gauged.Set(
        vector_potential_coefficient + ng.grad(gauge_potential),
        use_simd=False,
    )

    corrected = np.empty((station_count, 3))
    for index, point in enumerate(positions):
        corrected[index] = _finite_array(
            gauged(mesh(*point)),
            f"gauged vector potential at station {index}",
            shape=(3,),
        )
    As_after = np.einsum("ij,ij->i", corrected, frame.tangent)
    Ay_after = np.einsum("ij,ij->i", corrected, frame.vertical)
    maximum_residual = max(
        np.max(np.abs(As_after), initial=0.0),
        np.max(np.abs(Ay_after), initial=0.0),
    )
    if maximum_residual > tolerance:
        raise RuntimeError(
            "design-orbit gauge did not reach gauge_tolerance; increase HCurl "
            "order or refine the curvilinear mesh"
        )

    curl_change = float("nan")
    if verify_curl:
        difference = ng.curl(gauged) - ng.curl(ungauged)
        curl_change = float(
            ng.sqrt(
                ng.Integrate(
                    ng.InnerProduct(difference, difference),
                    mesh,
                    order=2 * p + 2,
                )
            )
        )
    return OrbitGaugeVectorPotential(
        vector_potential=gauged,
        ungauged_vector_potential=ungauged,
        gauge_potential=gauge_potential,
        order=p,
        frame=frame,
        As_before_t_m=As,
        Ay_before_t_m=Ay,
        As_after_t_m=As_after,
        Ay_after_t_m=Ay_after,
        rbf_width_m=width,
        rbf_condition=condition,
        curl_change_l2_t_m32=curl_change,
    )


__all__ = [
    "BishopRMFFrame",
    "CurvilinearBeamMesh",
    "EarlyTimesHCurlFieldCertificate",
    "OrbitGaugeVectorPotential",
    "TransverseVectorPotentialData",
    "TransverseVectorPotentialPolynomialFit",
    "bishop_rmf_frame",
    "build_curvilinear_beam_mesh",
    "certify_earlytimes_hcurl_vector_potential",
    "fit_transverse_vector_potential_polynomials",
    "project_design_orbit_gauge",
    "sample_transverse_vector_potential",
]
