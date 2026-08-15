"""Electromagnetic force post-processing from solved SI field samples.

This module provides the small, solver-independent force kernels shared by
Radia workflows.  NGSolve, Radia, BEM, and reduced-order solvers remain
responsible for producing the fields and quadrature samples; these functions
turn those samples into force with one explicit convention:

* vectors use Cartesian components on the last axis;
* magnetic flux density is in tesla;
* current density is in ampere per square metre;
* volume and area weights are in cubic and square metres;
* returned force is in newtons and torque is in newton metres;
* phasor routines require an explicit peak or RMS amplitude convention.

The Maxwell-stress surface must lie wholly in air or vacuum and its normals
must point out of the body whose force is requested.  For magnetic material,
prefer a weighted-stress or virtual-work calculation and use a contour/surface
integral as an independent sensitivity check.
"""

from __future__ import annotations

import math

import numpy as np

MU0 = 4.0e-7 * math.pi

__all__ = [
    "MU0",
    "air_gap_holding_force",
    "air_gap_maxwell_pressure",
    "air_gap_shear_stress",
    "air_gap_shear_torque",
    "air_gap_shear_torque_from_angle_samples",
    "coenergy_torque_from_angle_samples",
    "force_torque_result",
    "integrate_lorentz_force",
    "integrate_lorentz_force_and_torque",
    "integrate_lorentz_torque",
    "integrate_maxwell_surface_force",
    "integrate_maxwell_surface_force_and_torque",
    "integrate_maxwell_surface_torque",
    "integrate_time_average_lorentz_force",
    "integrate_time_average_lorentz_force_and_torque",
    "integrate_time_average_lorentz_torque",
    "integrate_time_average_maxwell_surface_force",
    "integrate_time_average_maxwell_surface_force_and_torque",
    "integrate_time_average_maxwell_surface_torque",
    "lorentz_force_density",
    "maxwell_stress_tensor_air",
    "maxwell_traction_air",
    "time_average_air_gap_shear_stress",
    "time_average_air_gap_shear_torque_from_angle_samples",
    "time_average_lorentz_force_density",
    "time_average_maxwell_stress_tensor_air",
    "time_average_maxwell_traction_air",
    "virtual_work_force_from_displacement_samples",
]


def _real_vectors(values, name: str) -> np.ndarray:
    array = np.asarray(values)
    if np.iscomplexobj(array):
        raise ValueError(f"{name} must be real for a static force calculation")
    try:
        array = np.asarray(array, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain real numbers") from exc
    if array.ndim == 0 or array.shape[-1] != 3:
        raise ValueError(f"{name} must have shape (..., 3)")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _complex_vectors(values, name: str) -> np.ndarray:
    try:
        array = np.asarray(values, dtype=complex)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain numeric phasors") from exc
    if array.ndim == 0 or array.shape[-1] != 3:
        raise ValueError(f"{name} must have shape (..., 3)")
    if not np.all(np.isfinite(array.real)) or not np.all(np.isfinite(array.imag)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _positive_scalar(value, name: str) -> float:
    scalar = float(value)
    if not math.isfinite(scalar) or scalar <= 0.0:
        raise ValueError(f"{name} must be finite and > 0")
    return scalar


def _nonnegative_scalar(value, name: str) -> float:
    scalar = float(value)
    if not math.isfinite(scalar) or scalar < 0.0:
        raise ValueError(f"{name} must be finite and >= 0")
    return scalar


def _broadcast_vectors(first, first_name: str, second, second_name: str):
    a = _real_vectors(first, first_name)
    b = _real_vectors(second, second_name)
    try:
        return np.broadcast_arrays(a, b)
    except ValueError as exc:
        raise ValueError(
            f"{first_name} and {second_name} leading dimensions are not broadcastable"
        ) from exc


def _broadcast_complex_vectors(first, first_name: str, second, second_name: str):
    a = _complex_vectors(first, first_name)
    b = _complex_vectors(second, second_name)
    try:
        return np.broadcast_arrays(a, b)
    except ValueError as exc:
        raise ValueError(
            f"{first_name} and {second_name} leading dimensions are not broadcastable"
        ) from exc


def _quadrature_sum(values: np.ndarray, weights, name: str) -> np.ndarray:
    weight_array = np.asarray(weights)
    if np.iscomplexobj(weight_array):
        raise ValueError(f"{name} must be real")
    try:
        weight_array = np.asarray(weight_array, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain real numbers") from exc
    if not np.all(np.isfinite(weight_array)):
        raise ValueError(f"{name} must contain only finite values")
    if np.any(weight_array < 0.0):
        raise ValueError(f"{name} must be >= 0")
    try:
        weights_broadcast = np.broadcast_to(weight_array, values.shape[:-1])
    except ValueError as exc:
        raise ValueError(
            f"{name} must be scalar or broadcast to the sample dimensions "
            f"{values.shape[:-1]}"
        ) from exc
    weighted = values * weights_broadcast[..., np.newaxis]
    sample_axes = tuple(range(weighted.ndim - 1))
    if not sample_axes:
        return weighted.copy()
    return np.sum(weighted, axis=sample_axes)


def _phasor_average_factor(amplitude: str) -> float:
    key = str(amplitude).strip().lower()
    if key == "peak":
        return 0.5
    if key == "rms":
        return 1.0
    raise ValueError("amplitude must be 'peak' or 'rms'")


def _force_and_torque_from_samples(
    force_samples: np.ndarray,
    weights,
    sample_points_m,
    pivot_m,
    weight_name: str,
) -> tuple[np.ndarray, np.ndarray]:
    points = _real_vectors(sample_points_m, "sample_points_m")
    try:
        points, forces = np.broadcast_arrays(points, force_samples)
    except ValueError as exc:
        raise ValueError(
            "sample_points_m and force samples leading dimensions are not broadcastable"
        ) from exc
    if pivot_m is None:
        pivot = np.zeros(3, dtype=float)
    else:
        pivot = _real_vectors(pivot_m, "pivot_m")
        if pivot.ndim != 1:
            raise ValueError("pivot_m must be one Cartesian three-vector")
    force = _quadrature_sum(forces, weights, weight_name)
    torque_samples = np.cross(points - pivot, forces)
    torque = _quadrature_sum(torque_samples, weights, weight_name)
    return force, torque


def _real_table(values, name: str, minimum_size: int = 1) -> np.ndarray:
    try:
        array = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain real numbers") from exc
    if array.ndim != 1 or array.size < minimum_size:
        raise ValueError(f"{name} must be a one-dimensional table with at least {minimum_size} samples")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _strictly_increasing(values: np.ndarray, name: str) -> None:
    if np.any(np.diff(values) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")


def _energy_sign(energy_kind: str) -> float:
    key = str(energy_kind).strip().lower().replace("-", "_").replace(" ", "_")
    if key in {"coenergy", "magnetic_coenergy", "constant_current", "w_prime"}:
        return 1.0
    if key in {"stored_energy", "field_energy", "magnetic_energy", "constant_flux"}:
        return -1.0
    raise ValueError(
        "energy_kind must be 'coenergy'/'constant_current' or "
        "'stored_energy'/'constant_flux'"
    )


def lorentz_force_density(current_density_A_per_m2, magnetic_flux_density_T):
    """Return the magnetostatic Lorentz force density ``J x B`` in N/m^3.

    Inputs may be one three-vector or broadcast-compatible arrays whose last
    dimension is three.  ``B`` must be the field acting on the selected
    current distribution; omit a singular self-field when the discretization
    has not regularized it.
    """

    current_density, magnetic_flux_density = _broadcast_vectors(
        current_density_A_per_m2,
        "current_density_A_per_m2",
        magnetic_flux_density_T,
        "magnetic_flux_density_T",
    )
    return np.cross(current_density, magnetic_flux_density)


def integrate_lorentz_force(
    current_density_A_per_m2,
    magnetic_flux_density_T,
    volume_weights_m3,
):
    """Integrate ``J x B`` over volume and return a Cartesian force in N.

    ``volume_weights_m3`` contains the physical quadrature weights, including
    any Jacobian and symmetry-sector factor.  It may be scalar or broadcast to
    the leading sample dimensions of ``J`` and ``B``.
    """

    density = lorentz_force_density(
        current_density_A_per_m2,
        magnetic_flux_density_T,
    )
    return _quadrature_sum(density, volume_weights_m3, "volume_weights_m3")


def integrate_lorentz_force_and_torque(
    current_density_A_per_m2,
    magnetic_flux_density_T,
    volume_weights_m3,
    sample_points_m,
    pivot_m=None,
):
    """Integrate static Lorentz force and torque about ``pivot_m``.

    ``sample_points_m`` are the physical quadrature locations corresponding to
    the J/B rows.  Returns ``(force_N, torque_Nm)`` as Cartesian arrays.
    """

    density = lorentz_force_density(
        current_density_A_per_m2,
        magnetic_flux_density_T,
    )
    return _force_and_torque_from_samples(
        density,
        volume_weights_m3,
        sample_points_m,
        pivot_m,
        "volume_weights_m3",
    )


def integrate_lorentz_torque(
    current_density_A_per_m2,
    magnetic_flux_density_T,
    volume_weights_m3,
    sample_points_m,
    pivot_m=None,
):
    """Integrate static Lorentz torque in N m about ``pivot_m``."""

    return integrate_lorentz_force_and_torque(
        current_density_A_per_m2,
        magnetic_flux_density_T,
        volume_weights_m3,
        sample_points_m,
        pivot_m=pivot_m,
    )[1]


def time_average_lorentz_force_density(
    current_density_phasor_A_per_m2,
    magnetic_flux_density_phasor_T,
    amplitude="peak",
):
    """Return cycle-averaged Lorentz density from complex phasors.

    For peak phasors this is ``0.5*Re(J x conj(B))``; for RMS phasors the
    factor is one.  The returned array is real and has units N/m^3.
    """

    current_density, magnetic_flux_density = _broadcast_complex_vectors(
        current_density_phasor_A_per_m2,
        "current_density_phasor_A_per_m2",
        magnetic_flux_density_phasor_T,
        "magnetic_flux_density_phasor_T",
    )
    factor = _phasor_average_factor(amplitude)
    return factor * np.real(
        np.cross(current_density, np.conjugate(magnetic_flux_density))
    )


def integrate_time_average_lorentz_force(
    current_density_phasor_A_per_m2,
    magnetic_flux_density_phasor_T,
    volume_weights_m3,
    amplitude="peak",
):
    """Integrate cycle-averaged phasor Lorentz force in N."""

    density = time_average_lorentz_force_density(
        current_density_phasor_A_per_m2,
        magnetic_flux_density_phasor_T,
        amplitude=amplitude,
    )
    return _quadrature_sum(density, volume_weights_m3, "volume_weights_m3")


def integrate_time_average_lorentz_force_and_torque(
    current_density_phasor_A_per_m2,
    magnetic_flux_density_phasor_T,
    volume_weights_m3,
    sample_points_m,
    pivot_m=None,
    amplitude="peak",
):
    """Integrate cycle-averaged Lorentz force and torque from phasors."""

    density = time_average_lorentz_force_density(
        current_density_phasor_A_per_m2,
        magnetic_flux_density_phasor_T,
        amplitude=amplitude,
    )
    return _force_and_torque_from_samples(
        density,
        volume_weights_m3,
        sample_points_m,
        pivot_m,
        "volume_weights_m3",
    )


def integrate_time_average_lorentz_torque(
    current_density_phasor_A_per_m2,
    magnetic_flux_density_phasor_T,
    volume_weights_m3,
    sample_points_m,
    pivot_m=None,
    amplitude="peak",
):
    """Integrate cycle-averaged Lorentz torque in N m from phasors."""

    return integrate_time_average_lorentz_force_and_torque(
        current_density_phasor_A_per_m2,
        magnetic_flux_density_phasor_T,
        volume_weights_m3,
        sample_points_m,
        pivot_m=pivot_m,
        amplitude=amplitude,
    )[1]


def maxwell_stress_tensor_air(magnetic_flux_density_T, permeability_H_per_m=MU0):
    """Return the static magnetic Maxwell stress tensor in air, in Pa.

    The tensor is

    ``T = (B tensor B - 0.5 * (B dot B) * I) / permeability``.

    The result has shape ``(..., 3, 3)``.  This real-valued function is for
    static fields; cycle-averaged phasor forces require the conjugated
    time-average identity and an explicit peak/RMS convention.
    """

    magnetic_flux_density = _real_vectors(
        magnetic_flux_density_T,
        "magnetic_flux_density_T",
    )
    permeability = _positive_scalar(permeability_H_per_m, "permeability_H_per_m")
    outer = (
        magnetic_flux_density[..., :, np.newaxis]
        * magnetic_flux_density[..., np.newaxis, :]
    )
    magnitude_squared = np.sum(magnetic_flux_density**2, axis=-1)
    return (
        outer
        - 0.5
        * magnitude_squared[..., np.newaxis, np.newaxis]
        * np.eye(3, dtype=float)
    ) / permeability


def time_average_maxwell_stress_tensor_air(
    magnetic_flux_density_phasor_T,
    permeability_H_per_m=MU0,
    amplitude="peak",
):
    """Return cycle-averaged magnetic Maxwell stress in air, in Pa.

    The peak-phasor convention is
    ``0.5/mu * Re(B tensor conj(B) - 0.5*(B dot conj(B))*I)``.
    RMS phasors omit the leading 0.5.
    """

    magnetic_flux_density = _complex_vectors(
        magnetic_flux_density_phasor_T,
        "magnetic_flux_density_phasor_T",
    )
    permeability = _positive_scalar(permeability_H_per_m, "permeability_H_per_m")
    factor = _phasor_average_factor(amplitude)
    outer = np.real(
        magnetic_flux_density[..., :, np.newaxis]
        * np.conjugate(magnetic_flux_density[..., np.newaxis, :])
    )
    magnitude_squared = np.sum(np.abs(magnetic_flux_density) ** 2, axis=-1)
    return factor * (
        outer
        - 0.5
        * magnitude_squared[..., np.newaxis, np.newaxis]
        * np.eye(3, dtype=float)
    ) / permeability


def maxwell_traction_air(
    magnetic_flux_density_T,
    outward_normal,
    permeability_H_per_m=MU0,
):
    """Return Maxwell traction ``T n`` in air, in Pa.

    Normals are normalized sample by sample.  ``B`` and the normals may be one
    vector or broadcast-compatible arrays with Cartesian components on the
    last axis.
    """

    magnetic_flux_density, normals = _broadcast_vectors(
        magnetic_flux_density_T,
        "magnetic_flux_density_T",
        outward_normal,
        "outward_normal",
    )
    normal_magnitude = np.linalg.norm(normals, axis=-1)
    if np.any(normal_magnitude == 0.0):
        raise ValueError("outward_normal must be nonzero")
    unit_normals = normals / normal_magnitude[..., np.newaxis]
    stress = maxwell_stress_tensor_air(
        magnetic_flux_density,
        permeability_H_per_m=permeability_H_per_m,
    )
    return np.einsum("...ij,...j->...i", stress, unit_normals)


def time_average_maxwell_traction_air(
    magnetic_flux_density_phasor_T,
    outward_normal,
    permeability_H_per_m=MU0,
    amplitude="peak",
):
    """Return cycle-averaged magnetic Maxwell traction ``<T> n`` in Pa."""

    magnetic_flux_density = _complex_vectors(
        magnetic_flux_density_phasor_T,
        "magnetic_flux_density_phasor_T",
    )
    normals = _real_vectors(outward_normal, "outward_normal")
    try:
        magnetic_flux_density, normals = np.broadcast_arrays(
            magnetic_flux_density,
            normals,
        )
    except ValueError as exc:
        raise ValueError(
            "magnetic_flux_density_phasor_T and outward_normal leading dimensions "
            "are not broadcastable"
        ) from exc
    normal_magnitude = np.linalg.norm(normals, axis=-1)
    if np.any(normal_magnitude == 0.0):
        raise ValueError("outward_normal must be nonzero")
    unit_normals = normals / normal_magnitude[..., np.newaxis]
    stress = time_average_maxwell_stress_tensor_air(
        magnetic_flux_density,
        permeability_H_per_m=permeability_H_per_m,
        amplitude=amplitude,
    )
    return np.einsum("...ij,...j->...i", stress, unit_normals)


def integrate_maxwell_surface_force(
    magnetic_flux_density_T,
    outward_normal,
    area_weights_m2,
    permeability_H_per_m=MU0,
):
    """Integrate air-side Maxwell traction over a closed surface.

    ``area_weights_m2`` contains physical surface quadrature weights.  The
    integration surface must enclose the requested body, remain in air, and
    use outward normals.  The return value is ``[Fx, Fy, Fz]`` in newtons.
    """

    traction = maxwell_traction_air(
        magnetic_flux_density_T,
        outward_normal,
        permeability_H_per_m=permeability_H_per_m,
    )
    return _quadrature_sum(traction, area_weights_m2, "area_weights_m2")


def integrate_maxwell_surface_force_and_torque(
    magnetic_flux_density_T,
    outward_normal,
    area_weights_m2,
    sample_points_m,
    pivot_m=None,
    permeability_H_per_m=MU0,
):
    """Integrate static Maxwell surface force and torque about ``pivot_m``."""

    traction = maxwell_traction_air(
        magnetic_flux_density_T,
        outward_normal,
        permeability_H_per_m=permeability_H_per_m,
    )
    return _force_and_torque_from_samples(
        traction,
        area_weights_m2,
        sample_points_m,
        pivot_m,
        "area_weights_m2",
    )


def integrate_maxwell_surface_torque(
    magnetic_flux_density_T,
    outward_normal,
    area_weights_m2,
    sample_points_m,
    pivot_m=None,
    permeability_H_per_m=MU0,
):
    """Integrate static Maxwell surface torque in N m about ``pivot_m``."""

    return integrate_maxwell_surface_force_and_torque(
        magnetic_flux_density_T,
        outward_normal,
        area_weights_m2,
        sample_points_m,
        pivot_m=pivot_m,
        permeability_H_per_m=permeability_H_per_m,
    )[1]


def integrate_time_average_maxwell_surface_force(
    magnetic_flux_density_phasor_T,
    outward_normal,
    area_weights_m2,
    permeability_H_per_m=MU0,
    amplitude="peak",
):
    """Integrate cycle-averaged Maxwell surface force from phasor B."""

    traction = time_average_maxwell_traction_air(
        magnetic_flux_density_phasor_T,
        outward_normal,
        permeability_H_per_m=permeability_H_per_m,
        amplitude=amplitude,
    )
    return _quadrature_sum(traction, area_weights_m2, "area_weights_m2")


def integrate_time_average_maxwell_surface_force_and_torque(
    magnetic_flux_density_phasor_T,
    outward_normal,
    area_weights_m2,
    sample_points_m,
    pivot_m=None,
    permeability_H_per_m=MU0,
    amplitude="peak",
):
    """Integrate cycle-averaged Maxwell surface force and torque."""

    traction = time_average_maxwell_traction_air(
        magnetic_flux_density_phasor_T,
        outward_normal,
        permeability_H_per_m=permeability_H_per_m,
        amplitude=amplitude,
    )
    return _force_and_torque_from_samples(
        traction,
        area_weights_m2,
        sample_points_m,
        pivot_m,
        "area_weights_m2",
    )


def integrate_time_average_maxwell_surface_torque(
    magnetic_flux_density_phasor_T,
    outward_normal,
    area_weights_m2,
    sample_points_m,
    pivot_m=None,
    permeability_H_per_m=MU0,
    amplitude="peak",
):
    """Integrate cycle-averaged Maxwell surface torque in N m."""

    return integrate_time_average_maxwell_surface_force_and_torque(
        magnetic_flux_density_phasor_T,
        outward_normal,
        area_weights_m2,
        sample_points_m,
        pivot_m=pivot_m,
        permeability_H_per_m=permeability_H_per_m,
        amplitude=amplitude,
    )[1]


def air_gap_maxwell_pressure(
    magnetic_flux_density_normal_T,
    permeability_H_per_m=MU0,
):
    """Return uniform normal-field Maxwell pressure ``B_n^2/(2 mu)`` in Pa."""

    permeability = _positive_scalar(permeability_H_per_m, "permeability_H_per_m")
    field = float(magnetic_flux_density_normal_T)
    if not math.isfinite(field):
        raise ValueError("magnetic_flux_density_normal_T must be finite")
    return field * field / (2.0 * permeability)


def air_gap_holding_force(
    magnetic_flux_density_normal_T,
    active_area_m2,
    faces=1,
    permeability_H_per_m=MU0,
):
    """Return a uniform-gap holding-force estimate in N."""

    area = _nonnegative_scalar(active_area_m2, "active_area_m2")
    face_count = int(faces)
    if face_count < 1 or face_count != float(faces):
        raise ValueError("faces must be a positive integer")
    return (
        air_gap_maxwell_pressure(
            magnetic_flux_density_normal_T,
            permeability_H_per_m=permeability_H_per_m,
        )
        * area
        * face_count
    )


def air_gap_shear_stress(
    magnetic_flux_density_radial_T,
    magnetic_flux_density_tangential_T,
    permeability_H_per_m=MU0,
):
    """Return static cylindrical air-gap shear stress ``B_r B_t/mu`` in Pa."""

    permeability = _positive_scalar(permeability_H_per_m, "permeability_H_per_m")
    radial = float(magnetic_flux_density_radial_T)
    tangential = float(magnetic_flux_density_tangential_T)
    if not math.isfinite(radial) or not math.isfinite(tangential):
        raise ValueError("air-gap flux-density components must be finite")
    return radial * tangential / permeability


def time_average_air_gap_shear_stress(
    magnetic_flux_density_radial_phasor_T,
    magnetic_flux_density_tangential_phasor_T,
    permeability_H_per_m=MU0,
    amplitude="peak",
):
    """Return cycle-averaged air-gap shear stress from complex phasors."""

    permeability = _positive_scalar(permeability_H_per_m, "permeability_H_per_m")
    radial = complex(magnetic_flux_density_radial_phasor_T)
    tangential = complex(magnetic_flux_density_tangential_phasor_T)
    if not all(
        math.isfinite(value)
        for value in (radial.real, radial.imag, tangential.real, tangential.imag)
    ):
        raise ValueError("air-gap flux-density phasors must be finite")
    return (
        _phasor_average_factor(amplitude)
        * (radial * tangential.conjugate()).real
        / permeability
    )


def air_gap_shear_torque(
    magnetic_flux_density_radial_T,
    magnetic_flux_density_tangential_T,
    radius_m,
    axial_length_m=1.0,
    angle_rad=2.0 * math.pi,
    permeability_H_per_m=MU0,
):
    """Return torque in N m from uniform cylindrical air-gap shear stress."""

    radius = _nonnegative_scalar(radius_m, "radius_m")
    length = _nonnegative_scalar(axial_length_m, "axial_length_m")
    angle = _nonnegative_scalar(angle_rad, "angle_rad")
    shear = air_gap_shear_stress(
        magnetic_flux_density_radial_T,
        magnetic_flux_density_tangential_T,
        permeability_H_per_m=permeability_H_per_m,
    )
    return shear * radius * radius * angle * length


def _integrate_periodic_angle_samples(
    angles_rad: np.ndarray,
    values: np.ndarray,
    *,
    periodic: bool,
    period_rad: float,
) -> tuple[float, list[dict[str, float]]]:
    """Trapezoid-integrate scalar angle samples and retain segment evidence."""

    _strictly_increasing(angles_rad, "angles_rad")
    period = float(period_rad)
    if not math.isfinite(period):
        raise ValueError("period_rad must be finite")
    if periodic and period <= 0.0:
        raise ValueError("period_rad must be > 0 for periodic integration")
    if periodic and angles_rad[-1] - angles_rad[0] >= period:
        raise ValueError(
            "periodic angles must omit the duplicate endpoint and span less than period_rad"
        )
    segment_count = angles_rad.size if periodic else angles_rad.size - 1
    integral = 0.0
    rows: list[dict[str, float]] = []
    for index in range(segment_count):
        next_index = (index + 1) % angles_rad.size
        angle_start = float(angles_rad[index])
        angle_end = float(angles_rad[next_index])
        if periodic and next_index == 0:
            angle_end += period
        width = angle_end - angle_start
        if width <= 0.0:
            raise ValueError("angle segment width must be > 0")
        average = 0.5 * float(values[index] + values[next_index])
        integral += average * width
        rows.append(
            {
                "segment_index": index + 1,
                "angle_start_rad": angle_start,
                "angle_end_rad": angle_end,
                "angle_width_rad": width,
                "value_start": float(values[index]),
                "value_end": float(values[next_index]),
                "value_average": average,
            }
        )
    return integral, rows


def air_gap_shear_torque_from_angle_samples(
    angles_rad,
    magnetic_flux_density_radial_T,
    magnetic_flux_density_tangential_T,
    radius_m,
    axial_length_m=1.0,
    periodic=True,
    period_rad=2.0 * math.pi,
    permeability_H_per_m=MU0,
):
    """Integrate sampled cylindrical air-gap shear into torque.

    The samples omit a duplicate endpoint for a periodic contour. The returned
    JSON-friendly summary includes the trapezoid segments, tangential force,
    and torque, so sector scaling can be audited rather than hidden.
    """

    angles = _real_table(angles_rad, "angles_rad", minimum_size=2)
    radial = _real_table(
        magnetic_flux_density_radial_T,
        "magnetic_flux_density_radial_T",
        minimum_size=2,
    )
    tangential = _real_table(
        magnetic_flux_density_tangential_T,
        "magnetic_flux_density_tangential_T",
        minimum_size=2,
    )
    if angles.shape != radial.shape or angles.shape != tangential.shape:
        raise ValueError("angles_rad, radial B, and tangential B must have the same length")
    radius = _nonnegative_scalar(radius_m, "radius_m")
    length = _nonnegative_scalar(axial_length_m, "axial_length_m")
    permeability = _positive_scalar(permeability_H_per_m, "permeability_H_per_m")
    shear = radial * tangential / permeability
    integral, rows = _integrate_periodic_angle_samples(
        angles,
        shear,
        periodic=bool(periodic),
        period_rad=period_rad,
    )
    for index, row in enumerate(rows):
        next_index = (index + 1) % angles.size
        row.update(
            {
                "B_radial_start_T": float(radial[index]),
                "B_radial_end_T": float(radial[next_index]),
                "B_tangential_start_T": float(tangential[index]),
                "B_tangential_end_T": float(tangential[next_index]),
                "shear_start_Pa": row.pop("value_start"),
                "shear_end_Pa": row.pop("value_end"),
                "shear_average_Pa": row.pop("value_average"),
            }
        )
        row["tangential_force_N"] = (
            row["shear_average_Pa"] * radius * length * row["angle_width_rad"]
        )
        row["torque_Nm"] = row["tangential_force_N"] * radius
    integrated_angle = sum(row["angle_width_rad"] for row in rows)
    force = radius * length * integral
    torque = radius * force
    return {
        "n_samples": int(angles.size),
        "n_segments": len(rows),
        "periodic": bool(periodic),
        "period_rad": float(period_rad),
        "radius_m": radius,
        "axial_length_m": length,
        "permeability_H_per_m": permeability,
        "integrated_angle_rad": integrated_angle,
        "integral_shear_dtheta_Pa_rad": integral,
        "average_shear_stress_Pa": integral / integrated_angle,
        "tangential_force_N": force,
        "torque_Nm": torque,
        "torque_per_axial_length_N": torque / length if length > 0.0 else math.inf,
        "rows": rows,
    }


def time_average_air_gap_shear_torque_from_angle_samples(
    angles_rad,
    magnetic_flux_density_radial_phasor_T,
    magnetic_flux_density_tangential_phasor_T,
    radius_m,
    axial_length_m=1.0,
    periodic=True,
    period_rad=2.0 * math.pi,
    permeability_H_per_m=MU0,
    amplitude="peak",
):
    """Integrate cycle-averaged phasor air-gap shear into torque."""

    angles = _real_table(angles_rad, "angles_rad", minimum_size=2)
    radial = np.asarray(magnetic_flux_density_radial_phasor_T, dtype=complex)
    tangential = np.asarray(magnetic_flux_density_tangential_phasor_T, dtype=complex)
    if radial.ndim != 1 or tangential.ndim != 1:
        raise ValueError("radial and tangential B phasors must be one-dimensional tables")
    if angles.shape != radial.shape or angles.shape != tangential.shape:
        raise ValueError("angles_rad, radial B, and tangential B must have the same length")
    if not (
        np.all(np.isfinite(radial.real))
        and np.all(np.isfinite(radial.imag))
        and np.all(np.isfinite(tangential.real))
        and np.all(np.isfinite(tangential.imag))
    ):
        raise ValueError("air-gap flux-density phasors must be finite")
    permeability = _positive_scalar(permeability_H_per_m, "permeability_H_per_m")
    factor = _phasor_average_factor(amplitude)
    shear = factor * np.real(radial * np.conjugate(tangential)) / permeability
    radius = _nonnegative_scalar(radius_m, "radius_m")
    length = _nonnegative_scalar(axial_length_m, "axial_length_m")
    integral, rows = _integrate_periodic_angle_samples(
        angles,
        shear,
        periodic=bool(periodic),
        period_rad=period_rad,
    )
    for index, row in enumerate(rows):
        next_index = (index + 1) % angles.size
        row.update(
            {
                "B_radial_start_real_T": float(radial[index].real),
                "B_radial_start_imag_T": float(radial[index].imag),
                "B_radial_end_real_T": float(radial[next_index].real),
                "B_radial_end_imag_T": float(radial[next_index].imag),
                "B_tangential_start_real_T": float(tangential[index].real),
                "B_tangential_start_imag_T": float(tangential[index].imag),
                "B_tangential_end_real_T": float(tangential[next_index].real),
                "B_tangential_end_imag_T": float(tangential[next_index].imag),
                "shear_start_Pa": row.pop("value_start"),
                "shear_end_Pa": row.pop("value_end"),
                "shear_average_Pa": row.pop("value_average"),
            }
        )
        row["tangential_force_N"] = (
            row["shear_average_Pa"] * radius * length * row["angle_width_rad"]
        )
        row["torque_Nm"] = row["tangential_force_N"] * radius
    integrated_angle = sum(row["angle_width_rad"] for row in rows)
    force = radius * length * integral
    torque = radius * force
    return {
        "n_samples": int(angles.size),
        "n_segments": len(rows),
        "periodic": bool(periodic),
        "period_rad": float(period_rad),
        "radius_m": radius,
        "axial_length_m": length,
        "permeability_H_per_m": permeability,
        "phasor_amplitude": str(amplitude).strip().lower(),
        "integrated_angle_rad": integrated_angle,
        "integral_shear_dtheta_Pa_rad": integral,
        "average_shear_stress_Pa": integral / integrated_angle,
        "tangential_force_N": force,
        "torque_Nm": torque,
        "torque_per_axial_length_N": torque / length if length > 0.0 else math.inf,
        "rows": rows,
    }


def virtual_work_force_from_displacement_samples(
    positions_m,
    energy_J,
    energy_kind="coenergy",
):
    """Differentiate an energy table into force samples in N.

    Fixed-current coenergy uses ``F=dW'/dx``; fixed-flux stored energy uses
    ``F=-dW/dx``.  Interior points use centred differences and endpoints use
    one-sided differences, matching the legacy radia-ngsolve force contract.
    """

    positions = _real_table(positions_m, "positions_m", minimum_size=3)
    energy = _real_table(energy_J, "energy_J", minimum_size=3)
    if positions.shape != energy.shape:
        raise ValueError("positions_m and energy_J must have the same length")
    _strictly_increasing(positions, "positions_m")
    derivative = np.empty_like(energy)
    derivative[0] = (energy[1] - energy[0]) / (positions[1] - positions[0])
    derivative[-1] = (energy[-1] - energy[-2]) / (positions[-1] - positions[-2])
    derivative[1:-1] = (energy[2:] - energy[:-2]) / (positions[2:] - positions[:-2])
    return _energy_sign(energy_kind) * derivative


def coenergy_torque_from_angle_samples(
    angles_rad,
    coenergy_J,
    periodic=False,
    period_rad=2.0 * math.pi,
):
    """Differentiate fixed-current coenergy into torque samples in N m."""

    angles = _real_table(angles_rad, "angles_rad", minimum_size=3)
    coenergy = _real_table(coenergy_J, "coenergy_J", minimum_size=3)
    if angles.shape != coenergy.shape:
        raise ValueError("angles_rad and coenergy_J must have the same length")
    _strictly_increasing(angles, "angles_rad")
    torque = np.empty_like(coenergy)
    if periodic:
        period = _positive_scalar(period_rad, "period_rad")
        count = angles.size
        for index in range(count):
            minus = (index - 1) % count
            plus = (index + 1) % count
            angle_minus = angles[minus] - (period if minus > index else 0.0)
            angle_plus = angles[plus] + (period if plus < index else 0.0)
            torque[index] = (
                (coenergy[plus] - coenergy[minus]) / (angle_plus - angle_minus)
            )
    else:
        torque[0] = (coenergy[1] - coenergy[0]) / (angles[1] - angles[0])
        torque[-1] = (coenergy[-1] - coenergy[-2]) / (angles[-1] - angles[-2])
        torque[1:-1] = (
            (coenergy[2:] - coenergy[:-2]) / (angles[2:] - angles[:-2])
        )
    return torque


def force_torque_result(
    force_N,
    torque_Nm=None,
    *,
    method,
    frame="global_cartesian",
    pivot_m=None,
    field_convention="static",
    amplitude=None,
    dimensionality="3d",
    per_unit_depth=False,
):
    """Build the shared JSON-serializable force/torque result contract."""

    force = None
    if force_N is not None:
        force_array = _real_vectors(force_N, "force_N")
        if force_array.ndim != 1:
            raise ValueError("force_N must be one Cartesian three-vector")
        force = force_array.tolist()
    torque = None
    if torque_Nm is not None:
        torque_array = _real_vectors(torque_Nm, "torque_Nm")
        if torque_array.ndim != 1:
            raise ValueError("torque_Nm must be one Cartesian three-vector")
        torque = torque_array.tolist()
    if force is None and torque is None:
        raise ValueError("at least one of force_N or torque_Nm must be provided")
    pivot = None
    if pivot_m is not None:
        pivot_array = _real_vectors(pivot_m, "pivot_m")
        if pivot_array.ndim != 1:
            raise ValueError("pivot_m must be one Cartesian three-vector")
        pivot = pivot_array.tolist()
    convention = str(field_convention).strip().lower()
    if convention not in {"static", "time_average_phasor"}:
        raise ValueError("field_convention must be 'static' or 'time_average_phasor'")
    if convention == "time_average_phasor":
        _phasor_average_factor(amplitude)
    elif amplitude is not None:
        raise ValueError("amplitude applies only to time_average_phasor results")
    dimension = str(dimensionality).strip().lower()
    if dimension not in {"3d", "2d_planar", "axisymmetric"}:
        raise ValueError("dimensionality must be '3d', '2d_planar', or 'axisymmetric'")
    method_name = str(method).strip()
    frame_name = str(frame).strip()
    if not method_name:
        raise ValueError("method must be nonempty")
    if not frame_name:
        raise ValueError("frame must be nonempty")
    return {
        "schema": "radia.force-result/v1",
        "method": method_name,
        "frame": frame_name,
        "dimensionality": dimension,
        "per_unit_depth": bool(per_unit_depth),
        "field_convention": convention,
        "phasor_amplitude": amplitude if convention == "time_average_phasor" else None,
        "pivot_m": pivot,
        "force_N": force,
        "torque_Nm": torque,
        "units": {
            "force": "N/m" if per_unit_depth else "N",
            "torque": "N" if per_unit_depth else "N m",
            "pivot": "m",
        },
    }
