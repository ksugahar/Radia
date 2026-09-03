"""Readable analytic references for plane-wave scattering by spheres.

These partial-wave series deliberately remain independent of Radia's native
kernels. They use SciPy's special functions and serve as validation oracles
for the NGSolve FEM/BEM acoustic application lane. The outgoing-wave
convention is ``exp(+i k r)``; fluid sound speed and density are normalized to
one.
"""

from __future__ import annotations

import numpy as np
from scipy.special import eval_legendre, spherical_jn, spherical_yn

_MAXIMUM_TERMS = 512


def _positive(value, name):
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def _term_count(value, *, automatic=False):
    if value is None and automatic:
        return -1
    count = int(value)
    if automatic and count == -1:
        return count
    if count < 0 or count > _MAXIMUM_TERMS:
        suffix = " or -1 for automatic selection" if automatic else ""
        raise ValueError(f"terms must lie between 0 and {_MAXIMUM_TERMS}{suffix}")
    return count


def _jn(order, value):
    return spherical_jn(order, value)


def _jn_d(order, value):
    return spherical_jn(order, value, derivative=True)


def _jn_dd(order, value):
    return -(2.0 / value) * _jn_d(order, value) - (
        1.0 - order * (order + 1) / value**2
    ) * _jn(order, value)


def _h1(order, value):
    return spherical_jn(order, value) + 1j * spherical_yn(order, value)


def _h1_d(order, value):
    return spherical_jn(order, value, derivative=True) + 1j * spherical_yn(
        order, value, derivative=True
    )


def _prepare_points(points, radius, *, allow_interior=False):
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    if len(points) == 0 or not np.all(np.isfinite(points)):
        raise ValueError("points must be a nonempty N-by-3 array of finite values")
    distance = np.linalg.norm(points, axis=1)
    if not allow_interior and np.any(distance < radius * (1.0 - 1.0e-9)):
        raise ValueError(
            "evaluation points must lie on or outside the sphere r >= R"
        )
    safe_distance = np.maximum(distance, 1.0e-30)
    cosine = points[:, 2] / safe_distance
    return points, distance, safe_distance, cosine


def _incident(wavenumber, points):
    return np.exp(1j * wavenumber * points[:, 2])


def soft_sphere_scattering(wavenumber, radius, points, terms=None):
    """Return the sound-soft sphere partial-wave reference."""
    wavenumber = _positive(wavenumber, "wavenumber")
    radius = _positive(radius, "radius")
    requested = _term_count(terms, automatic=True)
    points, distance, _, cosine = _prepare_points(points, radius)
    count = int(np.ceil(wavenumber * radius)) + 12 if requested < 0 else requested
    _term_count(count)

    scattered = np.zeros(distance.shape, dtype=complex)
    last = np.zeros(distance.shape, dtype=complex)
    for order in range(count + 1):
        coefficient = (
            -(1j**order)
            * (2 * order + 1)
            * _jn(order, wavenumber * radius)
            / _h1(order, wavenumber * radius)
        )
        last = (
            coefficient
            * _h1(order, wavenumber * distance)
            * eval_legendre(order, cosine)
        )
        scattered += last
    incident = _incident(wavenumber, points)
    return {
        "backend": "scipy-reference",
        "kind": "soft_sphere_plane_wave_scattering_series",
        "wavenumber": wavenumber,
        "radius": radius,
        "terms": count,
        "truncation_tail": float(np.max(np.abs(last))),
        "scattered": scattered,
        "incident": incident,
        "total": incident + scattered,
    }


def rigid_sphere_scattering(wavenumber, radius, points, terms=1):
    """Return the sound-hard sphere partial-wave reference."""
    wavenumber = _positive(wavenumber, "wavenumber")
    radius = _positive(radius, "radius")
    requested = _term_count(terms)
    points, distance, _, cosine = _prepare_points(points, radius)
    count = max(
        requested,
        int(np.ceil(wavenumber * max(radius, float(np.max(distance))))) + 12,
    )
    _term_count(count)

    scattered = np.zeros(distance.shape, dtype=complex)
    last = np.zeros(distance.shape, dtype=complex)
    boundary = wavenumber * radius
    for order in range(count + 1):
        coefficient = (
            -(1j**order)
            * (2 * order + 1)
            * _jn_d(order, boundary)
            / _h1_d(order, boundary)
        )
        last = (
            coefficient
            * _h1(order, wavenumber * distance)
            * eval_legendre(order, cosine)
        )
        scattered += last
    incident = _incident(wavenumber, points)
    return {
        "backend": "scipy-reference",
        "kind": "rigid_sphere_plane_wave_scattering_series",
        "wavenumber": wavenumber,
        "radius": radius,
        "terms": count,
        "truncation_tail": float(np.max(np.abs(last))),
        "scattered": scattered,
        "incident": incident,
        "total": incident + scattered,
    }


def fluid_sphere_scattering(
    wavenumber,
    radius,
    points,
    interior_wavenumber=None,
    density_ratio=1.0,
    terms=None,
):
    """Return the Anderson penetrable-fluid-sphere reference."""
    wavenumber = _positive(wavenumber, "wavenumber")
    radius = _positive(radius, "radius")
    interior_wavenumber = _positive(
        wavenumber if interior_wavenumber is None else interior_wavenumber,
        "interior_wavenumber",
    )
    density_ratio = _positive(density_ratio, "density_ratio")
    requested = _term_count(terms, automatic=True)
    points, distance, safe_distance, cosine = _prepare_points(
        points, radius, allow_interior=True
    )
    count = max(
        0 if requested < 0 else requested,
        int(
            np.ceil(
                max(
                    wavenumber * max(radius, float(np.max(distance))),
                    interior_wavenumber * radius,
                )
            )
        )
        + 12,
    )
    _term_count(count)

    inside = distance <= radius * (1.0 + 1.0e-12)
    x0 = wavenumber * radius
    x1 = interior_wavenumber * radius
    total = np.zeros(distance.shape, dtype=complex)
    last = np.zeros(distance.shape, dtype=complex)
    for order in range(count + 1):
        legendre = eval_legendre(order, cosine)
        incident_coefficient = (1j**order) * (2 * order + 1)
        j0 = _jn(order, x0)
        h0 = _h1(order, x0)
        j1 = _jn(order, x1)
        beta = (
            (interior_wavenumber / density_ratio)
            * _jn_d(order, x1)
            / j1
        )
        scattered_coefficient = -incident_coefficient * (
            wavenumber * _jn_d(order, x0) - beta * j0
        ) / (wavenumber * _h1_d(order, x0) - beta * h0)
        interior_coefficient = (
            incident_coefficient * j0 + scattered_coefficient * h0
        ) / j1

        mode = np.zeros(distance.shape, dtype=complex)
        mode[inside] = (
            interior_coefficient
            * _jn(order, interior_wavenumber * safe_distance[inside])
            * legendre[inside]
        )
        mode[~inside] = (
            incident_coefficient
            * _jn(order, wavenumber * safe_distance[~inside])
            + scattered_coefficient
            * _h1(order, wavenumber * safe_distance[~inside])
        ) * legendre[~inside]
        total += mode
        last = mode

    return {
        "backend": "scipy-reference",
        "kind": "fluid_sphere_transmission_scattering_series",
        "wavenumber": wavenumber,
        "interior_wavenumber": interior_wavenumber,
        "density_ratio": density_ratio,
        "radius": radius,
        "terms": count,
        "truncation_tail": float(np.max(np.abs(last))),
        "incident": _incident(wavenumber, points),
        "total": total,
        "inside_mask": inside,
    }


def _elastic_coefficient(
    order,
    wavenumber,
    radius,
    longitudinal_speed,
    shear_speed,
    density_ratio,
):
    omega = wavenumber
    k_longitudinal = omega / longitudinal_speed
    x = wavenumber * radius
    xl = k_longitudinal * radius
    fluid_factor = wavenumber / omega**2
    mu = density_ratio * shear_speed**2
    lame_lambda = density_ratio * (
        longitudinal_speed**2 - 2.0 * shear_speed**2
    )

    if shear_speed == 0.0:
        matrix = np.array(
            [
                [
                    fluid_factor * _h1_d(order, x),
                    -k_longitudinal * _jn_d(order, xl),
                ],
                [
                    _h1(order, x),
                    -lame_lambda * k_longitudinal**2 * _jn(order, xl),
                ],
            ],
            dtype=complex,
        )
        rhs = np.array(
            [-fluid_factor * _jn_d(order, x), -_jn(order, x)], dtype=complex
        )
        return np.linalg.solve(matrix, rhs)[0]

    k_transverse = omega / shear_speed
    xt = k_transverse * radius
    angular = order * (order + 1)
    ur_a = k_longitudinal * _jn_d(order, xl)
    ur_b = angular / radius * _jn(order, xt)
    dur_a = k_longitudinal**2 * _jn_dd(order, xl)
    dur_b = angular * (
        k_transverse * _jn_d(order, xt) / radius
        - _jn(order, xt) / radius**2
    )
    srr_a = (
        -lame_lambda * k_longitudinal**2 * _jn(order, xl) + 2.0 * mu * dur_a
    )
    srr_b = 2.0 * mu * dur_b
    va_a = _jn(order, xl) / radius
    va_b = _jn(order, xt) / radius + k_transverse * _jn_d(order, xt)
    vp_a = -_jn(order, xl) / radius**2 + k_longitudinal * _jn_d(
        order, xl
    ) / radius
    vp_b = (
        -_jn(order, xt) / radius**2
        + k_transverse * _jn_d(order, xt) / radius
        + k_transverse**2 * _jn_dd(order, xt)
    )
    srt_a = mu * (ur_a / radius + vp_a - va_a / radius)
    srt_b = mu * (ur_b / radius + vp_b - va_b / radius)
    matrix = np.array(
        [
            [fluid_factor * _h1_d(order, x), -ur_a, -ur_b],
            [_h1(order, x), srr_a, srr_b],
            [0.0, srt_a, srt_b],
        ],
        dtype=complex,
    )
    rhs = np.array(
        [-fluid_factor * _jn_d(order, x), -_jn(order, x), 0.0], dtype=complex
    )
    return np.linalg.solve(matrix, rhs)[0]


def elastic_sphere_scattering(
    wavenumber,
    radius,
    points,
    longitudinal_speed=2.0,
    shear_speed=1.0,
    density_ratio=1.5,
    terms=0,
):
    """Return the Faran elastic-solid-sphere reference."""
    wavenumber = _positive(wavenumber, "wavenumber")
    radius = _positive(radius, "radius")
    longitudinal_speed = _positive(longitudinal_speed, "longitudinal_speed")
    shear_speed = float(shear_speed)
    if not np.isfinite(shear_speed) or shear_speed < 0.0:
        raise ValueError("shear_speed must be nonnegative and finite")
    density_ratio = _positive(density_ratio, "density_ratio")
    requested = _term_count(terms)
    points, distance, _, cosine = _prepare_points(points, radius)
    count = (
        requested
        if requested > 0
        else int(np.ceil(wavenumber * radius)) + 10
    )
    _term_count(count)

    coefficients = [
        _elastic_coefficient(
            order,
            wavenumber,
            radius,
            longitudinal_speed,
            shear_speed,
            density_ratio,
        )
        for order in range(count + 1)
    ]
    scattered = np.zeros(distance.shape, dtype=complex)
    last = np.zeros(distance.shape, dtype=complex)
    for order, coefficient in enumerate(coefficients):
        last = (
            (1j**order)
            * (2 * order + 1)
            * coefficient
            * _h1(order, wavenumber * distance)
            * eval_legendre(order, cosine)
        )
        scattered += last
    incident = _incident(wavenumber, points)
    return {
        "backend": "scipy-reference",
        "kind": "elastic_solid_sphere_faran_scattering_series",
        "wavenumber": wavenumber,
        "radius": radius,
        "longitudinal_speed": longitudinal_speed,
        "shear_speed": shear_speed,
        "density_ratio": density_ratio,
        "terms": count,
        "truncation_tail": float(np.max(np.abs(last))),
        "incident": incident,
        "scattered": scattered,
        "total": incident + scattered,
    }
