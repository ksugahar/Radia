r"""Closed-form acoustic radiation helpers for readable FEM/BEM validation.

The module starts with the canonical free-field exterior problem that every
acoustic FEM/BEM code should be able to explain before it can be trusted on
general geometries: a uniformly pulsating sphere in an infinite fluid.

For peak phasors with time factor ``exp(+i omega t)``, an outgoing spherical
wave is

    p(r) = A exp(-i k r) / r,

and Euler's equation gives the local specific acoustic impedance

    z(r) = p / v_r = rho c / (1 - i/(k r)).

At the sphere surface ``r=a`` with prescribed radial velocity ``v_a``, this
gives the radiation resistance and reactance directly:

    Re(z_a)/(rho c) = (ka)^2 / (1 + (ka)^2),
    Im(z_a)/(rho c) =  ka    / (1 + (ka)^2).

These are pure analytic helpers. They are useful as low-frequency BEM gates:
the resistance scales as ``(ka)^2`` while the reactive near field scales as
``ka``, so cancellation and sign conventions show up immediately.
"""

from __future__ import annotations

import cmath
import math


def helmholtz_green_3d(distance, wavenumber):
    r"""Outgoing 3D scalar Helmholtz Green function.

    With the module's ``exp(+i omega t)`` convention, the outgoing free-space
    kernel is

        G_k(r) = exp(-i k r) / (4 pi r).

    This is the point-source kernel behind acoustic single-layer BEM.  Use
    :func:`helmholtz_green_low_frequency_series` when studying the low-frequency
    split into the singular Laplace kernel plus smooth corrections.
    """

    r = float(distance)
    if r <= 0.0:
        raise ValueError("distance must be > 0")
    k = complex(wavenumber)
    return cmath.exp(-1j * k * r) / (4.0 * math.pi * r)


def helmholtz_green_low_frequency_series(distance, wavenumber, order=6):
    r"""Low-frequency series/split of the outgoing 3D Helmholtz Green function.

    The expansion

        exp(-i k r)/(4 pi r)
          = 1/(4 pi r) - i k/(4 pi) - k^2 r/(8 pi)
            + i k^3 r^2/(24 pi) + ...

    cleanly separates the static Laplace singularity from a smooth regular
    remainder.  That split is the readable low-frequency BEM gate: the singular
    quadrature is the same as electrostatics, while the frequency-dependent
    corrections are regular panel integrals.

    Returns a dictionary with the complex ``terms`` through ``order``, the
    ``laplace_term`` (n=0), the ``regular_part`` (n>=1), ``approx``, ``exact``,
    and absolute error.  ``order`` is the highest Taylor index retained.
    """

    r = float(distance)
    if r <= 0.0:
        raise ValueError("distance must be > 0")
    nmax = int(order)
    if nmax < 0:
        raise ValueError("order must be >= 0")
    k = complex(wavenumber)
    terms = []
    for n in range(nmax + 1):
        term = ((-1j * k) ** n) * (r ** (n - 1)) / (4.0 * math.pi * math.factorial(n))
        terms.append(term)
    approx = sum(terms)
    exact = helmholtz_green_3d(r, k)
    return {
        "distance": r,
        "wavenumber": k,
        "order": nmax,
        "kr_abs": abs(k * r),
        "terms": terms,
        "laplace_term": terms[0],
        "regular_part": sum(terms[1:]) if len(terms) > 1 else 0.0j,
        "approx": approx,
        "exact": exact,
        "abs_error": abs(approx - exact),
    }


def _specific_spherical_impedance(k_radius, rho, c):
    kr = float(k_radius)
    if kr <= 0.0:
        raise ValueError("k_radius must be > 0")
    denom = 1.0 + kr * kr
    return rho * c * complex(kr * kr / denom, kr / denom)


def pulsating_sphere_radiation(
    radius,
    frequency,
    surface_velocity,
    rho=1.2041,
    c=343.0,
    sample_radius=None,
):
    r"""Radiation of a uniformly pulsating sphere in an infinite fluid.

    Parameters
    ----------
    radius : float
        Sphere radius ``a`` [m], > 0.
    frequency : float
        Frequency [Hz], > 0.
    surface_velocity : complex
        Peak radial surface velocity phasor ``v_a`` [m/s].
    rho : float, default 1.2041
        Fluid density [kg/m^3].
    c : float, default 343.0
        Sound speed [m/s].
    sample_radius : float, optional
        Radius where the outgoing pressure and radial velocity are reported.
        Defaults to ``10*a``. Must be >= ``a``.

    Returns
    -------
    dict
        Frequencies, ``ka``, surface impedance, volume-velocity impedance,
        radiated active power, and one exact spherical-wave sample point.

    Notes
    -----
    Peak phasors are used, so active power is
    ``0.5 * Re(p * conj(v))`` integrated over the sphere. If your solver uses
    RMS phasors, omit the factor of 0.5 when comparing powers.
    """

    a = float(radius)
    f = float(frequency)
    rrho = float(rho)
    cc = float(c)
    if a <= 0.0:
        raise ValueError("radius must be > 0")
    if f <= 0.0:
        raise ValueError("frequency must be > 0")
    if rrho <= 0.0:
        raise ValueError("rho must be > 0")
    if cc <= 0.0:
        raise ValueError("c must be > 0")

    r = float(10.0 * a if sample_radius is None else sample_radius)
    if r < a:
        raise ValueError("sample_radius must be >= radius")

    v_surface = complex(surface_velocity)
    omega = 2.0 * math.pi * f
    k = omega / cc
    ka = k * a
    area = 4.0 * math.pi * a * a

    z_surface = _specific_spherical_impedance(ka, rrho, cc)
    p_surface = z_surface * v_surface
    volume_velocity = area * v_surface
    z_volume_velocity = z_surface / area
    radiated_power = 0.5 * area * z_surface.real * abs(v_surface) ** 2

    phase = cmath.exp(-1j * k * (r - a))
    p_sample = p_surface * (a / r) * phase
    z_sample = _specific_spherical_impedance(k * r, rrho, cc)
    v_sample = p_sample / z_sample
    intensity_sample = 0.5 * (p_sample * v_sample.conjugate()).real
    power_from_sample = 4.0 * math.pi * r * r * intensity_sample
    plane_wave_intensity_sample = abs(p_sample) ** 2 / (2.0 * rrho * cc)

    return {
        "radius": a,
        "frequency": f,
        "omega": omega,
        "wavenumber": k,
        "ka": ka,
        "rho": rrho,
        "c": cc,
        "surface_area": area,
        "surface_velocity": v_surface,
        "volume_velocity": volume_velocity,
        "specific_impedance": z_surface,
        "specific_resistance": z_surface.real,
        "specific_reactance": z_surface.imag,
        "radiation_efficiency": z_surface.real / (rrho * cc),
        "reactance_ratio": z_surface.imag / (rrho * cc),
        "volume_velocity_impedance": z_volume_velocity,
        "radiated_power": radiated_power,
        "sample_radius": r,
        "sample_pressure": p_sample,
        "sample_radial_velocity": v_sample,
        "sample_specific_impedance": z_sample,
        "sample_intensity": intensity_sample,
        "sample_power": power_from_sample,
        "sample_plane_wave_intensity": plane_wave_intensity_sample,
    }
