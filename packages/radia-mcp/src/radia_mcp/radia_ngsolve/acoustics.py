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


def spherical_hankel2(degree, argument):
    r"""Spherical Hankel function ``h_l^(2)(z)`` for outgoing waves.

    The module uses the ``exp(+i omega t)`` convention, so outgoing scalar
    Helmholtz waves are represented by ``h_l^(2)(k r)``.  A small recurrence is
    enough for low-order educational FEM/BEM checks and avoids a SciPy
    dependency in the public helper.
    """

    ell = int(degree)
    if ell < 0:
        raise ValueError("degree must be >= 0")
    z = complex(argument)
    if z == 0.0:
        raise ValueError("argument must be nonzero")

    h0 = 1j * cmath.exp(-1j * z) / z
    if ell == 0:
        return h0

    h1 = -cmath.exp(-1j * z) * (1.0 / z - 1j / (z * z))
    if ell == 1:
        return h1

    prev, current = h0, h1
    for n in range(1, ell):
        nxt = (2 * n + 1) * current / z - prev
        prev, current = current, nxt
    return current


def spherical_helmholtz_dtn_eigenvalue(radius, wavenumber, degree):
    r"""Exterior spherical Helmholtz DtN eigenvalue for one angular degree.

    For a pressure trace ``p(a) Y_l^m`` on a sphere of radius ``a``, the
    outgoing exterior field is proportional to ``h_l^(2)(k r) Y_l^m``.  The
    exact Dirichlet-to-Neumann eigenvalue is

        lambda_l = (partial_r p / p)|_{r=a}
                 = k h_l^(2)'(k a) / h_l^(2)(k a).

    It is a compact analytic gate for acoustic FEM/BEM coupling: the FEM trace
    supplies pressure on the sphere, while the exterior BEM or radiation
    condition supplies the normal derivative.  The outward normal is the
    increasing-radius direction.
    """

    a = float(radius)
    if a <= 0.0:
        raise ValueError("radius must be > 0")
    ell = int(degree)
    if ell < 0:
        raise ValueError("degree must be >= 0")
    k = complex(wavenumber)
    z = k * a
    if z == 0.0:
        raise ValueError("wavenumber * radius must be nonzero")

    if ell == 0:
        derivative_ratio = -1j - 1.0 / z
    else:
        h_l = spherical_hankel2(ell, z)
        derivative_ratio = spherical_hankel2(ell - 1, z) / h_l - (ell + 1) / z
    return k * derivative_ratio


def spherical_mode_radiation_impedance(
    radius,
    frequency,
    degree,
    rho=1.2041,
    c=343.0,
):
    r"""Radiation impedance of one outgoing spherical acoustic mode.

    The returned ``specific_impedance`` is ``p / v_n`` on the spherical
    boundary for one ``Y_l^m`` pressure/normal-velocity mode.  With
    ``exp(+i omega t)``, Euler's equation gives

        v_n = i * lambda_l * p / (omega rho),
        p / v_n = -i * omega rho / lambda_l.

    Degree zero matches :func:`pulsating_sphere_radiation`.  Higher degrees are
    useful as readable FEM/BEM gates because each spherical-harmonic trace mode
    has an exact exterior DtN value.
    """

    a = float(radius)
    f = float(frequency)
    rrho = float(rho)
    cc = float(c)
    ell = int(degree)
    if a <= 0.0:
        raise ValueError("radius must be > 0")
    if f <= 0.0:
        raise ValueError("frequency must be > 0")
    if rrho <= 0.0:
        raise ValueError("rho must be > 0")
    if cc <= 0.0:
        raise ValueError("c must be > 0")
    if ell < 0:
        raise ValueError("degree must be >= 0")

    omega = 2.0 * math.pi * f
    k = omega / cc
    dtn = spherical_helmholtz_dtn_eigenvalue(a, k, ell)
    z_specific = -1j * omega * rrho / dtn
    return {
        "radius": a,
        "frequency": f,
        "omega": omega,
        "wavenumber": k,
        "ka": k * a,
        "degree": ell,
        "rho": rrho,
        "c": cc,
        "dtn_eigenvalue": dtn,
        "specific_impedance": z_specific,
        "radiation_efficiency": z_specific.real / (rrho * cc),
        "reactance_ratio": z_specific.imag / (rrho * cc),
    }


def planar_helmholtz_dtn_symbol(wavenumber, tangential_wavenumber=0.0):
    r"""Exterior half-space Helmholtz DtN symbol for one planar trace mode.

    On a flat boundary with outward normal ``n`` into the exterior half-space,
    a pressure trace Fourier mode with tangential wavenumber ``k_t`` has an
    outgoing/decaying exterior field

        p(x_t, n) = p_0 exp(i k_t x_t) exp(-i q n),

    where ``q^2 = k^2 - k_t^2``.  The Dirichlet-to-Neumann symbol is therefore

        lambda(k_t) = partial_n p / p = -i q.

    The square-root branch is chosen so propagating modes have ``Re(q) >= 0``
    and evanescent modes decay into the exterior (``Im(q) <= 0``).  This is the
    planar analogue of :func:`spherical_helmholtz_dtn_eigenvalue` and is a tiny
    readable FEM/BEM coupling gate: FEM pressure trace in, exterior normal
    derivative out.
    """

    k = complex(wavenumber)
    kt = float(tangential_wavenumber)
    if abs(k) <= 0.0:
        raise ValueError("wavenumber must be nonzero")
    if kt < 0.0:
        raise ValueError("tangential_wavenumber must be >= 0")

    q = cmath.sqrt(k * k - kt * kt)
    if q.imag > 0.0 or (abs(q.imag) <= 1.0e-15 and q.real < 0.0):
        q = -q

    if abs(k.imag) <= 1.0e-15 and k.real > 0.0:
        if abs(kt - k.real) <= 1.0e-14 * max(1.0, k.real):
            regime = "grazing"
        elif kt < k.real:
            regime = "propagating"
        else:
            regime = "evanescent"
    else:
        regime = "complex"

    dtn = -1j * q
    return {
        "wavenumber": k,
        "tangential_wavenumber": kt,
        "normal_wavenumber": q,
        "dtn_eigenvalue": dtn,
        "symbol_identity_residual": dtn * dtn - (kt * kt - k * k),
        "regime": regime,
    }


def planar_mode_radiation_impedance(
    frequency,
    tangential_wavenumber=None,
    incidence_angle_rad=None,
    rho=1.2041,
    c=343.0,
):
    r"""Specific acoustic impedance for a planar outgoing exterior mode.

    Exactly one of ``tangential_wavenumber`` or ``incidence_angle_rad`` must be
    provided.  For a propagating plane wave at angle ``theta`` from the normal,

        z_n = p / v_n = rho c / cos(theta),

    while evanescent modes have zero active radiation resistance and a purely
    reactive normal impedance.  The return dictionary includes the matching DtN
    symbol so FEM/BEM sign conventions can be checked in one place.
    """

    f = float(frequency)
    rrho = float(rho)
    cc = float(c)
    if f <= 0.0:
        raise ValueError("frequency must be > 0")
    if rrho <= 0.0:
        raise ValueError("rho must be > 0")
    if cc <= 0.0:
        raise ValueError("c must be > 0")
    if (tangential_wavenumber is None) == (incidence_angle_rad is None):
        raise ValueError("provide exactly one of tangential_wavenumber or incidence_angle_rad")

    omega = 2.0 * math.pi * f
    k = omega / cc
    angle = None
    if incidence_angle_rad is not None:
        angle = float(incidence_angle_rad)
        if abs(angle) >= 0.5 * math.pi:
            raise ValueError("incidence_angle_rad must be strictly between -pi/2 and pi/2")
        kt = k * abs(math.sin(angle))
    else:
        kt = float(tangential_wavenumber)

    symbol = planar_helmholtz_dtn_symbol(k, kt)
    q = symbol["normal_wavenumber"]
    if abs(q) <= 0.0:
        raise ValueError("grazing modes have infinite normal impedance")

    z_specific = omega * rrho / q
    return {
        "frequency": f,
        "omega": omega,
        "wavenumber": k,
        "tangential_wavenumber": kt,
        "incidence_angle_rad": angle,
        "rho": rrho,
        "c": cc,
        "regime": symbol["regime"],
        "normal_wavenumber": q,
        "dtn_eigenvalue": symbol["dtn_eigenvalue"],
        "specific_impedance": z_specific,
        "normalized_impedance": z_specific / (rrho * cc),
        "radiation_efficiency": z_specific.real / (rrho * cc),
        "reactance_ratio": z_specific.imag / (rrho * cc),
    }


def acoustic_dtn_from_impedance(
    frequency,
    specific_impedance=None,
    specific_admittance=None,
    rho=1.2041,
):
    r"""Convert acoustic impedance/admittance to a Helmholtz DtN coefficient.

    With the module's ``exp(+i omega t)`` convention, Euler's equation gives

        v_n = i (partial_n p) / (omega rho).

    For a boundary specific impedance ``z = p/v_n`` or admittance ``Y=v_n/p``,
    the equivalent scalar Helmholtz Robin/DtN coefficient is

        partial_n p = lambda p,     lambda = -i omega rho / z = -i omega rho Y.

    This is the tiny conversion bridge between FEM impedance-boundary rows and
    exterior BEM/DtN operators.  Exactly one of ``specific_impedance`` or
    ``specific_admittance`` must be supplied.
    """

    f = float(frequency)
    rrho = float(rho)
    if f <= 0.0:
        raise ValueError("frequency must be > 0")
    if rrho <= 0.0:
        raise ValueError("rho must be > 0")
    if (specific_impedance is None) == (specific_admittance is None):
        raise ValueError("provide exactly one of specific_impedance or specific_admittance")

    omega = 2.0 * math.pi * f
    if specific_impedance is not None:
        z = complex(specific_impedance)
        if z == 0.0:
            raise ValueError("specific_impedance must be nonzero")
        y = 1.0 / z
    else:
        y = complex(specific_admittance)
    dtn = -1j * omega * rrho * y
    return {
        "frequency": f,
        "omega": omega,
        "rho": rrho,
        "specific_impedance": math.inf if y == 0.0 else 1.0 / y,
        "specific_admittance": y,
        "dtn_eigenvalue": dtn,
        "robin_coefficient": dtn,
    }


def acoustic_impedance_from_dtn(frequency, dtn_eigenvalue, rho=1.2041):
    r"""Convert a Helmholtz DtN/Robin coefficient to acoustic impedance.

    This is the inverse of :func:`acoustic_dtn_from_impedance`:

        z = -i omega rho / lambda,     Y = i lambda / (omega rho).

    It is useful for checking whether a boundary operator is active/radiating
    (``Re(z)>0``) or purely reactive.
    """

    f = float(frequency)
    rrho = float(rho)
    if f <= 0.0:
        raise ValueError("frequency must be > 0")
    if rrho <= 0.0:
        raise ValueError("rho must be > 0")
    lam = complex(dtn_eigenvalue)
    if lam == 0.0:
        raise ValueError("dtn_eigenvalue must be nonzero")
    omega = 2.0 * math.pi * f
    z = -1j * omega * rrho / lam
    return {
        "frequency": f,
        "omega": omega,
        "rho": rrho,
        "dtn_eigenvalue": lam,
        "specific_impedance": z,
        "specific_admittance": 1.0 / z,
    }


def baffled_circular_piston_radiation(
    radius,
    frequency,
    surface_velocity=1.0,
    rho=1.2041,
    c=343.0,
):
    r"""Radiation impedance of a uniformly vibrating circular piston in an infinite baffle.

    A flat circular piston of radius ``a`` with uniform normal velocity ``v0``
    is the canonical acoustic FEM/BEM boundary example for a baffled speaker,
    transducer, or duct opening.  The average specific radiation impedance is

        z / (rho c) = 1 - J_1(2ka)/(ka) + i H_1(2ka)/(ka),

    where ``J_1`` is a Bessel function and ``H_1`` is a Struve function.  At low
    frequency, the radiation resistance scales as ``(ka)^2/2`` and the
    reactance as ``8 ka/(3 pi)``; at high frequency the resistance tends to the
    plane-wave value ``rho c``.  Peak phasors are used, so active power is
    ``0.5 * area * Re(z) * |v0|^2``.
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

    from scipy.special import j1, struve

    velocity = complex(surface_velocity)
    omega = 2.0 * math.pi * f
    k = omega / cc
    ka = k * a
    area = math.pi * a * a
    resistance_ratio = 1.0 - float(j1(2.0 * ka)) / ka
    reactance_ratio = float(struve(1, 2.0 * ka)) / ka
    z_specific = rrho * cc * complex(resistance_ratio, reactance_ratio)
    volume_velocity = area * velocity
    z_volume_velocity = z_specific / area
    radiated_power = 0.5 * area * z_specific.real * abs(velocity) ** 2
    return {
        "radius": a,
        "frequency": f,
        "omega": omega,
        "wavenumber": k,
        "ka": ka,
        "rho": rrho,
        "c": cc,
        "surface_area": area,
        "surface_velocity": velocity,
        "volume_velocity": volume_velocity,
        "specific_impedance": z_specific,
        "specific_resistance": z_specific.real,
        "specific_reactance": z_specific.imag,
        "radiation_efficiency": resistance_ratio,
        "reactance_ratio": reactance_ratio,
        "volume_velocity_impedance": z_volume_velocity,
        "radiated_power": radiated_power,
        "low_ka_resistance_asymptote": 0.5 * ka * ka,
        "low_ka_reactance_asymptote": 8.0 * ka / (3.0 * math.pi),
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
