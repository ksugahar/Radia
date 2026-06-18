r"""Antenna directivity / radiation-resistance / array-factor / aperture closed forms
for radia-ngsolve -- the radiation (far-field) counterpart of the static and
cavity solvers in this package.

Where the cavity/waveguide modules solve interior eigenproblems, an antenna
RADIATES into unbounded space: the relevant figures of merit are closed forms in
the canonical electromagnetics texts. This module collects the textbook
expressions for the elementary radiators and apertures:

    * the Hertzian (short) dipole and the small magnetic loop -- the two dual
      elementary sources, each with directivity D0 = 3/2;
    * the resonant half-wave dipole input impedance (Z = 73.08 + j42.5 Ohm) and
      its directivity D0 = 4/Cin(2 pi) = 1.6409;
    * the finite-length centre-fed dipole pattern and its directivity by direct
      numerical integration of the radiation intensity;
    * the uniform N-element linear array factor (with progressive scan phase),
      the discrete analogue of the continuous aperture; and
    * the uniformly illuminated circular aperture (Airy) pattern, with its
      2 J1(u)/u far field, first-null beamwidth and aperture directivity.

All quantities are normalized to wavelength (length / lambda, circumference /
lambda, diameter / lambda), so the formulas are frequency-independent. These are
pure analytic helpers (no field solve); the underlying derivations are standard.

References:
    C. A. Balanis, "Antenna Theory: Analysis and Design", 4th ed.,
        Ch. 4 (linear wire antennas, short and finite dipole),
        Ch. 5 (loop antennas), Ch. 6 (arrays), Ch. 12 (aperture antennas).
    J. D. Kraus, "Antennas", 2nd ed. (radiation resistance, directivity).
    S. Ramo, J. R. Whinnery, T. Van Duzer, "Fields and Waves in
        Communication Electronics", 3rd ed. (far-field radiation).
    J. D. Jackson, "Classical Electrodynamics", 3rd ed. (multipole radiation).
"""
import cmath
import math

EPS0 = 8.8541878128e-12      # vacuum permittivity [F/m]
MU0 = 4.0e-7 * math.pi       # vacuum permeability [H/m]
C0 = 299792458.0             # speed of light in vacuum [m/s]
ETA0 = 376.730313668         # intrinsic impedance of free space [Ohm] (= MU0*C0)

# Cosine-integral constant Cin(2 pi) = gamma + ln(2 pi) - Ci(2 pi); it sets both
# the half-wave dipole radiation resistance and its directivity (Balanis Ch.4).
CIN_2PI = 2.4376576957


def short_dipole_radiation_resistance(length_over_lambda):
    r"""Radiation resistance of a Hertzian (short, uniform-current) dipole.

    For a dipole of length L << lambda carrying a uniform current, the radiated
    power gives (Balanis Ch.4; Kraus):

        R_rad = 20 * pi^2 * (L/lambda)^2   [Ohm].

    The directivity is independent of length, D0 = 3/2 (the sin^2(theta) pattern).

    Parameters
    ----------
    length_over_lambda : float
        Electrical length L/lambda (> 0; the short-dipole result assumes << 1).

    Returns
    -------
    dict
        ``R_rad`` [Ohm] and ``directivity`` (= 1.5).
    """
    L = float(length_over_lambda)
    if L <= 0.0:
        raise ValueError("length_over_lambda must be > 0")
    R = 20.0 * math.pi ** 2 * L ** 2
    return {"R_rad": R, "directivity": 1.5}


def small_loop_radiation_resistance(circumference_over_lambda, n_turns=1):
    r"""Radiation resistance of a small (electrically small) magnetic loop.

    For a loop of circumference C << lambda with ``n_turns`` turns the radiated
    power gives (Balanis Ch.5):

        R_rad = 20 * pi^2 * (C/lambda)^4 * n^2   [Ohm].

    The fourth-power dependence (vs. the dipole's square) is the hallmark of a
    magnetic-dipole radiator. The directivity equals that of the elementary
    electric dipole, D0 = 3/2 (the loop is its magnetic dual).

    Parameters
    ----------
    circumference_over_lambda : float
        Electrical circumference C/lambda (> 0; small-loop result assumes << 1).
    n_turns : int
        Number of turns (>= 1).

    Returns
    -------
    dict
        ``R_rad`` [Ohm] and ``directivity`` (= 1.5).
    """
    C = float(circumference_over_lambda)
    n = int(n_turns)
    if C <= 0.0:
        raise ValueError("circumference_over_lambda must be > 0")
    if n < 1:
        raise ValueError("n_turns must be >= 1")
    R = 20.0 * math.pi ** 2 * C ** 4 * n ** 2
    return {"R_rad": R, "directivity": 1.5}


def half_wave_dipole_impedance():
    r"""Input impedance and directivity of the resonant half-wave dipole.

    For the centre-fed L = lambda/2 dipole the radiation resistance follows from
    the cosine integral (Balanis Ch.4; Kraus):

        R_rad = (eta0 / (4 pi)) * Cin(2 pi) ~= 73.08 Ohm,

    with Cin(2 pi) = 2.4376576957. The (thin-wire) reactance is X ~= +42.5 Ohm,
    so Z ~= 73.08 + j 42.5 Ohm; the slightly shorter resonant length removes the
    reactance. The directivity is

        D0 = 4 / Cin(2 pi) = 1.640922   (2.151 dBi).

    Returns
    -------
    dict
        ``R_rad`` [Ohm], ``X`` [Ohm], ``Z`` (complex impedance) and
        ``directivity``.
    """
    R = (ETA0 / (4.0 * math.pi)) * CIN_2PI
    X = 42.5
    return {
        "R_rad": R,
        "X": X,
        "Z": complex(R, X),
        "directivity": 4.0 / CIN_2PI,
    }


def finite_dipole_directivity(length_over_lambda):
    r"""Directivity of a finite-length centre-fed (sinusoidal-current) dipole.

    The far-field pattern of a centre-fed dipole of length L is (Balanis Ch.4):

        F(theta) = [cos((pi L/lambda) cos theta) - cos(pi L/lambda)] / sin theta,

    and the directivity is

        D0 = 2 * F_max^2 / INT_0^pi F(theta)^2 sin theta d theta,

    evaluated here by midpoint numerical integration. For L/lambda = 0.5 this
    reproduces the half-wave-dipole directivity D0 = 1.6409.

    Parameters
    ----------
    length_over_lambda : float
        Electrical length L/lambda (> 0).

    Returns
    -------
    dict
        ``directivity`` and ``directivity_dBi`` (10 log10 D0).
    """
    L = float(length_over_lambda)
    if L <= 0.0:
        raise ValueError("length_over_lambda must be > 0")
    beta_l = math.pi * L  # (pi L/lambda); = k L / 2 with k = 2 pi/lambda

    def f(theta):
        s = math.sin(theta)
        if s == 0.0:
            return 0.0
        return (math.cos(beta_l * math.cos(theta)) - math.cos(beta_l)) / s

    n = 2000
    integral = 0.0
    fmax2 = 0.0
    for i in range(n):
        theta = (i + 0.5) * math.pi / n           # midpoint in (0, pi)
        val = f(theta)
        fmax2 = max(fmax2, val * val)
        integral += val * val * math.sin(theta)
    integral *= math.pi / n
    if integral <= 0.0 or fmax2 <= 0.0:
        raise ValueError("degenerate pattern (zero radiated power)")
    D0 = 2.0 * fmax2 / integral
    return {"directivity": D0, "directivity_dBi": 10.0 * math.log10(D0)}


def uniform_linear_array_factor(n_elements, spacing_over_lambda, theta_deg,
                                scan_deg=0.0):
    r"""Array factor of a uniform N-element linear array (progressive scan phase).

    For N equal-amplitude elements spaced d apart, with a progressive phase that
    steers the beam to ``scan_deg``, the electrical phase per element is
    (Balanis Ch.6):

        psi = 2 pi (d/lambda) (cos theta - cos theta_scan),

    and the (un-normalized) array factor magnitude is

        |AF| = |sin(N psi/2) / sin(psi/2)|     -> N   as psi -> 0   (L'Hopital),

    with the normalized factor AF/N (unity at the main beam). theta is measured
    from the array axis.

    Parameters
    ----------
    n_elements : int
        Number of elements (>= 1).
    spacing_over_lambda : float
        Element spacing d/lambda (> 0).
    theta_deg : float
        Observation angle from the array axis [deg].
    scan_deg : float
        Beam scan (main-beam) angle from the array axis [deg].

    Returns
    -------
    dict
        ``AF`` (un-normalized), ``AF_normalized`` (= AF/N) and ``psi`` [rad].
    """
    N = int(n_elements)
    d = float(spacing_over_lambda)
    if N < 1:
        raise ValueError("n_elements must be >= 1")
    if d <= 0.0:
        raise ValueError("spacing_over_lambda must be > 0")
    theta = math.radians(theta_deg)
    scan = math.radians(scan_deg)
    psi = 2.0 * math.pi * d * (math.cos(theta) - math.cos(scan))
    half = psi / 2.0
    s = math.sin(half)
    if abs(s) < 1e-12:                              # L'Hopital limit psi -> 0
        AF = float(N)
    else:
        AF = abs(math.sin(N * half) / s)
    return {"AF": AF, "AF_normalized": AF / N, "psi": psi}


def _bessel_j1(x):
    try:
        from scipy.special import j1
    except ModuleNotFoundError:
        pass
    else:
        return float(j1(x))

    # Power series: J1(x) = sum_k (-1)^k (x/2)^(2k+1)/(k!(k+1)!).
    term = 0.5 * x
    total = term
    x2_over_4 = 0.25 * x * x
    for k in range(1, 80):
        term *= -x2_over_4 / (k * (k + 1))
        total += term
        if abs(term) <= 1e-16 * max(1.0, abs(total)):
            break
    return total


def circular_aperture_pattern(diameter_over_lambda, theta_deg):
    r"""Far-field pattern of a uniformly illuminated circular aperture (Airy).

    For a circular aperture of diameter D with uniform illumination the
    normalized far field is (Balanis Ch.12):

        E(theta) = 2 J1(u) / u ,   u = pi (D/lambda) sin theta,

    which tends to 1 as u -> 0 and has its first null at u = 3.8317 (the first
    zero of J1). The uniform-aperture directivity is

        D0 = (pi D/lambda)^2 = (4 pi / lambda^2) * Area.

    Parameters
    ----------
    diameter_over_lambda : float
        Electrical diameter D/lambda (> 0).
    theta_deg : float
        Observation angle from the aperture boresight [deg].

    Returns
    -------
    dict
        ``E_normalized`` (= 2 J1(u)/u), ``u``, ``directivity_uniform`` and
        ``first_null_deg`` (the first-null angle, or ``None`` if the null is
        beyond 90 deg for this aperture).
    """
    D = float(diameter_over_lambda)
    if D <= 0.0:
        raise ValueError("diameter_over_lambda must be > 0")

    theta = math.radians(theta_deg)
    u = math.pi * D * math.sin(theta)
    if abs(u) < 1e-12:
        E = 1.0
    else:
        E = 2.0 * _bessel_j1(u) / u

    # first null: u = 3.8317  ->  sin(theta_null) = 3.8317 / (pi D/lambda)
    U1 = 3.8317059702
    s_null = U1 / (math.pi * D)
    first_null_deg = math.degrees(math.asin(s_null)) if s_null <= 1.0 else None

    return {
        "E_normalized": E,
        "u": u,
        "directivity_uniform": (math.pi * D) ** 2,
        "first_null_deg": first_null_deg,
    }
