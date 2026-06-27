"""Waveguide / cavity cutoff modes for radia-ngsolve -- the 2D Helmholtz
(transverse Laplacian) EIGENVALUE problem, the wave-physics counterpart of the
static BVP solvers.

A hollow metallic waveguide supports propagating modes whose transverse pattern
solves the scalar Helmholtz eigenproblem on the cross-section:

    -nabla_t^2 psi = k_c^2 psi ,

with the boundary condition set by the mode family (PEC walls):

    * TM modes : E_z = 0 on the wall  ->  DIRICHLET Laplacian
    * TE modes : dH_z/dn = 0 on the wall  ->  NEUMANN  Laplacian

The eigenvalues k_c^2 are the squared cutoff wavenumbers; below the cutoff
frequency f_c = c k_c/(2 pi) the mode is evanescent (does not propagate). The
SAME eigenproblem gives the TM resonances of a 2D cavity and the modes of a
vibrating membrane (drum) -- only the physical reading of k_c changes.

The Neumann problem always has a trivial constant (k_c = 0) eigenmode (no field);
:func:`helmholtz_cutoff_wavenumbers_2d` discards it and returns the lowest
PHYSICAL cutoff wavenumbers. Validated against the exact rectangular-guide
spectrum.
"""
import cmath
import math

from ngsolve import (H1, HCurl, BilinearForm, GridFunction, grad, curl, dx, ArnoldiSolver)

C0 = 299792458.0          # speed of light in vacuum [m/s]
MU0 = 4.0e-7 * math.pi     # vacuum permeability [H/m]  (eta0 = MU0*C0 = 376.73 Ohm)


def rectangular_waveguide_cutoff(a, b, m, n):
    """Exact cutoff FREQUENCY [Hz] of the TE_mn / TM_mn mode of a rectangular metallic
    waveguide of inner width ``a`` and height ``b`` [m]:

        f_c,mn = (c/2) sqrt((m/a)^2 + (n/b)^2),

    i.e. k_c = sqrt((m pi/a)^2 + (n pi/b)^2). Allowed indices: TE needs m,n >= 0 (not both 0,
    dominant TE10 -> f_c = c/(2a)); TM needs m,n >= 1 (lowest TM11). The empty band between the
    dominant TE10 and the next mode sets the single-mode bandwidth (e.g. WR-90: 6.56 -> 13.1 GHz).
    """
    return 0.5 * C0 * math.hypot(m / a, n / b)


def rectangular_waveguide_mode_table(a, b, max_m=3, max_n=3, families=("TE", "TM"), c=C0):
    """Sorted cutoff-mode table for a PEC rectangular waveguide.

    ``TE`` modes allow ``m,n >= 0`` except ``TE00``; ``TM`` modes require
    ``m,n >= 1``.  Each returned row contains ``family``, ``mode``, ``m``,
    ``n``, ``cutoff_frequency`` [Hz], and ``cutoff_wavenumber`` [1/m].  The
    table is intentionally small and transparent: it is the analytic checklist
    for deciding whether a port solve is below cutoff, single-mode, or already
    multi-mode.
    """
    aa, bb = float(a), float(b)
    if aa <= 0.0 or bb <= 0.0:
        raise ValueError("waveguide dimensions must be positive")
    mm, nn = int(max_m), int(max_n)
    if mm < 0 or nn < 0:
        raise ValueError("max_m and max_n must be >= 0")

    selected = tuple(str(f).upper() for f in families)
    unknown = [f for f in selected if f not in ("TE", "TM")]
    if unknown:
        raise ValueError("families may contain only 'TE' and 'TM'")

    modes = []
    for family in selected:
        for m in range(mm + 1):
            for n in range(nn + 1):
                if family == "TE":
                    if m == 0 and n == 0:
                        continue
                elif m == 0 or n == 0:
                    continue
                fc = 0.5 * c * math.hypot(m / aa, n / bb)
                modes.append({
                    "family": family,
                    "mode": f"{family}{m}{n}",
                    "m": m,
                    "n": n,
                    "cutoff_frequency": fc,
                    "cutoff_wavenumber": 2.0 * math.pi * fc / c,
                })
    modes.sort(key=lambda row: (
        row["cutoff_frequency"],
        0 if row["family"] == "TE" else 1,
        row["m"],
        row["n"],
    ))
    return modes


def rectangular_waveguide_band_summary(a, b, frequency, max_m=3, max_n=3, c=C0):
    """Classify a rectangular guide at one frequency from its cutoff table.

    Modes with ``frequency > cutoff_frequency`` are listed as propagating; the
    next higher cutoff is reported for margin checks.  ``single_mode`` is true
    only when the sole propagating mode is the dominant ``TE10``.  This mirrors
    the practical RF-port workflow without needing any solver-specific output.
    """
    f = float(frequency)
    if f <= 0.0:
        raise ValueError("frequency must be positive")
    table = rectangular_waveguide_mode_table(a, b, max_m=max_m, max_n=max_n, c=c)
    propagating = [row for row in table if f > row["cutoff_frequency"]]
    next_rows = [row for row in table if f <= row["cutoff_frequency"]]
    dominant = table[0] if table else None
    next_mode = next_rows[0] if next_rows else None
    single_mode = (
        len(propagating) == 1
        and propagating[0]["family"] == "TE"
        and propagating[0]["m"] == 1
        and propagating[0]["n"] == 0
    )
    return {
        "frequency": f,
        "dominant_mode": dominant,
        "propagating_modes": propagating,
        "n_propagating": len(propagating),
        "next_mode": next_mode,
        "below_dominant_cutoff": len(propagating) == 0,
        "single_mode": single_mode,
        "single_mode_upper_cutoff": table[1]["cutoff_frequency"] if len(table) > 1 else None,
    }


def cutoff_frequency(kc, c=C0):
    """Cutoff frequency f_c = c k_c/(2 pi) [Hz] from a cutoff wavenumber k_c [1/m]."""
    return c * kc / (2.0 * math.pi)


def rectangular_cavity_modes(a, b, d, max_index=3, c=C0, limit=None):
    """Sorted closed-form mode table for a PEC rectangular cavity.

    Enumerates index triples ``(m, n, p)`` with at least two non-zero indices and
    frequency

        f_mnp = (c/2) sqrt((m/a)^2 + (n/b)^2 + (p/d)^2).

    The "at least two" rule keeps the physical TE/TM cavity family used by the
    public regression tests (e.g. TE101, TM110, TE011) while excluding the
    one-dimensional gradient/null family that appears as spurious near-zero
    modes in shifted HCurl eigenvalue solves. Returns dictionaries sorted by
    ``frequency`` then by indices; ``limit`` optionally truncates the list.
    """
    aa, bb, dd = float(a), float(b), float(d)
    if aa <= 0.0 or bb <= 0.0 or dd <= 0.0:
        raise ValueError("cavity dimensions must be positive")
    mi = int(max_index)
    if mi < 1:
        raise ValueError("max_index must be >= 1")
    if limit is not None and int(limit) < 1:
        raise ValueError("limit must be >= 1 when supplied")

    modes = []
    for m in range(mi + 1):
        for n in range(mi + 1):
            for p in range(mi + 1):
                if (m > 0) + (n > 0) + (p > 0) < 2:
                    continue
                f = rectangular_cavity_frequency(aa, bb, dd, m, n, p, c)
                modes.append({"m": m, "n": n, "p": p,
                              "indices": (m, n, p), "frequency": f})
    modes.sort(key=lambda row: (row["frequency"], row["m"], row["n"], row["p"]))
    return modes if limit is None else modes[:int(limit)]


def circular_waveguide_cutoff(radius, mode, m, n):
    """Exact cutoff FREQUENCY [Hz] of the TE_mn / TM_mn mode of a CIRCULAR metallic waveguide of
    inner ``radius`` a -- the Bessel-function counterpart of the rectangular
    :func:`rectangular_waveguide_cutoff`:

        TM_mn: k_c = j_mn / a    (j_mn  = n-th positive zero of J_m),
        TE_mn: k_c = j'_mn / a   (j'_mn = n-th positive zero of J'_m),    f_c = c k_c/(2 pi).

    The dominant mode is TE11 (j'_11 = 1.8412 -> f_c = 0.293 c/a); the next is TM01 (j_01 = 2.4048).
    Circular guides have the famous degeneracy TE_0n / TM_1n (because j'_0n = j_1n). ``mode`` is
    'TE' or 'TM'; m = azimuthal index (>=0), n = radial index (>=1)."""
    from scipy.special import jn_zeros, jnp_zeros
    if mode.upper() == "TE":
        z = jnp_zeros(m, n)[n - 1]
    elif mode.upper() == "TM":
        z = jn_zeros(m, n)[n - 1]
    else:
        raise ValueError("mode must be 'TE' or 'TM'")
    return C0 * z / (2.0 * math.pi * radius)


def circular_waveguide_mode_table(radius, max_m=3, max_n=2, families=("TE", "TM"), c=C0):
    """Sorted cutoff-mode table for a PEC circular waveguide.

    ``TE`` modes use zeros of ``J'_m`` and ``TM`` modes use zeros of ``J_m``.  The azimuthal
    degeneracy is reported explicitly: modes with ``m=0`` have one angular pattern, while
    ``m>0`` has the usual twofold ``cos(m phi)`` / ``sin(m phi)`` degeneracy.  Each row contains
    ``family``, ``mode``, ``m``, ``n``, ``cutoff_frequency`` [Hz], ``cutoff_wavenumber`` [1/m],
    and ``angular_degeneracy``.
    """
    r = float(radius)
    if r <= 0.0:
        raise ValueError("radius must be positive")
    mm, nn = int(max_m), int(max_n)
    if mm < 0:
        raise ValueError("max_m must be >= 0")
    if nn < 1:
        raise ValueError("max_n must be >= 1")

    selected = tuple(str(f).upper() for f in families)
    unknown = [f for f in selected if f not in ("TE", "TM")]
    if unknown:
        raise ValueError("families may contain only 'TE' and 'TM'")

    modes = []
    for family in selected:
        for m in range(mm + 1):
            for n in range(1, nn + 1):
                fc = circular_waveguide_cutoff(r, family, m, n)
                modes.append({
                    "family": family,
                    "mode": f"{family}{m}{n}",
                    "m": m,
                    "n": n,
                    "cutoff_frequency": fc,
                    "cutoff_wavenumber": 2.0 * math.pi * fc / c,
                    "angular_degeneracy": 1 if m == 0 else 2,
                })
    modes.sort(key=lambda row: (
        row["cutoff_frequency"],
        0 if row["family"] == "TE" else 1,
        row["m"],
        row["n"],
    ))
    return modes


def circular_waveguide_band_summary(radius, frequency, max_m=3, max_n=2, c=C0):
    """Classify a circular guide at one frequency from its cutoff table.

    Modes with ``frequency > cutoff_frequency`` are listed as propagating.  ``single_mode`` is
    true when the only propagating row is the dominant ``TE11`` pair.  The row still reports
    ``angular_degeneracy=2``, so callers can distinguish the two polarizations from the next
    physical cutoff.
    """
    f = float(frequency)
    if f <= 0.0:
        raise ValueError("frequency must be positive")
    table = circular_waveguide_mode_table(radius, max_m=max_m, max_n=max_n, c=c)
    propagating = [row for row in table if f > row["cutoff_frequency"]]
    next_rows = [row for row in table if f <= row["cutoff_frequency"]]
    dominant = table[0] if table else None
    next_mode = next_rows[0] if next_rows else None
    single_mode = (
        len(propagating) == 1
        and propagating[0]["family"] == "TE"
        and propagating[0]["m"] == 1
        and propagating[0]["n"] == 1
    )
    return {
        "frequency": f,
        "dominant_mode": dominant,
        "propagating_modes": propagating,
        "n_propagating_rows": len(propagating),
        "n_propagating_with_degeneracy": sum(row["angular_degeneracy"] for row in propagating),
        "next_mode": next_mode,
        "below_dominant_cutoff": len(propagating) == 0,
        "single_mode": single_mode,
        "single_mode_upper_cutoff": table[1]["cutoff_frequency"] if len(table) > 1 else None,
    }


def waveguide_dispersion(frequency, fc, c=C0):
    """Dispersion of a hollow metallic waveguide mode ABOVE cutoff -- the propagation sequel to the
    cutoff eigenvalue (:func:`rectangular_waveguide_cutoff` / :func:`circular_waveguide_cutoff` /
    :func:`helmholtz_cutoff_wavenumbers_2d`). A mode of cutoff frequency ``fc`` [Hz] driven at
    ``frequency`` f > fc propagates with axial constant

        beta = (2 pi/c) sqrt(f^2 - fc^2) = k0 sqrt(1 - (fc/f)^2),    k0 = 2 pi f/c,

    so the GUIDE WAVELENGTH lambda_g = 2 pi/beta = lambda0 / sqrt(1-(fc/f)^2) exceeds the free-space
    lambda0 = c/f and obeys  1/lambda_g^2 = 1/lambda0^2 - 1/lambda_c^2. The PHASE velocity
    v_p = omega/beta = c/sqrt(1-(fc/f)^2) > c carries no energy; the GROUP (energy/signal) velocity
    v_g = d omega/d beta = c sqrt(1-(fc/f)^2) < c does, and the two satisfy the reciprocal identity

        v_p * v_g = c^2 .

    At cutoff (f -> fc+) beta -> 0, lambda_g -> inf, v_p -> inf, v_g -> 0 (standing wave, no axial
    transport); for f >> fc both velocities -> c (TEM-like). Below cutoff the mode is evanescent --
    see :func:`waveguide_evanescent_attenuation`.

    Returns a dict: ``beta`` [1/m], ``lambda_g`` [m], ``v_phase`` & ``v_group`` [m/s],
    ``fc_over_f``. Raises ValueError at or below cutoff (use the evanescent branch there)."""
    if frequency <= fc:
        raise ValueError("frequency must exceed the cutoff fc (mode is evanescent below it)")
    ratio = fc / frequency
    s = math.sqrt(1.0 - ratio * ratio)
    beta = (2.0 * math.pi * frequency / c) * s
    return {"beta": beta, "lambda_g": 2.0 * math.pi / beta,
            "v_phase": c / s, "v_group": c * s, "fc_over_f": ratio}


def guide_wavelength(frequency, fc, c=C0):
    """Guide wavelength lambda_g = lambda0/sqrt(1-(fc/f)^2) [m] (always > the free-space lambda0
    = c/f). The transverse standing-wave pattern stretches the axial period; component lengths
    (slots, irises, lambda_g/4 transformers) are set by lambda_g, not lambda0."""
    return waveguide_dispersion(frequency, fc, c)["lambda_g"]


def waveguide_wave_impedance(frequency, fc, mode="TE", eta=MU0 * C0):
    """Wave impedance [Ohm] of a propagating hollow-guide mode above cutoff.

    For a mode with cutoff ``fc`` driven at ``frequency > fc`` and intrinsic
    medium impedance ``eta``:

        Z_TE = eta / sqrt(1 - (fc/f)^2)
        Z_TM = eta * sqrt(1 - (fc/f)^2)

    Thus TE impedance tends to infinity at cutoff while TM impedance tends to
    zero; both tend to ``eta`` far above cutoff, and ``Z_TE * Z_TM = eta^2``.
    These impedances set the port normalization and the interface reflection
    coefficient used by dielectric-slab / cascade S-parameter helpers.
    """
    f = float(frequency)
    if f <= fc:
        raise ValueError("frequency must exceed cutoff for propagating wave impedance")
    s = math.sqrt(1.0 - (fc / f) ** 2)
    m = str(mode).upper()
    if m == "TE":
        z = eta / s
    elif m == "TM":
        z = eta * s
    else:
        raise ValueError("mode must be 'TE' or 'TM'")
    return {"mode": m, "frequency": f, "fc": fc, "eta": eta, "Z": z, "fc_over_f": fc / f}


def rectangular_waveguide_te10_port_normalization(frequency, width_a, height_b,
                                                  power_w=1.0, c=C0):
    """Peak field amplitudes for a power-normalised rectangular-guide TE10 port.

    With the usual lossless TE10 field convention

        E_y = E0 sin(pi x/a),   H_x = E_y / Z_TE,

    the time-average power through the cross-section is

        P = (1/2) integral(E_y H_x) dA = a b E0^2 / (4 Z_TE).

    This helper turns a requested port power into the corresponding peak
    electric and magnetic field amplitudes.  It is intentionally a closed-form
    port-normalisation checklist: the returned ``poynting_power_W`` recomputes
    the power integral from the amplitudes so examples/tests can catch unit or
    peak/RMS mistakes before a full-wave solve is trusted.
    """
    f = float(frequency)
    a = float(width_a)
    b = float(height_b)
    pwr = float(power_w)
    if f <= 0.0:
        raise ValueError("frequency must be positive")
    if a <= 0.0 or b <= 0.0:
        raise ValueError("width_a and height_b must be positive")
    if pwr < 0.0:
        raise ValueError("power_w must be non-negative")

    fc = 0.5 * c / a
    disp = waveguide_dispersion(f, fc, c)
    z_te = waveguide_wave_impedance(f, fc, "TE", eta=MU0 * c)["Z"]
    sin2_area = 0.5 * a * b
    e0 = math.sqrt(2.0 * pwr * z_te / sin2_area) if pwr else 0.0
    hx0 = e0 / z_te
    kc = math.pi / a
    omega = 2.0 * math.pi * f
    hz0 = e0 * kc / (omega * MU0)
    poynting = 0.5 * e0 * hx0 * sin2_area
    return {
        "frequency": f,
        "width_a": a,
        "height_b": b,
        "power_w": pwr,
        "fc": fc,
        "beta": disp["beta"],
        "lambda_g": disp["lambda_g"],
        "v_group": disp["v_group"],
        "Z_TE_ohm": z_te,
        "sin2_area_integral_m2": sin2_area,
        "E_y_peak_V_per_m": e0,
        "H_x_peak_A_per_m": hx0,
        "H_z_wall_peak_A_per_m": hz0,
        "H_z_over_H_x_peak": hz0 / hx0 if hx0 else math.inf,
        "poynting_power_W": poynting,
        "poynting_abs_error_W": abs(poynting - pwr),
    }


def rectangular_waveguide_te10_conductor_loss(frequency, width_a, height_b, sigma,
                                              length=None, mu_r=1.0, c=C0):
    """Conductor-loss attenuation of the TE10 mode in a rectangular metal waveguide.

    For an air-filled PEC-shaped guide with finite-conductivity walls, the good-conductor
    surface resistance ``R_s`` dissipates wall power ``(R_s/2)|H_t|^2``.  Integrating the
    TE10 fields over the four walls and normalising by transmitted power gives the
    amplitude attenuation constant

        alpha_c = R_s (2 b k_c^2 + a k_0^2) / (eta k_0 beta a b)   [Np/m],

    where ``a`` is width, ``b`` is height, ``k_c=pi/a``, ``k_0=2*pi*f/c``, and
    ``beta=sqrt(k_0^2-k_c^2)``.  The loss diverges near cutoff because group
    velocity and transmitted power collapse; far above cutoff it approaches the
    broad-wall skin-loss scale.  If ``length`` is supplied, the returned row also
    includes ``S21_mag=exp(-alpha_c length)`` and insertion loss in dB.
    """
    f = float(frequency)
    a = float(width_a)
    b = float(height_b)
    sig = float(sigma)
    if f <= 0.0:
        raise ValueError("frequency must be positive")
    if a <= 0.0 or b <= 0.0:
        raise ValueError("width_a and height_b must be positive")
    if sig <= 0.0:
        raise ValueError("sigma must be positive")
    if mu_r <= 0.0:
        raise ValueError("mu_r must be positive")

    fc = 0.5 * c / a
    if f <= fc:
        raise ValueError("frequency must exceed TE10 cutoff")
    k0 = 2.0 * math.pi * f / c
    kc = math.pi / a
    beta = math.sqrt(k0 * k0 - kc * kc)
    eta = MU0 * c
    rs = surface_resistance(f, sig, mu_r=mu_r)
    delta = skin_depth(f, sig, mu_r=mu_r)
    alpha = rs * (2.0 * b * kc * kc + a * k0 * k0) / (eta * k0 * beta * a * b)
    out = {
        "frequency": f,
        "width_a": a,
        "height_b": b,
        "sigma": sig,
        "mu_r": float(mu_r),
        "fc": fc,
        "k0": k0,
        "kc": kc,
        "beta": beta,
        "surface_resistance_ohm": rs,
        "skin_depth_m": delta,
        "alpha_np_per_m": alpha,
        "alpha_db_per_m": 20.0 * math.log10(math.e) * alpha,
    }
    if length is not None:
        ell = float(length)
        if ell < 0.0:
            raise ValueError("length must be non-negative")
        s21 = math.exp(-alpha * ell)
        out.update({
            "length_m": ell,
            "S21_mag": s21,
            "insertion_loss_db": -20.0 * math.log10(s21),
            "power_transmission_fraction": s21 * s21,
            "power_loss_fraction": 1.0 - s21 * s21,
        })
    return out


def waveguide_evanescent_attenuation(frequency, fc, c=C0):
    """BELOW cutoff (f < fc) the mode does NOT propagate; its amplitude decays as exp(-alpha z) with

        alpha = (2 pi/c) sqrt(fc^2 - f^2)   [Np/m]  (purely reactive -- lossless 'cutoff' attenuator).

    The waveguide is a high-pass filter: alpha -> 0 as f -> fc-, and alpha -> 2 pi fc/c = k_c (the
    cutoff wavenumber) as f -> 0. Used for cutoff attenuators and below-cutoff isolation. Raises
    ValueError at or above fc (use :func:`waveguide_dispersion` there)."""
    if frequency >= fc:
        raise ValueError("frequency must be below the cutoff fc (use waveguide_dispersion above it)")
    return (2.0 * math.pi / c) * math.sqrt(fc * fc - frequency * frequency)


def waveguide_dielectric_slab_sparams(frequency, width_a, eps_r, slab_length, c=C0):
    r"""The 2-port SCATTERING MATRIX (S11, S21) of a rectangular waveguide TE10 section loaded with a
    DIELECTRIC SLAB (relative permittivity ``eps_r``, axial length ``slab_length``, air on both sides) --
    the MISMATCHED-2-port sequel to the matched-line propagation (:func:`waveguide_dispersion`): a real
    reflection appears where the wave impedance steps.

    The TE10 transverse cutoff wavenumber ``k_c = pi/a`` is set by the width ``a`` alone, so the axial
    propagation constants in the air and slab regions are

        beta_air  = sqrt(k0^2 - k_c^2) = (2 pi/c) sqrt(f^2 - fc^2),   fc = c/(2a),
        beta_slab = sqrt(eps_r k0^2 - k_c^2),                          k0 = 2 pi f/c.

    The TE wave impedance Z_TE = omega mu0 / beta scales as 1/beta, so the air<->slab interface reflects
    with ``Gamma = (Z_slab - Z_air)/(Z_slab + Z_air) = (beta_air - beta_slab)/(beta_air + beta_slab)``.  For
    a finite slab the two interfaces + the in-slab phase ``theta = beta_slab * slab_length`` cascade to the
    exact transmission-line S-matrix (reference planes at the slab faces)

        S11 = Gamma (1 - e^{-2 i theta}) / (1 - Gamma^2 e^{-2 i theta}),
        S21 = (1 - Gamma^2) e^{-i theta} / (1 - Gamma^2 e^{-2 i theta}).

    LOSSLESS UNITARITY ``|S11|^2 + |S21|^2 = 1`` holds exactly.  Limits: slab_length -> 0 or eps_r -> 1
    give ``S11 = 0, |S21| = 1`` (matched); a QUARTER-WAVE slab (theta = pi/2) maximises the reflection at
    ``|S11| = 2|Gamma|/(1+Gamma^2)``.  Returns ``{frequency, fc, eps_r, slab_length, beta_air, beta_slab,
    gamma, theta, S11, S21, S11_mag, S21_mag, unitarity}``.  Textbook transmission-line theory; the closed
    form is the analytic gate a full-wave port (S-parameter) solve is checked against.  Raises ValueError
    at/below the TE10 cutoff."""
    f = float(frequency); a = float(width_a); er = float(eps_r); d = float(slab_length)
    if a <= 0.0:
        raise ValueError("width_a must be > 0 (got %r)" % (width_a,))
    if d < 0.0:
        raise ValueError("slab_length must be >= 0 (got %r)" % (slab_length,))
    if er <= 0.0:
        raise ValueError("eps_r must be > 0 (got %r)" % (eps_r,))
    fc = 0.5 * c / a                                            # TE10 cutoff (air)
    if f <= fc:
        raise ValueError("frequency must exceed the TE10 cutoff fc=c/2a (the mode is evanescent in air)")
    k0 = 2.0 * math.pi * f / c
    kc = math.pi / a                                           # transverse cutoff wavenumber
    beta_air = math.sqrt(k0 * k0 - kc * kc)
    arg = er * k0 * k0 - kc * kc
    if arg <= 0.0:                                             # slab below cutoff (only if eps_r small)
        raise ValueError("frequency below the slab cutoff fc/sqrt(eps_r); the slab mode is evanescent")
    beta_slab = math.sqrt(arg)
    gamma = (beta_air - beta_slab) / (beta_air + beta_slab)    # TE interface reflection (Z ~ 1/beta)
    theta = beta_slab * d
    e2 = cmath.exp(-2j * theta)
    den = 1.0 - gamma * gamma * e2
    S11 = gamma * (1.0 - e2) / den
    S21 = (1.0 - gamma * gamma) * cmath.exp(-1j * theta) / den
    return {"frequency": f, "fc": fc, "eps_r": er, "slab_length": d,
            "beta_air": beta_air, "beta_slab": beta_slab, "gamma": gamma, "theta": theta,
            "S11": S11, "S21": S21, "S11_mag": abs(S11), "S21_mag": abs(S21),
            "unitarity": abs(S11) ** 2 + abs(S21) ** 2}


def waveguide_cascade_sparams(frequency, width_a, sections, c=C0):
    r"""The 2-port SCATTERING MATRIX (S11, S21) of a rectangular waveguide TE10 line loaded with an
    ARBITRARY CASCADE of uniform dielectric sections -- the N-section generalisation of the single
    :func:`waveguide_dielectric_slab_sparams` (which is exactly the one-section case). The input and
    output ports are air-filled, so the reference impedance is the air-guide TE10 impedance.

    Each section ``(length_i, eps_r_i)`` is a uniform TE10 transmission line of axial constant
    ``beta_i = sqrt(eps_r_i k0^2 - k_c^2)`` (k_c = pi/a, k0 = 2 pi f/c) and wave impedance
    ``Z_i ~ 1/beta_i``; its ABCD (chain) matrix is

        [[cos theta_i,        j Z_i sin theta_i],
         [j sin theta_i/Z_i,  cos theta_i      ]],     theta_i = beta_i * length_i.

    The whole structure is the ORDERED PRODUCT of the section ABCDs, converted to S-parameters with
    the air reference impedance ``Z0 ~ 1/beta_air`` on both ports:

        S11 = (A + B/Z0 - C Z0 - D)/(A + B/Z0 + C Z0 + D),   S21 = 2/(A + B/Z0 + C Z0 + D).

    A section BELOW its own cutoff (eps_r_i k0^2 < k_c^2) is handled automatically as an EVANESCENT
    (below-cutoff) attenuator: ``beta_i`` becomes imaginary and the cos/sin turn into cosh/sinh. The
    ABCD entries are INVARIANT under the sqrt branch sign (cos is even; the B and C entries each pair
    sin with a compensating 1/beta or beta), so no branch bookkeeping is needed. The result is
    reciprocal (``A D - B C = 1``) and, for purely propagating lossless sections, unitary
    (|S11|^2 + |S21|^2 = 1). With a single section it reproduces
    :func:`waveguide_dielectric_slab_sparams` to machine precision (its internal consistency gate),
    and a short-terminated equivalent reproduces :func:`waveguide_offset_short_s11`.

    ``sections`` is a non-empty list of ``(length_m, eps_r)`` pairs (the air ports are implicit).
    Returns ``{frequency, fc, n_sections, S11, S21, S11_mag, S21_mag, unitarity, abcd_det}``.
    Textbook transmission-line (ABCD-cascade) theory; the closed form is the analytic gate a
    full-wave port solve of a MULTI-section guide is checked against. Raises ValueError at/below the
    air TE10 cutoff or for a malformed section list."""
    f = float(frequency); a = float(width_a)
    if a <= 0.0:
        raise ValueError("width_a must be > 0 (got %r)" % (width_a,))
    if not sections:
        raise ValueError("sections must be a non-empty list of (length, eps_r) pairs")
    fc = 0.5 * c / a                                           # air TE10 cutoff
    if f <= fc:
        raise ValueError("frequency must exceed the air TE10 cutoff fc=c/2a")
    k0 = 2.0 * math.pi * f / c
    kc = math.pi / a
    beta_air = math.sqrt(k0 * k0 - kc * kc)
    z0 = 1.0 / beta_air                                        # air-guide reference impedance (~1/beta)
    A, B, Cc, D = 1.0 + 0j, 0j, 0j, 1.0 + 0j                   # total chain matrix (start = identity)
    for (length_i, eps_i) in sections:
        li = float(length_i); eri = float(eps_i)
        if li < 0.0:
            raise ValueError("section length must be >= 0 (got %r)" % (length_i,))
        if eri <= 0.0:
            raise ValueError("section eps_r must be > 0 (got %r)" % (eps_i,))
        beta = cmath.sqrt(eri * k0 * k0 - kc * kc)             # complex => evanescent section
        zc = 1.0 / beta
        th = beta * li
        ct = cmath.cos(th); st = cmath.sin(th)
        a11 = ct; a12 = 1j * zc * st; a21 = 1j * st / zc; a22 = ct
        nA = A * a11 + B * a21                                 # [A B; C D] := [A B; C D] * section
        nB = A * a12 + B * a22
        nC = Cc * a11 + D * a21
        nD = Cc * a12 + D * a22
        A, B, Cc, D = nA, nB, nC, nD
    den = A + B / z0 + Cc * z0 + D
    S11 = (A + B / z0 - Cc * z0 - D) / den
    S21 = 2.0 / den
    return {"frequency": f, "fc": fc, "n_sections": len(sections),
            "S11": S11, "S21": S21, "S11_mag": abs(S11), "S21_mag": abs(S21),
            "unitarity": abs(S11) ** 2 + abs(S21) ** 2, "abcd_det": A * D - B * Cc}


def waveguide_offset_short_s11(frequency, width_a, offset_length, c=C0):
    r"""The 1-port reflection S11 of an OFFSET SHORT -- a length ``offset_length`` of air-filled
    rectangular TE10 waveguide terminated in a PEC short (the classic 1-port calibration standard;
    the reflection sibling of the antenna 1-port and the transmission slab 2-port
    :func:`waveguide_dielectric_slab_sparams`).

    A short has load reflection ``Gamma_L = -1``; an air section of length d transforms it to the
    port reference plane by the round-trip phase:

        S11 = -exp(-2 j beta_air d),   beta_air = sqrt(k0^2 - (pi/a)^2),   k0 = 2 pi f/c.

    Lossless, so ``|S11| = 1`` exactly; the OBSERVABLE is the reflection PHASE
    ``arg(S11) = pi - 2 beta_air d``, whose frequency slope is the round-trip GROUP DELAY
    ``tau = 2 d / v_group``. Returns ``{frequency, fc, beta_air, offset_length, S11, S11_mag,
    phase_rad, group_delay}``. The analytic gate a full-wave 1-port (port + PEC end) solve is checked
    against. Raises ValueError at/below the TE10 cutoff."""
    f = float(frequency); a = float(width_a); d = float(offset_length)
    if a <= 0.0:
        raise ValueError("width_a must be > 0 (got %r)" % (width_a,))
    if d < 0.0:
        raise ValueError("offset_length must be >= 0 (got %r)" % (offset_length,))
    fc = 0.5 * c / a
    if f <= fc:
        raise ValueError("frequency must exceed the TE10 cutoff fc=c/2a")
    k0 = 2.0 * math.pi * f / c
    kc = math.pi / a
    beta_air = math.sqrt(k0 * k0 - kc * kc)
    S11 = -cmath.exp(-2j * beta_air * d)
    vg = c * math.sqrt(1.0 - (fc / f) ** 2)                    # TE10 group velocity
    return {"frequency": f, "fc": fc, "beta_air": beta_air, "offset_length": d,
            "S11": S11, "S11_mag": abs(S11), "phase_rad": cmath.phase(S11),
            "group_delay": 2.0 * d / vg}


def waveguide_offset_short_length_from_group_delay(frequency, width_a, group_delay, c=C0):
    r"""Recover offset-short length [m] from measured/simulated reflection group delay.

    For a shorted air-filled rectangular TE10 guide section,

        tau_g = 2 d / v_group,     v_group = c sqrt(1 - (fc/f)^2),     fc = c/(2a),

    so the physical offset length is ``d = tau_g v_group / 2``.  This is the
    inverse post-processing companion to :func:`waveguide_offset_short_s11` and
    is useful when a port trace is available but the reference-plane offset is
    what should be calibrated. Raises ValueError at/below TE10 cutoff.
    """
    f = float(frequency); a = float(width_a); tau = float(group_delay)
    if a <= 0.0:
        raise ValueError("width_a must be > 0 (got %r)" % (width_a,))
    if tau < 0.0:
        raise ValueError("group_delay must be >= 0 (got %r)" % (group_delay,))
    fc = 0.5 * c / a
    disp = waveguide_dispersion(f, fc, c)
    length = 0.5 * tau * disp["v_group"]
    return {"frequency": f, "fc": fc, "group_delay": tau,
            "v_group": disp["v_group"], "offset_length": length}


def reflection_metrics(s11):
    r"""Common scalar readouts of a one-port reflection coefficient ``S11``.

    CST/VNA-style post-processing usually shows the same complex reflection in several units:

    - ``gamma = |S11|``
    - return loss ``RL = -20 log10(gamma)`` [dB]
    - reflected power fraction ``gamma^2``
    - delivered/matched power fraction ``1 - gamma^2``
    - mismatch loss ``ML = -10 log10(1 - gamma^2)`` [dB]
    - voltage standing-wave ratio ``VSWR = (1+gamma)/(1-gamma)``

    ``gamma=0`` gives infinite return loss and ``VSWR=1``.  ``gamma=1`` gives zero return loss,
    infinite VSWR, and infinite mismatch loss.
    """
    gamma = abs(complex(s11))
    if gamma > 1.0 + 1e-12:
        raise ValueError("|S11| must be <= 1 for a passive one-port")
    gamma = min(gamma, 1.0)
    reflected = gamma * gamma
    delivered = max(0.0, 1.0 - reflected)
    return_loss = math.inf if gamma == 0.0 else (-20.0 * math.log10(gamma))
    mismatch_loss = math.inf if delivered == 0.0 else (-10.0 * math.log10(delivered))
    return {
        "gamma": gamma,
        "return_loss_db": 0.0 if gamma == 1.0 else return_loss,
        "reflected_power_fraction": reflected,
        "delivered_power_fraction": delivered,
        "mismatch_loss_db": 0.0 if delivered == 1.0 else mismatch_loss,
        "vswr": math.inf if gamma == 1.0 else (1.0 + gamma) / (1.0 - gamma),
    }


def sparameter_group_delay(frequencies, s_values):
    r"""Group delay [s] from a sampled complex S-parameter trace.

    The RF/VNA definition is

        tau_g = - d arg(S) / d omega,

    with the phase unwrapped before differentiating. Returns one group-delay
    value per sample, using one-sided differences at the ends and central
    differences in the interior. This is the post-processing companion to
    :func:`waveguide_offset_short_s11`, dielectric-slab S-parameters, and measured
    or simulated RF traces.
    """
    freqs = [float(f) for f in frequencies]
    vals = [complex(s) for s in s_values]
    if len(freqs) != len(vals):
        raise ValueError("frequencies and s_values must have the same length")
    if len(freqs) < 2:
        raise ValueError("at least two samples are required")
    if any(b <= a for a, b in zip(freqs, freqs[1:])):
        raise ValueError("frequencies must be strictly increasing")

    phases = [cmath.phase(v) for v in vals]
    unwrapped = [phases[0]]
    offset = 0.0
    for phase in phases[1:]:
        value = phase + offset
        while value - unwrapped[-1] > math.pi:
            offset -= 2.0 * math.pi
            value = phase + offset
        while value - unwrapped[-1] <= -math.pi:
            offset += 2.0 * math.pi
            value = phase + offset
        unwrapped.append(value)

    omegas = [2.0 * math.pi * f for f in freqs]
    delays = []
    for i in range(len(freqs)):
        if i == 0:
            slope = (unwrapped[1] - unwrapped[0]) / (omegas[1] - omegas[0])
        elif i == len(freqs) - 1:
            slope = (unwrapped[-1] - unwrapped[-2]) / (omegas[-1] - omegas[-2])
        else:
            slope = (unwrapped[i + 1] - unwrapped[i - 1]) / (omegas[i + 1] - omegas[i - 1])
        delays.append(-slope)
    return delays


def helmholtz_cutoff_wavenumbers_2d(mesh, n_modes, bc="neumann", wall="wall",
                                    order=3, shift=-1.0):
    """Lowest cutoff wavenumbers k_c [1/m] of a waveguide cross-section by solving the 2D
    Helmholtz eigenproblem  -nabla_t^2 psi = k_c^2 psi  on ``mesh``.

    Generalised-eigenvalue solve (stiffness ``grad.grad`` vs mass) via NGSolve's
    :func:`ArnoldiSolver` with a small NEGATIVE shift, so the shifted operator (A - shift*M =
    A + |shift|*M) is SPD and never singular -- the eigenvalues nearest the shift are then the
    algebraically smallest. Works for ANY cross-section (rectangular, ridged, circular,
    L-shaped); only the geometry/mesh changes.

    Args:
        mesh    : 2D NGSolve mesh of the guide cross-section.
        n_modes : number of (physical) cutoff wavenumbers to return, ascending.
        bc      : 'neumann' for TE modes (dH_z/dn = 0), 'dirichlet' for TM modes (E_z = 0).
        wall    : boundary name of the metal wall (used as Dirichlet for the TM case).
        order   : H1 polynomial order (default 3 -- cutoffs converge fast).
        shift    : ArnoldiSolver shift (default -1.0; keep negative/below the spectrum).

    Returns the ascending list of k_c [1/m]; the trivial Neumann constant mode (k_c ~ 0) is
    dropped. Convert to frequency with :func:`cutoff_frequency`.
    """
    fes = H1(mesh, order=order, dirichlet=(wall if bc == "dirichlet" else ""))
    u, v = fes.TnT()
    A = BilinearForm(grad(u) * grad(v) * dx).Assemble()
    M = BilinearForm(u * v * dx).Assemble()
    nev = n_modes + (4 if bc == "neumann" else 2)          # pad for the dropped ~0 mode + spares
    gf = GridFunction(fes, multidim=nev)
    lams = ArnoldiSolver(A.mat, M.mat, fes.FreeDofs(), list(gf.vecs), shift=shift)
    kc2 = sorted(l.real for l in lams)
    kc2 = [x for x in kc2 if x > 1e-6 * max(kc2)]          # drop the constant (Neumann) null mode
    return [math.sqrt(x) for x in kc2[:n_modes]]


def laplace_dirichlet_eigenvalues(mesh, n_modes, order=2, shift=-1.0, dirichlet=".*"):
    """Lowest Laplace-Dirichlet eigenvalues of  -nabla^2 u = lambda u  with u = 0 on the
    Dirichlet boundary, by the generalised eigenproblem (grad.grad vs mass) via ArnoldiSolver
    with a small negative shift (so A + |shift| M is SPD).

    The 3-D SCALAR companion of :func:`helmholtz_cutoff_wavenumbers_2d` (2-D cross-section) and
    :func:`maxwell_cavity_modes_3d` (3-D vector). These eigenvalues are ALSO the modal decay
    rates of the transient heat equation (T_n ~ exp(-alpha lambda_n t), see
    ``multiphysics.solve_heat_transient``) and the squared acoustic/quantum eigen-wavenumbers.
    Closed forms: a BALL of radius R has lowest lambda = (pi/R)^2 (radial s-mode sin(pi r/R)/r),
    next (4.493409/R)^2 (l=1, x3); a BOX a x b x c has lambda = pi^2(p^2/a^2+q^2/b^2+r^2/c^2).

    Works on ANY 3-D mesh -- a Netgen tet mesh (call mesh.Curve for curved boundaries) OR a Cubit
    high-order curved-HEX ``.vol`` loaded as-is (do NOT mesh.Curve a loaded high-order .vol).
    Returns the ascending list of eigenvalues (length ``n_modes``)."""
    fes = H1(mesh, order=order, dirichlet=dirichlet)
    u, v = fes.TnT()
    A = BilinearForm(grad(u) * grad(v) * dx).Assemble()
    M = BilinearForm(u * v * dx).Assemble()
    gf = GridFunction(fes, multidim=n_modes + 4)
    lams = ArnoldiSolver(A.mat, M.mat, fes.FreeDofs(), list(gf.vecs), shift=shift)
    vals = sorted(l.real for l in lams)
    vals = [x for x in vals if x > 1e-6 * max(vals)]
    return vals[:n_modes]


def rectangular_cavity_frequency(a, b, d, m, n, p, c=C0):
    """Exact resonant FREQUENCY [Hz] of the TE_mnp / TM_mnp mode of a rectangular metallic CAVITY
    (a closed PEC box a x b x d):

        f_mnp = (c/2) sqrt((m/a)^2 + (n/b)^2 + (p/d)^2).

    The 3-D / closed-box counterpart of the (open-ended) waveguide cutoff
    :func:`rectangular_waveguide_cutoff` -- a third index p appears for the standing wave along the
    box length. The dominant mode (for a > d > b) is TE101. Used for microwave cavity filters,
    resonators, and accelerator cavities."""
    return 0.5 * c * math.sqrt((m / a) ** 2 + (n / b) ** 2 + (p / d) ** 2)


def cylindrical_cavity_frequency(radius, length, mode, m, n, p):
    """Exact resonant FREQUENCY [Hz] of a cylindrical metallic CAVITY.

    For a closed PEC cylinder of radius ``a`` and length ``L``:

        TM_mnp: k_r = j_mn/a       (zeros of J_m)
        TE_mnp: k_r = j'_mn/a      (zeros of J'_m)
        f = c/(2*pi) * sqrt(k_r^2 + (p*pi/L)^2)

    ``m`` is the azimuthal index, ``n`` the radial root index, and ``p`` the axial
    half-wave index.  The familiar pillbox accelerator mode TM010 is ``mode='TM',
    m=0, n=1, p=0``: it is independent of cavity length and equals the circular
    waveguide TM01 cutoff.  More generally, ``p=0`` cylindrical-cavity frequencies
    reduce to the corresponding circular-waveguide cutoffs; ``p>0`` adds the
    standing-wave term along the closed cavity axis.
    """
    if radius <= 0.0 or length <= 0.0:
        raise ValueError("radius and length must be positive")
    if m < 0 or n < 1 or p < 0:
        raise ValueError("indices require m>=0, n>=1, p>=0")
    from scipy.special import jn_zeros, jnp_zeros
    mm, nn, pp = int(m), int(n), int(p)
    md = str(mode).upper()
    if md == "TM":
        kr = jn_zeros(mm, nn)[nn - 1] / radius
    elif md == "TE":
        kr = jnp_zeros(mm, nn)[nn - 1] / radius
    else:
        raise ValueError("mode must be 'TE' or 'TM'")
    kz = pp * math.pi / length
    return C0 * math.hypot(kr, kz) / (2.0 * math.pi)


def maxwell_cavity_modes_3d(mesh, n_modes, shift, order=2, pec="pec"):
    """Lowest resonant frequencies [Hz] of a 3-D PEC cavity by the FULL-WAVE MAXWELL eigenproblem

        curl curl E = (omega/c)^2 E ,   n x E = 0 on the PEC walls,

    on an ``HCurl`` (edge / Nedelec) space -- the genuine vector-electromagnetic cavity (not the
    scalar Helmholtz box). The curl-curl operator has a LARGE gradient kernel (curl grad = 0, the
    spurious zero modes); they are suppressed by giving ``ArnoldiSolver`` a ``shift`` near the
    target k^2 (the shift-invert then amplifies the PHYSICAL resonances and starves the far-away
    null space). Estimate ``shift`` from the expected fundamental, e.g.
    ``(pi/a)**2 + (pi/d)**2`` for TE101.

    Args:
        mesh    : 3-D NGSolve mesh of the cavity; PEC walls named ``pec``.
        n_modes : number of resonant frequencies to return, ascending.
        shift    : ArnoldiSolver shift ~ the squared wavenumber k^2 of interest [1/m^2].
        order   : HCurl order (default 2).
        pec     : boundary name of the metal walls (tangential E = 0).

    Returns the ascending list of resonant frequencies [Hz] (the residual near-zero gradient modes
    are dropped). Validated against the exact box spectrum."""
    fes = HCurl(mesh, order=order, dirichlet=pec)
    u, v = fes.TnT()
    A = BilinearForm(curl(u) * curl(v) * dx).Assemble()
    M = BilinearForm(u * v * dx).Assemble()
    gf = GridFunction(fes, multidim=n_modes + 4)
    lams = ArnoldiSolver(A.mat, M.mat, fes.FreeDofs(), list(gf.vecs), shift=shift)
    k2 = sorted(l.real for l in lams if l.real > 1e-3 * shift)     # drop the gradient null space
    return [C0 * math.sqrt(x) / (2.0 * math.pi) for x in k2[:n_modes]]


def skin_depth(frequency, sigma, mu_r=1.0):
    """Skin depth delta = sqrt(2/(omega mu sigma)) [m] -- the 1/e penetration of an AC field into a
    conductor (mu = mu0 mu_r). Sets the surface resistance of cavity / waveguide walls and lines; at
    10 GHz copper delta ~ 0.66 um."""
    omega = 2.0 * math.pi * frequency
    return math.sqrt(2.0 / (omega * MU0 * mu_r * sigma))


def surface_resistance(frequency, sigma, mu_r=1.0):
    """Surface resistance R_s = sqrt(omega mu/(2 sigma)) = 1/(sigma*delta) [Ohm] of a good conductor --
    the real part of the surface impedance. The time-average wall loss per unit area is
    (R_s/2)|H_tan|^2, so R_s sets the cavity/line loss; it grows as sqrt(frequency)."""
    omega = 2.0 * math.pi * frequency
    return math.sqrt(omega * MU0 * mu_r / (2.0 * sigma))


def rectangular_cavity_q(a, b, d, sigma, mu_r=1.0):
    """Wall-loss (unloaded) quality factor Q of the dominant TE101 mode of a rectangular cavity
    (a x b x d; a along x, d along z, b the height) with conducting walls -- the LOSS sequel to the
    lossless resonance :func:`rectangular_cavity_frequency` / :func:`maxwell_cavity_modes_3d` (#59).
    Q = omega * U / P_wall, with U the stored energy and P_wall the wall dissipation; for TE101

        Q = (k a d)^3 b eta / (2 pi^2 R_s (2 a^3 b + 2 b d^3 + a^3 d + a d^3)),

    k = omega/c, eta = mu0 c the free-space impedance, R_s = :func:`surface_resistance`. Q ~
    (volume/surface)/delta ~ sqrt(sigma): a better conductor (smaller surface resistance / skin depth)
    gives a sharper resonance -- the bridge from the wave family to the skin depth (#45/#55).
    Validated against the numerically-integrated TE101 field energy and wall loss."""
    f = rectangular_cavity_frequency(a, b, d, 1, 0, 1)
    k = 2.0 * math.pi * f / C0
    Rs = surface_resistance(f, sigma, mu_r)
    eta = MU0 * C0
    num = (k * a * d) ** 3 * b * eta
    den = 2.0 * math.pi ** 2 * Rs * (2 * a ** 3 * b + 2 * b * d ** 3 + a ** 3 * d + a * d ** 3)
    return num / den
