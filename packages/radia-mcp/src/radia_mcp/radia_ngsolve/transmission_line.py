"""TEM transmission-line per-unit-length parameters -- closed-form analytic helpers.

A two-conductor line carrying a transverse-electromagnetic (TEM) wave is fully
described by its per-unit-length capacitance C [F/m] and inductance L [H/m]. For
ANY lossless TEM line filled with a homogeneous medium (eps_r, mu_r) the two are
linked by the geometry-independent identity

    L * C = mu0 * mu_r * eps0 * eps_r ,

so the line's characteristic impedance and phase velocity follow as

    Z0 = sqrt(L / C) ,    vp = 1 / sqrt(L * C) = c / sqrt(eps_r * mu_r) .

Only the cross-section geometry enters C and L individually (through a single
logarithm or inverse-hyperbolic-cosine factor); their product cancels it. This
module collects the textbook closed forms for the canonical cross-sections
(coaxial, two-wire, microstrip, wire-over-ground-plane).

References:
    * D. M. Pozar, "Microwave Engineering," 4th ed., Sec. 2.2 (coax) and
      Sec. 3.8 (microstrip).
    * S. Ramo, J. R. Whinnery, T. Van Duzer, "Fields and Waves in Communication
      Electronics," 3rd ed. (two-wire line, capacitance of conductor pairs).
    * E. O. Hammerstad and O. Jensen, "Accurate Models for Microstrip
      Computer-Aided Design" (Hammerstad-Wheeler microstrip synthesis).
    * W. R. Smythe, "Static and Dynamic Electricity" (capacitance by images).
"""
import math

EPS0 = 8.8541878128e-12       # vacuum permittivity [F/m]
MU0 = 4.0e-7 * math.pi        # vacuum permeability [H/m]
C0 = 299792458.0              # speed of light in vacuum [m/s]
ETA0 = 376.730313668          # impedance of free space [Ohm] = sqrt(MU0/EPS0)


def tem_lc_identity_summary(C_per_m, L_per_m, eps_r=1.0, mu_r=1.0):
    """Summarize the lossless TEM ``L*C`` identity for a homogeneous line.

    A quasi-static electrostatic solve gives ``C_per_m`` and an independent
    magnetostatic solve gives ``L_per_m``.  For any homogeneous lossless TEM
    two-conductor line their product must satisfy

        L C = mu0 mu_r eps0 eps_r,

    independent of the cross-section geometry.  This helper packages the
    characteristic impedance, phase velocity, and relative identity errors in a
    form suitable for CST/radia/analytic cross-validation tables.
    """
    if C_per_m <= 0.0 or L_per_m <= 0.0:
        raise ValueError("C_per_m and L_per_m must be positive")
    if eps_r <= 0.0 or mu_r <= 0.0:
        raise ValueError("eps_r and mu_r must be positive")

    lc = L_per_m * C_per_m
    expected_lc = MU0 * mu_r * EPS0 * eps_r
    vp = 1.0 / math.sqrt(lc)
    vp_expected = C0 / math.sqrt(eps_r * mu_r)
    z0 = math.sqrt(L_per_m / C_per_m)
    return {
        "C_per_m": C_per_m,
        "L_per_m": L_per_m,
        "Z0": z0,
        "vp": vp,
        "LC_product": lc,
        "LC_expected": expected_lc,
        "LC_relative_error": (lc - expected_lc) / expected_lc,
        "vp_expected": vp_expected,
        "vp_relative_error": (vp - vp_expected) / vp_expected,
    }


def coaxial_line_parameters(inner_radius_a, outer_radius_b, eps_r=1.0, mu_r=1.0):
    """Per-unit-length parameters of a coaxial line (inner radius a, outer radius b).

    Closed form (Pozar, "Microwave Engineering," Sec. 2.2):

        C = 2*pi*eps0*eps_r / ln(b/a)          [F/m]
        L = mu0*mu_r * ln(b/a) / (2*pi)        [H/m]
        Z0 = sqrt(L/C) = (eta0/(2*pi)) sqrt(mu_r/eps_r) ln(b/a)
        vp = 1/sqrt(L*C) = c / sqrt(eps_r*mu_r)

    The single geometry factor is ln(b/a); it cancels in L*C, so vp is
    geometry-independent. Requires b > a > 0 and positive material constants.

    Returns dict: C_per_m, L_per_m, Z0, vp, ln_ratio.
    """
    if inner_radius_a <= 0.0 or outer_radius_b <= 0.0:
        raise ValueError("coaxial radii must be positive")
    if outer_radius_b <= inner_radius_a:
        raise ValueError("outer_radius_b must exceed inner_radius_a")
    if eps_r <= 0.0 or mu_r <= 0.0:
        raise ValueError("eps_r and mu_r must be positive")

    ln_ratio = math.log(outer_radius_b / inner_radius_a)
    C = 2.0 * math.pi * EPS0 * eps_r / ln_ratio
    L = MU0 * mu_r * ln_ratio / (2.0 * math.pi)
    Z0 = math.sqrt(L / C)
    vp = 1.0 / math.sqrt(L * C)
    return {"C_per_m": C, "L_per_m": L, "Z0": Z0, "vp": vp, "ln_ratio": ln_ratio}


def two_wire_line_parameters(wire_radius_a, separation_D, eps_r=1.0, mu_r=1.0):
    """Per-unit-length parameters of a two-wire (parallel-wire) line.

    Two identical conductors of radius a with center-to-center separation D
    (Ramo, Whinnery, Van Duzer, "Fields and Waves in Communication
    Electronics"):

        C = pi*eps0*eps_r / acosh(D/(2a))      [F/m]
        L = mu0*mu_r * acosh(D/(2a)) / pi      [H/m]
        Z0 = sqrt(L/C) = (eta0/pi) sqrt(mu_r/eps_r) acosh(D/(2a))
        vp = 1/sqrt(L*C) = c / sqrt(eps_r*mu_r)

    The geometry factor acosh(D/(2a)) cancels in L*C. Requires D > 2a (wires
    must not overlap) and positive material constants.

    Returns dict: C_per_m, L_per_m, Z0, vp, acosh_arg.
    """
    if wire_radius_a <= 0.0 or separation_D <= 0.0:
        raise ValueError("wire radius and separation must be positive")
    if separation_D <= 2.0 * wire_radius_a:
        raise ValueError("separation_D must exceed 2*wire_radius_a (wires overlap)")
    if eps_r <= 0.0 or mu_r <= 0.0:
        raise ValueError("eps_r and mu_r must be positive")

    acosh_arg = separation_D / (2.0 * wire_radius_a)
    factor = math.acosh(acosh_arg)
    C = math.pi * EPS0 * eps_r / factor
    L = MU0 * mu_r * factor / math.pi
    Z0 = math.sqrt(L / C)
    vp = 1.0 / math.sqrt(L * C)
    return {"C_per_m": C, "L_per_m": L, "Z0": Z0, "vp": vp, "acosh_arg": acosh_arg}


def microstrip_line_parameters(width_w, substrate_h, eps_r):
    """Microstrip effective permittivity and characteristic impedance.

    Hammerstad-Wheeler closed-form analysis (Pozar, "Microwave Engineering,"
    Sec. 3.8; Hammerstad & Jensen). With normalized width u = w/h:

        eps_eff = (eps_r+1)/2 + (eps_r-1)/2 * (1 + 12/u)**-0.5
                  (+ 0.04*(1-u)**2 correction term when u < 1)

        Z0 = (60/sqrt(eps_eff)) * ln(8/u + u/4)                       (u <= 1)
        Z0 = (120*pi/sqrt(eps_eff)) / (u + 1.393 + 0.667*ln(u+1.444)) (u > 1)

        vp = c / sqrt(eps_eff)

    Because the fields straddle the dielectric substrate and the air above, the
    effective permittivity lies between 1 and eps_r. Requires eps_r >= 1, w > 0,
    h > 0.

    Returns dict: eps_eff, Z0, vp, u.
    """
    if eps_r < 1.0:
        raise ValueError("eps_r must be >= 1")
    if width_w <= 0.0 or substrate_h <= 0.0:
        raise ValueError("width_w and substrate_h must be positive")

    u = width_w / substrate_h
    eps_eff = (eps_r + 1.0) / 2.0 + (eps_r - 1.0) / 2.0 * (1.0 + 12.0 / u) ** -0.5
    if u < 1.0:
        eps_eff += (eps_r - 1.0) / 2.0 * 0.04 * (1.0 - u) ** 2

    if u <= 1.0:
        Z0 = (60.0 / math.sqrt(eps_eff)) * math.log(8.0 / u + u / 4.0)
    else:
        Z0 = (120.0 * math.pi / math.sqrt(eps_eff)) / (
            u + 1.393 + 0.667 * math.log(u + 1.444)
        )

    vp = C0 / math.sqrt(eps_eff)
    return {"eps_eff": eps_eff, "Z0": Z0, "vp": vp, "u": u}


def wire_capacitance_per_length(geometry, radius, spacing, eps_r=1.0):
    """Per-unit-length capacitance C [F/m] for a single canonical geometry.

    Closed forms by the method of images / conformal mapping (Smythe, "Static
    and Dynamic Electricity"; Ramo, Whinnery, Van Duzer):

        'wire_plane' : wire of radius r at height h above a ground plane
                       C = 2*pi*eps / acosh(h/r)      (spacing = h)
        'coax'       : coax inner radius a, outer radius b
                       C = 2*pi*eps / ln(b/a)         (radius = a, spacing = b)
        'two_wire'   : two wires of radius r separated by D
                       C = pi*eps / acosh(D/(2r))     (spacing = D)

    where eps = eps0*eps_r. By image symmetry the wire-over-plane at height h
    equals one half of a two-wire line of separation D = 2h carrying twice the
    capacitance: C_plane(h) = 2 * C_two_wire(D = 2h).

    Returns dict: C_per_m, geometry.
    """
    if radius <= 0.0 or spacing <= 0.0:
        raise ValueError("radius and spacing must be positive")
    if eps_r <= 0.0:
        raise ValueError("eps_r must be positive")
    eps = EPS0 * eps_r

    if geometry == "wire_plane":
        if spacing <= radius:
            raise ValueError("height (spacing) must exceed wire radius")
        C = 2.0 * math.pi * eps / math.acosh(spacing / radius)
    elif geometry == "coax":
        if spacing <= radius:
            raise ValueError("outer radius (spacing) must exceed inner radius")
        C = 2.0 * math.pi * eps / math.log(spacing / radius)
    elif geometry == "two_wire":
        if spacing <= 2.0 * radius:
            raise ValueError("separation (spacing) must exceed 2*radius")
        C = math.pi * eps / math.acosh(spacing / (2.0 * radius))
    else:
        raise ValueError(
            "geometry must be one of 'wire_plane', 'coax', 'two_wire'"
        )

    return {"C_per_m": C, "geometry": geometry}
