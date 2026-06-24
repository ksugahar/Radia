r"""Electrostatics closed forms for radia-ngsolve -- capacitance, conduction
resistance, polarizability, and field-stress textbook results.

A small library of *analytic* electrostatic / steady-conduction formulae, the
electric-field counterparts of the magnetostatic helpers elsewhere in the
package.  Each routine is a pure closed form (no field solve), so it doubles as
a fast design tool and as the exact reference a numerical solver is checked
against.  The same Laplace operator underlies all of them; only the boundary
data and the constitutive constant (permittivity vs. conductivity) change,
which is why the resistance / capacitance duals share an algebraic skeleton.

References:
  * J. D. Jackson, "Classical Electrodynamics", 3rd ed., Secs. 2.5, 4.4.
  * D. J. Griffiths, "Introduction to Electrodynamics", 4th ed., Ex. 4.2, Ch. 2.
  * W. R. Smythe, "Static and Dynamic Electricity", 3rd ed.
  * R. Holm, "Electric Contacts: Theory and Application", 4th ed.
  * H. A. Pohl, "Dielectrophoresis", Cambridge Univ. Press, 1978.
  * J. C. Maxwell / K. W. Wagner interfacial-polarization layered dielectric.
"""
import math

EPS0 = 8.8541878128e-12          # vacuum permittivity [F/m] (CODATA)


def isolated_disk_capacitance(radius_a, eps_r=1.0):
    r"""Self-capacitance of an isolated thin conducting DISK of radius ``a`` [m]
    in a medium of relative permittivity ``eps_r``:

        C = 8 eps0 eps_r a .

    This is the flat-disk analogue of the isolated-sphere result C = 4 pi eps0 a;
    the disk holds 2/pi as much charge per volt as a sphere of the same radius.
    Returns ``{"C": ...}`` [F].  (Jackson; Smythe.)
    """
    if radius_a <= 0:
        raise ValueError("radius_a must be > 0")
    if eps_r <= 0:
        raise ValueError("eps_r must be > 0")
    return {"C": 8.0 * EPS0 * eps_r * radius_a}


def spreading_resistance_disk(resistivity_rho, radius_a):
    r"""Constriction (spreading) resistance of current entering a half-space
    through a circular contact spot of radius ``a`` [m], resistivity ``rho``
    [ohm.m].  The Holm constriction resistance for a flat circular a-spot is

        R = rho / (4 a) ,

    while current spreading from a HEMISPHERICAL contact of the same radius gives

        R_hemi = rho / (2 pi a) .

    Their ratio R / R_hemi = pi/2 is geometry-only (rho-independent).  Returns
    ``{"R": ..., "hemispherical_R": ...}`` [ohm].  (Holm, "Electric Contacts".)
    """
    if resistivity_rho <= 0:
        raise ValueError("resistivity_rho must be > 0")
    if radius_a <= 0:
        raise ValueError("radius_a must be > 0")
    return {
        "R": resistivity_rho / (4.0 * radius_a),
        "hemispherical_R": resistivity_rho / (2.0 * math.pi * radius_a),
    }


def coaxial_shell_resistance(resistivity_rho, r_inner, r_outer, length_L):
    r"""Radial DC resistance of a coaxial conducting shell (current flows
    radially through a cylindrical sleeve a<r<b, length ``L``):

        R = rho ln(b/a) / (2 pi L) .

    This is the RC dual of the coaxial CAPACITANCE C = 2 pi eps0 L / ln(b/a):
    the same ln(b/a) geometry factor governs both, so R * C = rho eps0 for a
    coax filled with one medium.  Returns ``{"R": ...}`` [ohm].  (Smythe.)
    """
    if resistivity_rho <= 0:
        raise ValueError("resistivity_rho must be > 0")
    if length_L <= 0:
        raise ValueError("length_L must be > 0")
    if not (r_outer > r_inner > 0):
        raise ValueError("require r_outer > r_inner > 0")
    return {"R": resistivity_rho * math.log(r_outer / r_inner) / (2.0 * math.pi * length_L)}


def layered_parallel_plate_capacitance(area, thicknesses, eps_r_layers):
    r"""Capacitance of a parallel-plate stack of dielectric layers normal to the
    plates (layers in SERIES on the field line):

        C = eps0 area / sum_i (d_i / eps_r_i) ,

    with effective permittivity (the series-average that reproduces a single slab
    of the total thickness)

        eps_eff = sum_i d_i / sum_i (d_i / eps_r_i) .

    Each layer carries the same normal D, so its slabs add reciprocally -- the
    Maxwell-Wagner interfacial-polarization stack.  ``thicknesses`` and
    ``eps_r_layers`` are equal-length sequences.  Returns
    ``{"C": ..., "eps_eff": ...}``.  (Maxwell-Wagner.)
    """
    if area <= 0:
        raise ValueError("area must be > 0")
    if len(thicknesses) != len(eps_r_layers):
        raise ValueError("thicknesses and eps_r_layers must have equal length")
    if len(thicknesses) == 0:
        raise ValueError("need at least one layer")
    if any(d <= 0 for d in thicknesses):
        raise ValueError("all thicknesses must be > 0")
    if any(e <= 0 for e in eps_r_layers):
        raise ValueError("all eps_r_layers must be > 0")
    series = sum(d / e for d, e in zip(thicknesses, eps_r_layers))
    total = sum(thicknesses)
    return {"C": EPS0 * area / series, "eps_eff": total / series}


def parallel_plate_capacitor_energy_force(eps_r, area, gap, voltage):
    r"""Capacitance, stored energy, and plate-attraction force of an ideal
    parallel-plate capacitor (plate ``area``, ``gap`` d, applied ``voltage`` V):

        C = eps0 eps_r area / d ,
        W = 1/2 C V^2 ,
        |F| = 1/2 eps0 eps_r area V^2 / d^2 = W / d .

    At fixed voltage the attractive force equals the energy density times the
    plate area, i.e. |F| = W/d, and scales as 1/d^2.  Returns
    ``{"C": ..., "energy": ..., "force": ...}``.  (Griffiths.)
    """
    if eps_r <= 0:
        raise ValueError("eps_r must be > 0")
    if area <= 0:
        raise ValueError("area must be > 0")
    if gap <= 0:
        raise ValueError("gap must be > 0")
    C = EPS0 * eps_r * area / gap
    energy = 0.5 * C * voltage * voltage
    force = 0.5 * EPS0 * eps_r * area * voltage * voltage / (gap * gap)
    return {"C": C, "energy": energy, "force": force}


def dielectric_sphere_polarizability(radius_a, eps_r):
    r"""Polarizability of a uniform dielectric sphere of radius ``a`` and
    relative permittivity ``eps_r`` in a uniform external field:

        alpha = 4 pi eps0 a^3 (eps_r - 1) / (eps_r + 2) ,

    and the volume-average induced dipole per applied field

        m_avg / E0 = 3 (eps_r - 1) / (eps_r + 2) .

    Vanishes for eps_r = 1 (no contrast) and saturates to the conducting-sphere
    value alpha -> 4 pi eps0 a^3 as eps_r -> infinity.  Returns
    ``{"alpha": ..., "M_avg_per_E0": ...}``.  (Jackson, Sec. 4.4.)
    """
    if radius_a <= 0:
        raise ValueError("radius_a must be > 0")
    if eps_r <= 0:
        raise ValueError("eps_r must be > 0")
    cm = (eps_r - 1.0) / (eps_r + 2.0)
    return {
        "alpha": 4.0 * math.pi * EPS0 * radius_a ** 3 * cm,
        "M_avg_per_E0": 3.0 * cm,
    }


def dielectric_sphere_uniform_field(radius_a, eps_r, applied_field_E0=1.0,
                                    sample_radius=None):
    r"""Dielectric sphere in a uniform applied field along the z-axis.

    A sphere of radius ``a`` and relative permittivity ``eps_r`` embedded in
    vacuum and driven by an applied field ``E0 zhat`` has a uniform interior
    field

        E_in = 3 E0 / (eps_r + 2),

    and an exterior dipole correction set by the Clausius-Mossotti factor
    ``K = (eps_r - 1)/(eps_r + 2)``.  On the axis and equator at radius ``r>a``:

        E_z,axis = E0 (1 + 2 K (a/r)^3),
        E_z,eq   = E0 (1 -   K (a/r)^3).

    Returns the interior field, polarizability, dipole moment, and, when
    ``sample_radius`` is supplied, the axial/equatorial exterior field samples.
    """
    if radius_a <= 0:
        raise ValueError("radius_a must be > 0")
    if eps_r <= 0:
        raise ValueError("eps_r must be > 0")
    a = float(radius_a)
    er = float(eps_r)
    e0 = float(applied_field_E0)
    cm = (er - 1.0) / (er + 2.0)
    alpha = 4.0 * math.pi * EPS0 * a ** 3 * cm
    out = {
        "clausius_mossotti": cm,
        "interior_field": 3.0 * e0 / (er + 2.0),
        "interior_field_factor": 3.0 / (er + 2.0),
        "polarizability": alpha,
        "dipole_moment": alpha * e0,
    }
    if sample_radius is not None:
        r = float(sample_radius)
        if r <= a:
            raise ValueError("sample_radius must be > radius_a for exterior field samples")
        q = (a / r) ** 3
        out.update({
            "sample_radius": r,
            "axial_field": e0 * (1.0 + 2.0 * cm * q),
            "equatorial_field": e0 * (1.0 - cm * q),
        })
    return out


def uniformly_polarized_sphere_field(polarization_P, radius_a):
    r"""Fields of a uniformly polarized dielectric sphere (frozen polarization
    ``P`` [C/m^2], radius ``a``):

        E_in = -P / (3 eps0)         (uniform inside),
        D_in = eps0 E_in + P = (2/3) P ,
        p    = (4/3) pi a^3 P         (exterior is an ideal point dipole).

    The interior depolarizing field is the 1/3 of the bound-charge field that
    makes the depolarization factor of a sphere exactly 1/3.  Returns
    ``{"E_in": ..., "D_in": ..., "dipole_moment": ...}``.  (Griffiths, Ex. 4.2.)
    """
    if radius_a <= 0:
        raise ValueError("radius_a must be > 0")
    return {
        "E_in": -polarization_P / (3.0 * EPS0),
        "D_in": (2.0 / 3.0) * polarization_P,
        "dipole_moment": (4.0 / 3.0) * math.pi * radius_a ** 3 * polarization_P,
    }


def charged_conducting_sphere_surface_stress(charge_Q, radius_a, eps_r=1.0):
    r"""Surface field and electrostatic (Maxwell) tension on a charged conducting
    sphere (total charge ``Q``, radius ``a``, surrounding medium ``eps_r``):

        E_surf  = Q / (4 pi eps0 eps_r a^2) ,
        p_out   = eps0 eps_r E_surf^2 / 2 = sigma^2 / (2 eps0 eps_r) ,

    where sigma = Q/(4 pi a^2) is the surface charge density.  The pressure is
    OUTWARD (the conductor is always pulled into the field, regardless of charge
    sign) and equals the local electrostatic energy density.  Returns
    ``{"surface_field": ..., "pressure": ...}``.  (Jackson, Sec. 2.5.)
    """
    if radius_a <= 0:
        raise ValueError("radius_a must be > 0")
    if eps_r <= 0:
        raise ValueError("eps_r must be > 0")
    E = charge_Q / (4.0 * math.pi * EPS0 * eps_r * radius_a ** 2)
    return {
        "surface_field": E,
        "pressure": 0.5 * EPS0 * eps_r * E * E,
    }


def sphere_above_plane_capacitance(radius_a, height_h, eps_r=1.0, n_terms=20):
    r"""Capacitance of a conducting sphere of radius ``a`` whose center sits at
    height ``h`` above a grounded conducting plane, via the classical image
    series (Smythe):

        C = 4 pi eps0 eps_r a sum_{n>=1} sinh(alpha) / sinh(n alpha),
        alpha = acosh(h/a) .

    Requires h > a (sphere clears the plane).  As h/a -> infinity every term
    past n=1 dies and C -> 4 pi eps0 eps_r a, the isolated-sphere value.  The
    series is truncated at ``n_terms`` (~20 converges for h >~ a).  Returns
    ``{"C": ..., "n_terms": ...}``.  (Smythe; method of images.)
    """
    if radius_a <= 0:
        raise ValueError("radius_a must be > 0")
    if eps_r <= 0:
        raise ValueError("eps_r must be > 0")
    if height_h <= radius_a:
        raise ValueError("require height_h > radius_a (sphere must clear the plane)")
    if n_terms < 1:
        raise ValueError("n_terms must be >= 1")
    alpha = math.acosh(height_h / radius_a)
    s = sum(math.sinh(alpha) / math.sinh(n * alpha) for n in range(1, n_terms + 1))
    return {"C": 4.0 * math.pi * EPS0 * eps_r * radius_a * s, "n_terms": n_terms}


def dielectrophoresis_force(radius_a, eps_medium_r, eps_particle_r, grad_E_squared):
    r"""Time-averaged dielectrophoretic (DEP) force on a small dielectric sphere
    of radius ``a`` suspended in a medium of permittivity ``eps_medium_r`` in a
    non-uniform field (gradient of E^2):

        F = 2 pi eps0 eps_m a^3 Re{K} grad(E^2),
        K = (eps_p - eps_m) / (eps_p + 2 eps_m)      (Clausius-Mossotti factor).

    The Clausius-Mossotti factor changes sign at eps_p = eps_m: a particle more
    polarizable than its medium is pulled toward field maxima (positive DEP),
    a less polarizable one is repelled (negative DEP).  Returns
    ``{"force": ..., "clausius_mossotti": ...}``.  (Pohl, "Dielectrophoresis".)
    """
    if radius_a <= 0:
        raise ValueError("radius_a must be > 0")
    if eps_medium_r <= 0:
        raise ValueError("eps_medium_r must be > 0")
    if eps_particle_r <= 0:
        raise ValueError("eps_particle_r must be > 0")
    cm = (eps_particle_r - eps_medium_r) / (eps_particle_r + 2.0 * eps_medium_r)
    return {
        "force": 2.0 * math.pi * EPS0 * eps_medium_r * radius_a ** 3 * cm * grad_E_squared,
        "clausius_mossotti": cm,
    }
