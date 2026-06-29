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


def coaxial_rc_duality_summary(conductivity_sigma, eps_r, r_inner, r_outer, length_L=1.0,
                               measured_resistance_ohm=None, measured_capacitance_F=None):
    r"""Summary gate for the EC/ES duality of a coaxial annulus.

    The same annular geometry solves two Laplace problems:

    * steady radial conduction: ``R = ln(b/a)/(2*pi*sigma*L)``
    * electrostatics: ``C = 2*pi*eps0*eps_r*L/ln(b/a)``

    Their product cancels the geometry and leaves the dielectric relaxation time
    ``R*C = eps/sigma``.  This is a compact teaching gate for checking that an
    EC solver and an ES solver used the same full circular boundaries and units.
    """

    if conductivity_sigma <= 0:
        raise ValueError("conductivity_sigma must be > 0")
    if eps_r <= 0:
        raise ValueError("eps_r must be > 0")
    if length_L <= 0:
        raise ValueError("length_L must be > 0")
    if not (r_outer > r_inner > 0):
        raise ValueError("require r_outer > r_inner > 0")

    sigma = float(conductivity_sigma)
    er = float(eps_r)
    length = float(length_L)
    resistance = coaxial_shell_resistance(1.0 / sigma, r_inner, r_outer, length)["R"]
    capacitance = 2.0 * math.pi * EPS0 * er * length / math.log(r_outer / r_inner)
    tau = resistance * capacitance
    expected_tau = EPS0 * er / sigma
    out = {
        "policy": "coaxial_annulus_ec_es_rc_duality",
        "R_ohm": resistance,
        "C_F": capacitance,
        "tau_s": tau,
        "expected_tau_s": expected_tau,
        "tau_rel_error": abs(tau - expected_tau) / expected_tau,
        "conductivity_sigma": sigma,
        "eps_r": er,
        "r_inner": float(r_inner),
        "r_outer": float(r_outer),
        "length_L": length,
        "checks": {
            "positive_R": resistance > 0.0,
            "positive_C": capacitance > 0.0,
            "geometry_cancels_in_RC": abs(tau - expected_tau) / expected_tau < 1.0e-12,
        },
    }
    if measured_resistance_ohm is not None:
        measured_r = float(measured_resistance_ohm)
        out["measured_R_ohm"] = measured_r
        out["measured_R_rel_error"] = abs(measured_r - resistance) / resistance
    if measured_capacitance_F is not None:
        measured_c = float(measured_capacitance_F)
        out["measured_C_F"] = measured_c
        out["measured_C_rel_error"] = abs(measured_c - capacitance) / capacitance
    if measured_resistance_ohm is not None and measured_capacitance_F is not None:
        measured_tau = float(measured_resistance_ohm) * float(measured_capacitance_F)
        out["measured_tau_s"] = measured_tau
        out["measured_tau_rel_error"] = abs(measured_tau - expected_tau) / expected_tau
    return out


def coaxial_capacitor_energy_force(eps_r, r_inner, r_outer, length_L, voltage):
    r"""Capacitance, energy, field pressure, and radial force of a coaxial capacitor.

    For a coaxial capacitor with inner radius ``a``, outer radius ``b``,
    length ``L``, dielectric ``eps = eps0 eps_r``, and voltage ``V``:

        C = 2 pi eps L / ln(b/a),
        E(r) = V / (r ln(b/a)),
        p(r) = 1/2 eps E(r)^2.

    The fixed-voltage generalized force for increasing the inner radius is
    positive and equals the Maxwell pressure integrated over the inner
    conductor surface:

        F_a = 1/2 V^2 dC/da = p(a) 2 pi a L.

    Increasing the outer radius lowers capacitance, so ``F_b`` is negative in
    the ``+b`` coordinate and ``|F_b| = p(b) 2 pi b L``.  A complete coax has
    zero net vector force by angular symmetry; these are radial coordinate
    forces / surface pressures for validation and actuator-style sweeps.
    """

    if eps_r <= 0:
        raise ValueError("eps_r must be > 0")
    if length_L <= 0:
        raise ValueError("length_L must be > 0")
    if not (r_outer > r_inner > 0):
        raise ValueError("require r_outer > r_inner > 0")
    eps = EPS0 * float(eps_r)
    a = float(r_inner)
    b = float(r_outer)
    length = float(length_L)
    v = float(voltage)
    log_ratio = math.log(b / a)
    capacitance = 2.0 * math.pi * eps * length / log_ratio
    energy = 0.5 * capacitance * v * v
    e_inner = v / (a * log_ratio)
    e_outer = v / (b * log_ratio)
    p_inner = 0.5 * eps * e_inner * e_inner
    p_outer = 0.5 * eps * e_outer * e_outer
    dC_da = 2.0 * math.pi * eps * length / (a * log_ratio * log_ratio)
    dC_db = -2.0 * math.pi * eps * length / (b * log_ratio * log_ratio)
    force_inner = 0.5 * v * v * dC_da
    force_outer = 0.5 * v * v * dC_db
    return {
        "C": capacitance,
        "energy": energy,
        "eps_r": float(eps_r),
        "r_inner": a,
        "r_outer": b,
        "length_L": length,
        "voltage": v,
        "log_radius_ratio": log_ratio,
        "electric_field_inner_V_per_m": e_inner,
        "electric_field_outer_V_per_m": e_outer,
        "pressure_inner_Pa": p_inner,
        "pressure_outer_Pa": p_outer,
        "dCdr_inner_F_per_m": dC_da,
        "dCdr_outer_F_per_m": dC_db,
        "inner_radius_force_N": force_inner,
        "outer_radius_force_N": force_outer,
        "inner_pressure_area_force_N": p_inner * 2.0 * math.pi * a * length,
        "outer_pressure_area_force_N": -p_outer * 2.0 * math.pi * b * length,
    }


def spherical_capacitor_energy_summary(eps_r, r_inner, r_outer, voltage, sample_radius=None,
                                       measured_capacitance_F=None,
                                       measured_energy_J=None,
                                       measured_sample_field_V_per_m=None):
    r"""Closed-form spherical-capacitor gate for 3D electrostatic solvers.

    A spherical shell with inner radius ``a``, outer radius ``b``, dielectric
    ``eps = eps0 eps_r``, inner voltage ``V`` and grounded outer conductor has

        C = 4 pi eps a b / (b-a),
        E(r) = V a b / ((b-a) r^2),
        W = 1/2 C V^2,
        p(r) = 1/2 eps E(r)^2.

    The pressure ratio is especially useful as a geometry-only check:

        p_inner / p_outer = (b/a)^4.

    Optional measured values are compared against the exact reference so a
    notebook, NGSolve validation, or commercial cross-check can write one compact
    JSON row without embedding solver-specific provenance in this public helper.
    """

    if eps_r <= 0:
        raise ValueError("eps_r must be > 0")
    if not (r_outer > r_inner > 0):
        raise ValueError("require r_outer > r_inner > 0")
    if sample_radius is not None and not (r_inner < sample_radius < r_outer):
        raise ValueError("sample_radius must lie between r_inner and r_outer")

    eps = EPS0 * float(eps_r)
    a = float(r_inner)
    b = float(r_outer)
    v = float(voltage)
    radius = math.sqrt(a * b) if sample_radius is None else float(sample_radius)
    capacitance = 4.0 * math.pi * eps * a * b / (b - a)
    energy = 0.5 * capacitance * v * v
    e_inner = v * a * b / ((b - a) * a * a)
    e_outer = v * a * b / ((b - a) * b * b)
    e_sample = v * a * b / ((b - a) * radius * radius)
    p_inner = 0.5 * eps * e_inner * e_inner
    p_outer = 0.5 * eps * e_outer * e_outer
    out = {
        "policy": "spherical_capacitor_energy_field_gate",
        "C_F": capacitance,
        "energy_J": energy,
        "eps_r": float(eps_r),
        "r_inner": a,
        "r_outer": b,
        "voltage": v,
        "sample_radius": radius,
        "electric_field_inner_V_per_m": e_inner,
        "electric_field_outer_V_per_m": e_outer,
        "electric_field_sample_V_per_m": e_sample,
        "pressure_inner_Pa": p_inner,
        "pressure_outer_Pa": p_outer,
        "pressure_ratio": p_inner / p_outer,
        "checks": {
            "positive_C": capacitance > 0.0,
            "energy_matches_half_CV2": abs(energy - 0.5 * capacitance * v * v) <= 1.0e-30,
            "field_scales_as_inverse_r_squared": abs((e_inner / e_outer) - (b / a) ** 2) / ((b / a) ** 2) < 1.0e-12,
            "pressure_ratio_geometry_only": abs((p_inner / p_outer) - (b / a) ** 4) / ((b / a) ** 4) < 1.0e-12,
        },
    }
    if measured_capacitance_F is not None:
        measured_c = float(measured_capacitance_F)
        out["measured_C_F"] = measured_c
        out["measured_C_rel_error"] = abs(measured_c - capacitance) / capacitance
    if measured_energy_J is not None:
        measured_w = float(measured_energy_J)
        out["measured_energy_J"] = measured_w
        out["measured_energy_rel_error"] = abs(measured_w - energy) / energy
    if measured_sample_field_V_per_m is not None:
        measured_e = abs(float(measured_sample_field_V_per_m))
        out["measured_sample_field_V_per_m"] = measured_e
        out["measured_sample_field_rel_error"] = abs(measured_e - e_sample) / e_sample
    return out


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


def layered_parallel_plate_stack_summary(
    area,
    thicknesses,
    eps_r_layers,
    voltage,
    measured_capacitance_F=None,
    measured_interface_voltages=None,
    measured_layer_fields_V_per_m=None,
    measured_energy_J=None,
):
    r"""Readable validation summary for a layered parallel-plate capacitor.

    Layers are stacked normal to the plates, so the normal displacement ``D`` is
    constant through the stack.  This gives

        C = eps0 A / sum_i(d_i/eps_ri),
        D = eps0 V / sum_i(d_i/eps_ri),
        E_i = D/(eps0 eps_ri),
        Delta V_i = E_i d_i.

    The helper records the capacitance, effective permittivity, layer fields,
    voltage drops, interface potentials from the grounded plate, and per-layer
    energy.  Optional measured values add residuals for solver artifacts.
    """

    base = layered_parallel_plate_capacitance(area, thicknesses, eps_r_layers)
    v = float(voltage)
    series = sum(d / e for d, e in zip(thicknesses, eps_r_layers))
    displacement = EPS0 * v / series
    fields = [displacement / (EPS0 * float(e)) for e in eps_r_layers]
    voltage_drops = [field * float(d) for field, d in zip(fields, thicknesses)]
    interface_voltages = []
    acc = 0.0
    for drop in voltage_drops[:-1]:
        acc += drop
        interface_voltages.append(acc)
    energies = [
        0.5 * EPS0 * float(e) * field * field * float(area) * float(d)
        for e, field, d in zip(eps_r_layers, fields, thicknesses)
    ]
    energy = sum(energies)

    out = {
        "policy": "layered_parallel_plate_series_dielectric",
        "C": base["C"],
        "eps_eff": base["eps_eff"],
        "area": float(area),
        "thicknesses": [float(d) for d in thicknesses],
        "eps_r_layers": [float(e) for e in eps_r_layers],
        "voltage": v,
        "normal_displacement_C_per_m2": displacement,
        "layer_fields_V_per_m": fields,
        "layer_voltage_drops_V": voltage_drops,
        "interface_voltages_V": interface_voltages,
        "layer_energies_J": energies,
        "energy_J": energy,
        "checks": {
            "voltage_drops_sum_to_drive": abs(sum(voltage_drops) - v) <= max(1.0e-12, abs(v) * 1.0e-12),
            "layer_energies_sum_to_total": abs(energy - 0.5 * base["C"] * v * v)
            <= max(1.0e-18, abs(energy) * 1.0e-12),
        },
    }
    if measured_capacitance_F is not None:
        measured_c = float(measured_capacitance_F)
        out["measured_capacitance_F"] = measured_c
        out["measured_capacitance_rel_error"] = abs(measured_c - base["C"]) / base["C"]
    if measured_interface_voltages is not None:
        measured = [float(x) for x in measured_interface_voltages]
        if len(measured) != len(interface_voltages):
            raise ValueError("measured_interface_voltages length must match n_layers - 1")
        out["measured_interface_voltages_V"] = measured
        out["measured_interface_voltage_abs_errors_V"] = [
            abs(a - b) for a, b in zip(measured, interface_voltages)
        ]
    if measured_layer_fields_V_per_m is not None:
        measured = [float(x) for x in measured_layer_fields_V_per_m]
        if len(measured) != len(fields):
            raise ValueError("measured_layer_fields_V_per_m length must match n_layers")
        out["measured_layer_fields_V_per_m"] = measured
        out["measured_layer_field_rel_errors"] = [
            abs(abs(a) - b) / b for a, b in zip(measured, fields)
        ]
    if measured_energy_J is not None:
        measured_e = float(measured_energy_J)
        out["measured_energy_J"] = measured_e
        out["measured_energy_rel_error"] = abs(measured_e - energy) / energy
    return out


def parallel_plate_capacitor_energy_force(eps_r, area, gap, voltage):
    r"""Capacitance, stored energy, and plate-attraction force of an ideal
    parallel-plate capacitor (plate ``area``, ``gap`` d, applied ``voltage`` V):

        C = eps0 eps_r area / d ,
        W = 1/2 C V^2 ,
        E = V / d ,
        p = 1/2 eps0 eps_r E^2 ,
        |F| = p area = 1/2 eps0 eps_r area V^2 / d^2 = W / d .

    At fixed voltage the attractive force equals the energy density times the
    plate area, i.e. |F| = W/d, and scales as 1/d^2.  The pressure is the
    normal electrostatic Maxwell traction for a field normal to the plate.
    Returns ``{"C": ..., "energy": ..., "force": ...}`` plus explicit field,
    pressure, and energy-density entries.  (Griffiths.)
    """
    if eps_r <= 0:
        raise ValueError("eps_r must be > 0")
    if area <= 0:
        raise ValueError("area must be > 0")
    if gap <= 0:
        raise ValueError("gap must be > 0")
    field = voltage / gap
    C = EPS0 * eps_r * area / gap
    energy_density = 0.5 * EPS0 * eps_r * field * field
    energy = 0.5 * C * voltage * voltage
    force = energy_density * area
    return {
        "C": C,
        "energy": energy,
        "force": force,
        "electric_field_V_per_m": field,
        "pressure_Pa": energy_density,
        "energy_density_J_per_m3": energy_density,
    }


def capacitance_gradient_force_summary(
    capacitance_F,
    dCdx_F_per_m,
    voltage_V=None,
    charge_C=None,
):
    r"""Electrostatic force from a capacitance gradient.

    For a generalized displacement coordinate ``x`` and capacitance ``C(x)``,
    the signed force along increasing ``x`` is

        fixed voltage:  F_x = 1/2 V^2 dC/dx ,
        fixed charge:   F_x = 1/2 Q^2 C^-2 dC/dx .

    The fixed-voltage expression is the coenergy/source-inclusive result; it is
    the sign convention used in MEMS capacitance-gradient actuators.  If ``x``
    is a gap or height, ``dC/dx`` is usually negative, so the returned force is
    negative (attractive, toward smaller gap).  Provide ``voltage_V``,
    ``charge_C``, or both.  Returns a JSON-friendly summary with whichever
    force routes were requested.
    """

    capacitance = float(capacitance_F)
    gradient = float(dCdx_F_per_m)
    if capacitance <= 0.0:
        raise ValueError("capacitance_F must be > 0")
    if voltage_V is None and charge_C is None:
        raise ValueError("provide voltage_V, charge_C, or both")

    out = {
        "capacitance_F": capacitance,
        "dCdx_F_per_m": gradient,
    }
    if voltage_V is not None:
        voltage = float(voltage_V)
        force = 0.5 * voltage * voltage * gradient
        out.update({
            "voltage_V": voltage,
            "fixed_voltage_force_N": force,
            "fixed_voltage_coenergy_gradient_N": force,
        })
    if charge_C is not None:
        charge = float(charge_C)
        force = 0.5 * charge * charge * gradient / (capacitance * capacitance)
        out.update({
            "charge_C": charge,
            "fixed_charge_force_N": force,
            "fixed_charge_energy_force_N": force,
        })
    if voltage_V is not None and charge_C is not None:
        consistent_charge = capacitance * float(voltage_V)
        out["charge_for_voltage_C"] = consistent_charge
        out["charge_consistency_error_C"] = float(charge_C) - consistent_charge
    return out


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
