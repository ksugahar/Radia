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


def acoustic_method_selection_manifest_gate(manifest, expected_problem_family=None):
    """Validate a public-safe acoustic method-selection lesson.

    The gate is intentionally solver-independent.  It captures reusable
    acoustic modeling choices learned from public acoustics literature and blog
    material, then applies the CAE-AI Lab wave-boundary policy: acoustic and
    electromagnetic exterior/absorbing boundaries use high-order surface
    impedance ``Zs`` and do not use PML as the default validation route.  BEM is
    a frequency-domain exterior-radiation/open-boundary lane; acoustic-structure
    interaction needs an explicit two-way interface; compact thermo/viscous
    submodels should travel as impedance plus power-balance metadata; room
    acoustics should record the modal/high-frequency split; and absorbing
    boundaries should keep local/extended reaction metadata as secondary
    physics, not as the lab open-boundary policy.
    """

    if not isinstance(manifest, dict):
        raise ValueError("manifest must be a dictionary")

    def get(*names, default=None):
        for name in names:
            if name in manifest and manifest[name] is not None:
                return manifest[name]
        return default

    def as_bool(*names):
        value = get(*names, default=False)
        if isinstance(value, str):
            return value.strip().lower() in {"1", "true", "yes", "y", "on"}
        return bool(value)

    def has_text(*names):
        value = get(*names)
        return value is not None and str(value).strip() != ""

    def lower(*names):
        value = get(*names, default="")
        return str(value).strip().lower()

    def positive_number(*names):
        value = get(*names)
        try:
            return value is not None and float(value) > 0.0
        except (TypeError, ValueError):
            return False

    family = lower("problem_family", "family")
    study_domain = lower("study_domain", "analysis_domain")
    primary_method = lower("primary_method", "method")
    exterior_method = lower("exterior_method", "open_boundary_method")
    boundary_model = lower("boundary_model", "boundary_condition")
    reaction_model = lower("reaction_model", "absorber_reaction_model")
    coupling_kind = lower("coupling_kind", "coupling")
    wave_family = lower("wave_family", "physics_family")
    open_boundary_policy = lower("open_boundary_policy", "absorbing_boundary_policy", "radiation_boundary_policy")

    uses_bem = primary_method == "bem" or exterior_method == "bem" or as_bool("uses_bem", "bem")
    uses_time_domain = study_domain == "time_domain" or as_bool("time_domain")
    uses_frequency_domain = study_domain == "frequency_domain" or as_bool("frequency_domain")
    is_wave_boundary = family in {"absorbing_boundary", "exterior_radiation", "wave_open_boundary"} or wave_family in {
        "acoustic",
        "electromagnetic",
        "em",
        "maxwell",
    }

    checks = {
        "problem_family_recorded": has_text("problem_family", "family"),
        "primary_method_recorded": has_text("primary_method", "method"),
        "study_domain_recorded": has_text("study_domain", "analysis_domain"),
        "result_artifact_id_recorded": has_text("result_artifact_id"),
        "result_output_schema_id_recorded": has_text("result_output_schema_id"),
    }
    if expected_problem_family is not None:
        checks["expected_problem_family_matches"] = family == str(expected_problem_family).strip().lower()

    if uses_bem:
        checks["bem_is_frequency_domain"] = uses_frequency_domain and not uses_time_domain
        checks["bem_open_exterior_or_surface_mesh_recorded"] = (
            as_bool("surface_only_boundary_mesh", "boundary_surface_mesh")
            or lower("domain_topology") in {"unbounded_exterior", "open_exterior"}
            or exterior_method == "bem"
        )

    if is_wave_boundary:
        checks["lab_wave_boundary_policy_is_high_order_zs"] = open_boundary_policy in {
            "high_order_zs",
            "high_order_surface_impedance",
            "higher_order_zs",
            "higher_order_surface_impedance",
        }
        checks["pml_not_used"] = not as_bool("uses_pml", "pml", "absorbing_layer")

    if family == "acoustic_structure_interaction":
        checks["asi_coupling_is_two_way"] = "two" in coupling_kind or coupling_kind in {
            "fem_bem",
            "solid_fluid",
            "fluid_structure",
        }
        checks["asi_structural_field_recorded"] = has_text("structural_field")
        checks["asi_acoustic_field_recorded"] = has_text("acoustic_field")

    if family == "impedance_lumping":
        checks["impedance_kind_recorded"] = has_text("impedance_kind", "impedance_model")
        checks["frequency_dependent_impedance_recorded"] = as_bool(
            "frequency_dependent_impedance", "impedance_frequency_dependent"
        )
        checks["power_balance_observable_recorded"] = has_text("power_balance_observable")

    if family == "room_acoustics":
        checks["schroeder_frequency_recorded"] = positive_number("schroeder_frequency_hz")
        checks["low_frequency_wave_method_recorded"] = has_text("low_frequency_method")
        checks["high_frequency_method_recorded"] = has_text("high_frequency_method")

    if family == "absorbing_boundary":
        checks["absorber_reaction_model_recorded"] = reaction_model in {"local", "extended"}
        checks["angle_dependency_recorded"] = as_bool("angle_dependency_recorded", "incident_angle_dependency")
        checks["boundary_model_is_high_order_zs"] = boundary_model in {
            "high_order_zs",
            "high_order_surface_impedance",
            "higher_order_zs",
            "higher_order_surface_impedance",
        }
        if reaction_model == "extended":
            checks["extended_reaction_metadata_recorded"] = as_bool(
                "frequency_dependent_boundary",
                "frequency_dependent_impedance",
            )

    return {
        "policy": "acoustic_method_selection_manifest_gate",
        "problem_family": family or None,
        "primary_method": primary_method or None,
        "exterior_method": exterior_method or None,
        "study_domain": study_domain or None,
        "boundary_model": boundary_model or None,
        "reaction_model": reaction_model or None,
        "wave_family": wave_family or None,
        "open_boundary_policy": open_boundary_policy or None,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "teaching_note": (
            "Use this before promoting public acoustic blog/literature lessons "
            "into FEM/BEM examples: method family, domain, coupling, impedance, "
            "absorber reaction, and output schema identity must travel together."
        ),
    }


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


def helmholtz_green_low_frequency_teaching_report(distance, wavenumber, order=6):
    """Explain the readable low-frequency Helmholtz BEM split.

    This wraps :func:`helmholtz_green_low_frequency_series` in a teaching-shaped
    dictionary.  The important lesson is not speed; it is that the singular
    Laplace kernel and the smooth Helmholtz correction should be kept visible
    when students read a FEM/BEM coupling code.
    """

    series = helmholtz_green_low_frequency_series(distance, wavenumber, order=order)
    laplace = series["laplace_term"]
    stable_correction = series["regular_part"]
    direct_correction = series["exact"] - laplace
    correction_scale = abs(stable_correction)
    cancellation_ratio = math.inf if correction_scale == 0.0 else abs(laplace) / correction_scale
    return {
        "kind": "low_frequency_helmholtz_teaching_report",
        "policy": "readable_bem_kernel_split_not_production_quadrature",
        "time_convention": "exp(+i omega t), outgoing exp(-i k r)",
        "distance": series["distance"],
        "wavenumber": series["wavenumber"],
        "order": series["order"],
        "kr_abs": series["kr_abs"],
        "laplace_term": laplace,
        "stable_correction": stable_correction,
        "direct_correction": direct_correction,
        "single_layer": series["approx"],
        "direct_green": series["exact"],
        "stable_error": series["abs_error"],
        "correction_agreement": abs(stable_correction - direct_correction),
        "cancellation_ratio": cancellation_ratio,
        "notes": [
            "G_k = G_0 + (exp(-1i*k*r)-1)/(4*pi*r)",
            "G_0 keeps the singular Laplace quadrature visible",
            "the correction is smooth and can be evaluated by a Taylor split at low k*r",
        ],
    }


def low_frequency_helmholtz_kernel_manifest_gate(
    report,
    expected_kernel_family=None,
    expected_low_frequency_strategy=None,
    expected_time_convention=None,
    max_kr_abs=None,
    max_stable_error=1.0e-12,
    max_correction_agreement=1.0e-12,
    min_cancellation_ratio=None,
):
    """Check a readable low-frequency Helmholtz BEM kernel report.

    This is a manifest gate for teaching code, not a production quadrature
    checker.  It verifies that the singular Laplace part, smooth Helmholtz
    correction, time convention, kernel family, and low-frequency strategy are
    recorded before a MATLAB/Gypsilab or NGSolve BEM row is reused.
    """

    if not isinstance(report, dict):
        raise ValueError("report must be a dictionary")

    def first(names):
        for name in names:
            if name in report and report[name] is not None:
                return report[name]
        return None

    def as_float(value, default=None):
        if value is None:
            return default
        return float(value)

    def has_value(value):
        return value is not None and str(value).strip() != ""

    kind = str(first(("kind", "report_kind")) or "")
    policy = str(first(("policy", "gate_policy")) or "")
    kernel_family = first(("kernel_family", "kernelFamily", "bem_kernel_family", "bemKernelFamily"))
    strategy = first((
        "low_frequency_strategy",
        "lowFrequencyStrategy",
        "stabilization_strategy",
        "kernel_split_strategy",
    ))
    time_convention = first(("time_convention", "timeConvention"))
    kr_abs = as_float(first(("kr_abs", "krAbs")), None)
    if kr_abs is None:
        kr = first(("kr", "kTimesR"))
        kr_abs = None if kr is None else abs(complex(kr))
    stable_error = as_float(first(("stable_error", "stableError")), math.inf)
    correction_agreement = as_float(first(("correction_agreement", "correctionAgreement")), math.inf)
    cancellation_ratio = as_float(first(("cancellation_ratio", "cancellationRatio")), 0.0)
    laplace_term = first(("laplace_term", "laplaceTerm"))
    stable_correction = first(("stable_correction", "stableCorrection"))
    direct_correction = first(("direct_correction", "directCorrection"))

    checks = {
        "kind_is_low_frequency_report": kind == "low_frequency_helmholtz_teaching_report",
        "policy_records_readable_split": "kernel_split" in policy or "bem_kernel_split" in policy,
        "kernel_family_recorded": has_value(kernel_family),
        "low_frequency_strategy_recorded": has_value(strategy),
        "time_convention_recorded": has_value(time_convention),
        "kr_abs_recorded": kr_abs is not None,
        "laplace_term_recorded": laplace_term is not None,
        "stable_correction_recorded": stable_correction is not None,
        "direct_correction_recorded": direct_correction is not None,
        "stable_error_within_tolerance": stable_error <= float(max_stable_error),
        "correction_agreement_within_tolerance": correction_agreement <= float(max_correction_agreement),
    }
    if expected_kernel_family is not None:
        checks["expected_kernel_family_matches"] = str(kernel_family or "") == str(expected_kernel_family)
    if expected_low_frequency_strategy is not None:
        checks["expected_low_frequency_strategy_matches"] = str(strategy or "") == str(expected_low_frequency_strategy)
    if expected_time_convention is not None:
        checks["expected_time_convention_matches"] = str(time_convention or "") == str(expected_time_convention)
    if max_kr_abs is not None:
        checks["kr_abs_within_low_frequency_limit"] = kr_abs is not None and kr_abs <= float(max_kr_abs)
    if min_cancellation_ratio is not None:
        checks["cancellation_ratio_large_enough"] = cancellation_ratio >= float(min_cancellation_ratio)

    return {
        "policy": "low_frequency_helmholtz_kernel_manifest_gate",
        "kind": kind,
        "kernel_family": None if kernel_family is None else str(kernel_family),
        "low_frequency_strategy": None if strategy is None else str(strategy),
        "time_convention": None if time_convention is None else str(time_convention),
        "kr_abs": kr_abs,
        "stable_error": stable_error,
        "correction_agreement": correction_agreement,
        "cancellation_ratio": cancellation_ratio,
        "expected_kernel_family": None if expected_kernel_family is None else str(expected_kernel_family),
        "expected_low_frequency_strategy": (
            None if expected_low_frequency_strategy is None else str(expected_low_frequency_strategy)
        ),
        "expected_time_convention": None if expected_time_convention is None else str(expected_time_convention),
        "max_kr_abs": None if max_kr_abs is None else float(max_kr_abs),
        "max_stable_error": float(max_stable_error),
        "max_correction_agreement": float(max_correction_agreement),
        "min_cancellation_ratio": None if min_cancellation_ratio is None else float(min_cancellation_ratio),
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use before low-frequency Helmholtz BEM teaching rows are reused: "
            "kernel family, time convention, Laplace singular term, smooth "
            "correction strategy, and kr limit must be explicit."
        ),
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


def acoustic_boundary_power_summary(
    pressure,
    normal_velocity,
    area=1.0,
    amplitude="peak",
):
    r"""Active/reactive acoustic power from boundary pressure and normal velocity.

    For complex pressure ``p`` and outward normal velocity ``v_n`` on a boundary
    patch, the complex normal intensity is

        I_n = alpha p conj(v_n),

    where ``alpha=0.5`` for peak phasors and ``alpha=1`` for RMS phasors.  The
    real part is active outward power density and the imaginary part is the
    reactive near-field exchange.  This is the mesh-postprocessing scalar that
    acoustic FEM/BEM models use after solving for pressure and boundary
    velocity traces.
    """

    p = complex(pressure)
    v = complex(normal_velocity)
    if not (
        math.isfinite(p.real)
        and math.isfinite(p.imag)
        and math.isfinite(v.real)
        and math.isfinite(v.imag)
    ):
        raise ValueError("pressure and normal_velocity must be finite")
    patch_area = float(area)
    if patch_area < 0.0:
        raise ValueError("area must be >= 0")
    if amplitude == "peak":
        factor = 0.5
    elif amplitude == "rms":
        factor = 1.0
    else:
        raise ValueError("amplitude must be 'peak' or 'rms'")

    intensity = factor * p * v.conjugate()
    return {
        "pressure": p,
        "normal_velocity": v,
        "area": patch_area,
        "amplitude": amplitude,
        "phasor_average_factor": factor,
        "complex_intensity": intensity,
        "active_intensity": intensity.real,
        "reactive_intensity": intensity.imag,
        "active_power": patch_area * intensity.real,
        "reactive_power": patch_area * intensity.imag,
        "apparent_intensity": abs(intensity),
        "apparent_power": patch_area * abs(intensity),
        "specific_impedance": None if v == 0.0 else p / v,
        "specific_admittance": None if p == 0.0 else v / p,
        "policy": "outward_active_power_positive_for_pressure_times_conjugate_normal_velocity",
    }


def acoustic_impedance_reflection_summary(
    specific_impedance,
    incidence_angle_rad=0.0,
    incident_pressure=1.0,
    rho=1.2041,
    c=343.0,
    amplitude="peak",
):
    r"""Plane-wave reflection and absorption at a local acoustic impedance load.

    This is the one-port companion to :func:`acoustic_boundary_power_summary`.
    A propagating plane wave is incident on a locally reacting boundary whose
    load impedance is measured with velocity positive *into* the load.  The
    pressure reflection coefficient is

        Gamma = (z_load - z_n) / (z_load + z_n),

    where ``z_n = rho*c/cos(theta)`` is the plane-wave normal impedance.  The
    active absorption coefficient is ``1 - |Gamma|^2``.  Peak phasors use the
    0.5 time-average factor; RMS phasors use 1.0.

    The return dictionary includes a direct power balance:

        incident_intensity - reflected_intensity == absorbed_intensity

    and the same absorbed complex intensity computed from total boundary
    pressure and load normal velocity.  This is a compact sign/conjugation
    check for acoustic FEM/BEM impedance boundaries.
    """

    rrho = float(rho)
    cc = float(c)
    theta = float(incidence_angle_rad)
    if rrho <= 0.0:
        raise ValueError("rho must be > 0")
    if cc <= 0.0:
        raise ValueError("c must be > 0")
    if abs(theta) >= 0.5 * math.pi:
        raise ValueError("incidence_angle_rad must be strictly between -pi/2 and pi/2")
    if amplitude == "peak":
        factor = 0.5
    elif amplitude == "rms":
        factor = 1.0
    else:
        raise ValueError("amplitude must be 'peak' or 'rms'")

    z_load = complex(specific_impedance)
    if not (
        math.isfinite(z_load.real)
        and math.isfinite(z_load.imag)
    ):
        raise ValueError("specific_impedance must be finite")
    if z_load == 0.0:
        gamma = -1.0 + 0.0j
    else:
        z_normal = rrho * cc / math.cos(theta)
        denom = z_load + z_normal
        if denom == 0.0:
            raise ValueError("specific_impedance + characteristic normal impedance must be nonzero")
        gamma = (z_load - z_normal) / denom

    z_normal = rrho * cc / math.cos(theta)
    p_inc = complex(incident_pressure)
    if not (math.isfinite(p_inc.real) and math.isfinite(p_inc.imag)):
        raise ValueError("incident_pressure must be finite")
    p_ref = gamma * p_inc
    p_total = p_inc + p_ref
    velocity_into_load = (p_inc - p_ref) / z_normal
    boundary_intensity = factor * p_total * velocity_into_load.conjugate()
    incident_intensity = factor * abs(p_inc) ** 2 / z_normal
    reflected_intensity = factor * abs(p_ref) ** 2 / z_normal
    absorbed_intensity = incident_intensity - reflected_intensity
    absorption = 1.0 - abs(gamma) ** 2

    return {
        "rho": rrho,
        "c": cc,
        "incidence_angle_rad": theta,
        "amplitude": amplitude,
        "phasor_average_factor": factor,
        "characteristic_normal_impedance": z_normal,
        "specific_impedance": z_load,
        "normalized_impedance": z_load / z_normal,
        "incident_pressure": p_inc,
        "reflected_pressure": p_ref,
        "total_boundary_pressure": p_total,
        "normal_velocity_into_load": velocity_into_load,
        "pressure_reflection_coefficient": gamma,
        "velocity_reflection_coefficient": -gamma,
        "power_reflection_coefficient": abs(gamma) ** 2,
        "absorption_coefficient": absorption,
        "incident_intensity": incident_intensity,
        "reflected_intensity": reflected_intensity,
        "absorbed_intensity": absorbed_intensity,
        "boundary_complex_intensity_into_load": boundary_intensity,
        "boundary_active_intensity_into_load": boundary_intensity.real,
        "boundary_reactive_intensity_into_load": boundary_intensity.imag,
        "power_balance_residual": incident_intensity - reflected_intensity - boundary_intensity.real,
        "policy": "velocity_positive_into_load_pressure_reflection_gamma_zload_minus_zn_over_zload_plus_zn",
    }


def acoustic_impedance_radiation_pressure_summary(
    specific_impedance,
    area=1.0,
    incidence_angle_rad=0.0,
    incident_pressure=1.0,
    rho=1.2041,
    c=343.0,
    amplitude="peak",
):
    r"""Normal acoustic momentum pressure from an impedance reflection summary.

    The incident intensity returned by :func:`acoustic_impedance_reflection_summary`
    is the normal energy flux into the boundary.  The corresponding normal
    momentum transfer to the load is

        pressure = (1 + R) I_inc / c = (A + 2 R) I_inc / c,

    where ``R=|Gamma|^2`` and ``A=1-R`` for a passive one-port load.  Thus a
    matched absorber gives ``I/c`` and a lossless reflector gives ``2 I/c``.
    """

    patch_area = float(area)
    if patch_area < 0.0:
        raise ValueError("area must be >= 0")
    reflection = acoustic_impedance_reflection_summary(
        specific_impedance,
        incidence_angle_rad=incidence_angle_rad,
        incident_pressure=incident_pressure,
        rho=rho,
        c=c,
        amplitude=amplitude,
    )
    speed = float(reflection["c"])
    incident_intensity = float(reflection["incident_intensity"])
    reflectance = float(reflection["power_reflection_coefficient"])
    absorption = float(reflection["absorption_coefficient"])
    pressure = (1.0 + reflectance) * incident_intensity / speed
    equivalent_pressure = (absorption + 2.0 * reflectance) * incident_intensity / speed
    absorbed_pressure = absorption * incident_intensity / speed
    reflected_pressure = 2.0 * reflectance * incident_intensity / speed
    return {
        "area": patch_area,
        "reflection": reflection,
        "incident_normal_intensity": incident_intensity,
        "power_reflection_coefficient": reflectance,
        "absorption_coefficient": absorption,
        "momentum_transfer_factor": 1.0 + reflectance,
        "absorber_reflector_equivalent_factor": absorption + 2.0 * reflectance,
        "absorbed_momentum_pressure_Pa": absorbed_pressure,
        "reflected_momentum_pressure_Pa": reflected_pressure,
        "normal_momentum_pressure_Pa": pressure,
        "normal_momentum_pressure_equivalent_Pa": equivalent_pressure,
        "normal_force_N": patch_area * pressure,
        "force_from_absorptance_reflectance_N": patch_area * equivalent_pressure,
        "force_balance_residual_N": patch_area * (pressure - equivalent_pressure),
        "passive_one_port": absorption >= -1.0e-12,
        "policy": "acoustic_impedance_momentum_pressure_from_absorptance_and_reflectance",
    }


def acoustic_impedance_reflection_sweep_summary(
    frequency_Hz,
    specific_impedance_values,
    area=1.0,
    incidence_angle_rad=0.0,
    incident_pressure=1.0,
    rho=1.2041,
    c=343.0,
    amplitude="peak",
    passivity_tolerance=1.0e-12,
):
    """Audit acoustic impedance reflection and momentum pressure over a sweep."""

    frequencies = [float(value) for value in frequency_Hz]
    impedances = [complex(value) for value in specific_impedance_values]
    if len(frequencies) != len(impedances):
        raise ValueError("frequency_Hz and specific_impedance_values must have the same length")
    if not frequencies:
        raise ValueError("at least one frequency sample is required")
    tolerance = float(passivity_tolerance)
    if tolerance < 0.0:
        raise ValueError("passivity_tolerance must be >= 0")

    rows = []
    violation_rows = []
    for idx, (frequency, impedance) in enumerate(zip(frequencies, impedances)):
        if not math.isfinite(frequency) or frequency < 0.0:
            raise ValueError("frequency samples must be finite and >= 0")
        if not math.isfinite(impedance.real) or not math.isfinite(impedance.imag):
            raise ValueError("specific impedance values must be finite")
        pressure = acoustic_impedance_radiation_pressure_summary(
            impedance,
            area=area,
            incidence_angle_rad=incidence_angle_rad,
            incident_pressure=incident_pressure,
            rho=rho,
            c=c,
            amplitude=amplitude,
        )
        reflection = pressure["reflection"]
        gamma = complex(reflection["pressure_reflection_coefficient"])
        absorption = float(pressure["absorption_coefficient"])
        row = {
            "index": idx,
            "frequency_Hz": frequency,
            "specific_impedance_real": impedance.real,
            "specific_impedance_imag": impedance.imag,
            "normalized_impedance_real": complex(reflection["normalized_impedance"]).real,
            "normalized_impedance_imag": complex(reflection["normalized_impedance"]).imag,
            "pressure_reflection_real": gamma.real,
            "pressure_reflection_imag": gamma.imag,
            "pressure_reflection_magnitude": abs(gamma),
            "pressure_reflection_phase_rad": math.atan2(gamma.imag, gamma.real),
            "power_reflection_coefficient": float(pressure["power_reflection_coefficient"]),
            "absorption_coefficient": absorption,
            "momentum_transfer_factor": float(pressure["momentum_transfer_factor"]),
            "incident_normal_intensity": float(pressure["incident_normal_intensity"]),
            "normal_momentum_pressure_Pa": float(pressure["normal_momentum_pressure_Pa"]),
            "normal_force_N": float(pressure["normal_force_N"]),
            "passivity_excess_absorption": max(0.0, -absorption),
            "passivity_ok": absorption >= -tolerance,
        }
        rows.append(row)
        if not row["passivity_ok"]:
            violation_rows.append(row)

    max_absorption_row = max(rows, key=lambda row: row["absorption_coefficient"])
    min_absorption_row = min(rows, key=lambda row: row["absorption_coefficient"])
    max_force_row = max(rows, key=lambda row: row["normal_force_N"])
    min_force_row = min(rows, key=lambda row: row["normal_force_N"])
    monotonic = all(
        frequencies[idx] < frequencies[idx + 1]
        for idx in range(len(frequencies) - 1)
    )

    return {
        "policy": "acoustic_impedance_reflection_sweep_momentum_audit",
        "n_points": len(rows),
        "frequency_min_Hz": min(frequencies),
        "frequency_max_Hz": max(frequencies),
        "frequency_monotonic_increasing": monotonic,
        "area": float(area),
        "incidence_angle_rad": float(incidence_angle_rad),
        "rho": float(rho),
        "c": float(c),
        "amplitude": amplitude,
        "passivity_tolerance": tolerance,
        "passivity_ok": not violation_rows,
        "passivity_violation_count": len(violation_rows),
        "max_passivity_excess_absorption": max(row["passivity_excess_absorption"] for row in rows),
        "mean_absorption_coefficient": sum(row["absorption_coefficient"] for row in rows) / len(rows),
        "max_absorption_coefficient": max_absorption_row["absorption_coefficient"],
        "max_absorption_frequency_Hz": max_absorption_row["frequency_Hz"],
        "min_absorption_coefficient": min_absorption_row["absorption_coefficient"],
        "min_absorption_frequency_Hz": min_absorption_row["frequency_Hz"],
        "mean_normal_force_N": sum(row["normal_force_N"] for row in rows) / len(rows),
        "max_normal_force_N": max_force_row["normal_force_N"],
        "max_force_frequency_Hz": max_force_row["frequency_Hz"],
        "min_normal_force_N": min_force_row["normal_force_N"],
        "min_force_frequency_Hz": min_force_row["frequency_Hz"],
        "force_span_N": max_force_row["normal_force_N"] - min_force_row["normal_force_N"],
        "max_force_row": max_force_row,
        "min_force_row": min_force_row,
        "max_absorption_row": max_absorption_row,
        "passivity_violation_rows": violation_rows,
        "status": "ok" if not violation_rows else "needs_attention",
        "rows": rows,
    }


def _bessel_j1_fallback(x):
    """Small dependency-free J1 approximation for acoustic piston gates."""

    value = float(x)
    ax = abs(value)
    if ax == 0.0:
        return 0.0
    if ax <= 20.0:
        term = value / 2.0
        total = term
        x2_over4 = value * value / 4.0
        for k in range(1, 80):
            term *= -x2_over4 / (k * (k + 1.0))
            total += term
            if abs(term) <= max(1.0e-18, 1.0e-16 * abs(total)):
                break
        return total

    phase = ax - 0.75 * math.pi
    root = math.sqrt(2.0 / (math.pi * ax))
    approx = root * (
        math.cos(phase)
        - 3.0 * math.sin(phase) / (8.0 * ax)
        + 15.0 * math.cos(phase) / (128.0 * ax * ax)
    )
    return -approx if value < 0.0 else approx


def _struve_h1_fallback(x):
    """Small dependency-free H1 approximation for acoustic piston gates."""

    value = float(x)
    ax = abs(value)
    if ax == 0.0:
        return 0.0
    if ax <= 20.0:
        term = 2.0 * ax * ax / (3.0 * math.pi)
        total = term
        x2_over4 = ax * ax / 4.0
        for k in range(80):
            term *= -x2_over4 / ((k + 1.5) * (k + 2.5))
            total += term
            if abs(term) <= max(1.0e-18, 1.0e-16 * abs(total)):
                break
        return total

    # H1(x) approaches Y1(x) + 2/pi.  This is sufficient for the high-ka
    # sanity gate where the piston resistance tends to rho*c and the reactance
    # is small compared with the plane-wave resistance.
    phase = ax - 0.75 * math.pi
    root = math.sqrt(2.0 / (math.pi * ax))
    y1 = root * (
        math.sin(phase)
        + 3.0 * math.cos(phase) / (8.0 * ax)
        - 15.0 * math.sin(phase) / (128.0 * ax * ax)
    )
    return 2.0 / math.pi + y1


def _baffled_piston_resistance_reactance_ratios(ka, prefer_scipy=True):
    """Return piston resistance/reactance ratios, with SciPy optional."""

    value = float(ka)
    if value <= 0.0:
        raise ValueError("ka must be > 0")
    x = 2.0 * value
    if prefer_scipy:
        try:
            from scipy.special import j1, struve

            return 1.0 - float(j1(x)) / value, float(struve(1, x)) / value, "scipy"
        except ModuleNotFoundError:
            pass

    return (
        1.0 - _bessel_j1_fallback(x) / value,
        _struve_h1_fallback(x) / value,
        "fallback",
    )


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

    velocity = complex(surface_velocity)
    omega = 2.0 * math.pi * f
    k = omega / cc
    ka = k * a
    area = math.pi * a * a
    resistance_ratio, reactance_ratio, special_function_source = (
        _baffled_piston_resistance_reactance_ratios(ka)
    )
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
        "special_function_source": special_function_source,
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
