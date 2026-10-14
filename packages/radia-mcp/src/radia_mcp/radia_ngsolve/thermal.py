"""Closed-form thermal post-processors for radia-ngsolve -- lumped thermal
networks, contact / constriction resistance, induction-heating surface power,
and convection / radiation heat-transfer coefficients.

These are the analytic (textbook) companions to a field solve: once a magnetic
or eddy-current FE run has produced a loss density, a lumped thermal network or
a closed-form heat-transfer coefficient turns that loss into a temperature rise
without a second mesh. All quantities are SI.

References
----------
* F. P. Incropera & D. P. DeWitt, *Fundamentals of Heat and Mass Transfer*
  (lumped capacitance Ch.5; external/internal convection Ch.7-8; radiation
  exchange Ch.13).
* M. M. Yovanovich, "Four decades of research on thermal contact, gap, and
  joint resistance in microelectronics," *IEEE Trans. CPMT* 28(2), 2005;
  H. S. Carslaw & J. C. Jaeger, *Conduction of Heat in Solids*, 2nd ed.
  (constriction / spreading resistance of a circular contact).
* V. Rudnev et al., *Handbook of Induction Heating*, 2nd ed. (skin depth,
  surface impedance, surface power density).
* F. W. Dittus & L. M. K. Boelter, *Univ. Calif. Publ. Eng.* 2, 443 (1930)
  (turbulent pipe-flow Nusselt number).
* S. W. Churchill & H. H. S. Chu, "Correlating equations for laminar and
  turbulent free convection from a vertical plate," *Int. J. Heat Mass
  Transfer* 18, 1323 (1975).
"""
from __future__ import annotations

import math

MU0 = 4.0e-7 * math.pi              # vacuum permeability [H/m]
SIGMA_SB = 5.670374419e-8           # Stefan-Boltzmann constant [W/(m^2 K^4)]


def lumped_thermal_time_constant(k, area, length, h_film, rho, c_p, power):
    """Lumped (single-node) thermal RC network for a conducting slab of
    cross-section ``area`` A [m^2] and ``length`` L [m] dissipating ``power``
    P [W], cooled through a surface film of coefficient ``h_film`` h [W/(m^2 K)].

    Conduction and film resistances add in series, the capacitance is the bulk
    heat store, and the response is first-order (lumped-capacitance model,
    Incropera Ch.5):

        R_cond = L / (k A) ,   R_film = 1 / (h A) ,   R_th = R_cond + R_film ,
        C_th   = rho c_p A L ,
        dT_ss  = P R_th ,       tau = R_th C_th ,
        T(t) - T_inf = dT_ss (1 - exp(-t/tau)) ,

    so after one time constant the rise has reached ``frac_at_tau = 1 - 1/e``
    (~63.2 %) of its steady value. ``k`` [W/(m K)], ``rho`` [kg/m^3],
    ``c_p`` [J/(kg K)].

    Returns ``{R_th, C_th, dT_steady, tau, frac_at_tau}``.
    """
    if k <= 0 or area <= 0 or length <= 0 or h_film <= 0 or rho <= 0 or c_p <= 0:
        raise ValueError("k, area, length, h_film, rho, c_p must all be positive")
    if power < 0:
        raise ValueError("power must be non-negative")
    r_cond = length / (k * area)
    r_film = 1.0 / (h_film * area)
    r_th = r_cond + r_film
    c_th = rho * c_p * area * length
    return {
        "R_th": r_th,
        "C_th": c_th,
        "dT_steady": power * r_th,
        "tau": r_th * c_th,
        "frac_at_tau": 1.0 - math.exp(-1.0),
    }


def thermal_constriction_resistance(k, radius_a, k2=None):
    """Spreading / constriction resistance of heat flowing into a half-space
    through a circular contact spot of radius ``radius_a`` a [m] in a body of
    conductivity ``k`` [W/(m K)] (Carslaw & Jaeger; Yovanovich 2005):

        isothermal contact : R = 1 / (4 k a)         (uniform-temperature disc)
        isoflux  contact   : R = 8 / (3 pi^2 k a)    (uniform-flux disc)

    The isoflux spot is slightly more resistive; the ratio is the pure number
    ``constriction_factor = R_isoflux / R_isothermal = 32 / (3 pi^2) ~ 1.081``.
    When a second body of conductivity ``k2`` is given, the two constriction
    resistances add in series across the interface (two-body joint):

        R_two_body = (1 / (4 a)) (1/k + 1/k2) .

    Returns ``{R_isothermal, R_isoflux, constriction_factor}`` plus
    ``R_two_body`` when ``k2`` is supplied.
    """
    if k <= 0 or radius_a <= 0:
        raise ValueError("k and radius_a must be positive")
    if k2 is not None and k2 <= 0:
        raise ValueError("k2 must be positive when given")
    r_iso = 1.0 / (4.0 * k * radius_a)
    r_flux = 8.0 / (3.0 * math.pi ** 2 * k * radius_a)
    out = {
        "R_isothermal": r_iso,
        "R_isoflux": r_flux,
        "constriction_factor": r_flux / r_iso,
    }
    if k2 is not None:
        out["R_two_body"] = (1.0 / (4.0 * radius_a)) * (1.0 / k + 1.0 / k2)
    return out


def induction_heating_surface_power(sigma, mu_r, freq, h_surface):
    """Surface power density absorbed by a thick conductor under an external
    tangential magnetic field of RMS amplitude ``h_surface`` H_s [A/m]
    (surface-impedance / standard-penetration model, Rudnev *Handbook of
    Induction Heating*):

        delta = sqrt(2 / (2 pi f mu0 mu_r sigma))     reference (skin) depth
        R_s   = 1 / (sigma delta)                     surface resistance
        P/A   = R_s H_s^2 = H_s^2 sqrt(pi f mu0 mu_r / sigma) .

    The two forms are algebraically identical (substituting delta into R_s),
    so the power density scales as sqrt(f): the second form makes the
    dependence explicit. ``sigma`` [S/m], ``freq`` f [Hz].

    Returns ``{power_per_area, reference_depth, Rs}``.
    """
    if sigma <= 0 or mu_r <= 0 or freq <= 0:
        raise ValueError("sigma, mu_r and freq must be positive")
    if h_surface < 0:
        raise ValueError("h_surface must be non-negative")
    delta = math.sqrt(2.0 / (2.0 * math.pi * freq * MU0 * mu_r * sigma))
    r_s = 1.0 / (sigma * delta)
    return {
        "power_per_area": r_s * h_surface ** 2,
        "reference_depth": delta,
        "Rs": r_s,
    }


def forced_convection_dittus_boelter(reynolds, prandtl, k_fluid,
                                     hydraulic_diameter, heating=True):
    """Dittus-Boelter turbulent internal-flow correlation for the convective
    coefficient in a duct of hydraulic diameter ``hydraulic_diameter`` D [m]
    (Dittus & Boelter 1930; Incropera Ch.8):

        Nu = 0.023 Re^0.8 Pr^n ,   n = 0.4 (heating)  or  0.3 (cooling),
        h  = Nu k_fluid / D .

    Valid for fully developed turbulent flow (Re >~ 10^4, 0.6 <~ Pr <~ 160);
    ``n`` is 0.4 when the wall heats the fluid and 0.3 when it cools the fluid.
    ``k_fluid`` [W/(m K)].

    Returns ``{Nu, h}``.
    """
    if reynolds <= 0 or prandtl <= 0 or k_fluid <= 0 or hydraulic_diameter <= 0:
        raise ValueError("reynolds, prandtl, k_fluid, hydraulic_diameter must be positive")
    n = 0.4 if heating else 0.3
    nu = 0.023 * reynolds ** 0.8 * prandtl ** n
    return {"Nu": nu, "h": nu * k_fluid / hydraulic_diameter}


def natural_convection_vertical_plate(rayleigh, prandtl, k_fluid, height):
    """Churchill-Chu correlation for free (natural) convection from an
    isothermal vertical plate of ``height`` L [m] (Churchill & Chu 1975),
    valid over the whole Rayleigh range:

        Nu = { 0.825 + 0.387 Ra^(1/6) /
               [ 1 + (0.492 / Pr)^(9/16) ]^(8/27) }^2 ,
        h  = Nu k_fluid / L .

    As Ra -> 0 the bracket collapses to the constant 0.825, so Nu -> 0.825^2
    ~ 0.68 (the conduction floor of the correlation). ``k_fluid`` [W/(m K)].

    Returns ``{Nu, h}``.
    """
    if rayleigh < 0 or prandtl <= 0 or k_fluid <= 0 or height <= 0:
        raise ValueError("rayleigh must be >=0; prandtl, k_fluid, height must be positive")
    denom = (1.0 + (0.492 / prandtl) ** (9.0 / 16.0)) ** (8.0 / 27.0)
    nu = (0.825 + 0.387 * rayleigh ** (1.0 / 6.0) / denom) ** 2
    return {"Nu": nu, "h": nu * k_fluid / height}


def radiation_exchange_two_surface(emissivity_1, area_1, T1, T2,
                                   emissivity_2=1.0, area_2=None,
                                   view_factor=1.0):
    """Net radiative heat exchange between two gray diffuse surfaces at absolute
    temperatures ``T1`` and ``T2`` [K] (Incropera Ch.13). For the general
    two-surface gray enclosure the net rate is

        Q = sigma_SB (T1^4 - T2^4) /
            [ (1-e1)/(e1 A1) + 1/(A1 F12) + (1-e2)/(e2 A2) ] .

    When ``area_2`` is omitted the second body is treated as a large isothermal
    enclosure (e2 -> effectively 1, A2 -> infinity, F12 = 1), reducing to the
    small-object limit

        Q = e1 sigma_SB A1 (T1^4 - T2^4) .

    The exchange is linearised about the mean temperature T_m = (T1+T2)/2 to a
    radiation film coefficient

        h_rad = 4 e1 sigma_SB T_m^3 ,

    so Q vanishes when T1 == T2. ``area_1``, ``area_2`` [m^2].

    Returns ``{Q, h_rad_linearized}``.
    """
    if area_1 <= 0:
        raise ValueError("area_1 must be positive")
    if not (0.0 < emissivity_1 <= 1.0) or not (0.0 < emissivity_2 <= 1.0):
        raise ValueError("emissivities must lie in (0, 1]")
    if T1 < 0 or T2 < 0:
        raise ValueError("absolute temperatures must be non-negative")
    if not (0.0 < view_factor <= 1.0):
        raise ValueError("view_factor must lie in (0, 1]")
    if area_2 is not None and area_2 <= 0:
        raise ValueError("area_2 must be positive when given")

    dpow = T1 ** 4 - T2 ** 4
    if area_2 is None:
        # small object in a large enclosure
        q = emissivity_1 * SIGMA_SB * area_1 * dpow
    else:
        denom = ((1.0 - emissivity_1) / (emissivity_1 * area_1)
                 + 1.0 / (area_1 * view_factor)
                 + (1.0 - emissivity_2) / (emissivity_2 * area_2))
        q = SIGMA_SB * dpow / denom
    t_m = 0.5 * (T1 + T2)
    return {"Q": q, "h_rad_linearized": 4.0 * emissivity_1 * SIGMA_SB * t_m ** 3}
