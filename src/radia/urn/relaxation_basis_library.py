"""
Relaxation Basis Function Library (Circuit-Compatible Version)

A curated collection of basis functions that have DIRECT circuit equivalents.
Every function has a documented RC/RL/RLC realization (Foster networks or
fixed log-spaced fractional-order approximations).

Target applications:
- SPICE model extraction
- Electromagnetic device modeling
- Power electronics simulation

Key principle: Every function here maps to a realizable circuit.

Circuit Equivalents:
- Debye: RC parallel
- Cole-Cole: RC ladder (Foster network approximation)
- CPE: RC ladder (Valsa/Charef approximation)
- Warburg: RC ladder (finite approximation)
- Skin effect: RL ladder (Dowell coefficients)
- RLC: Direct RLC elements
- Maxwell/Voigt: Spring-dashpot -> RL/RC analog

References:
[1] Valsa & Vlach, "RC models of fractional-order elements", IJCT, 2013
[2] Charef et al., "Fractional order systems approximation", IJCTA, 2006
[3] Dowell, "Effects of eddy currents in transformer windings", PROC IEE, 1966
[4] Foster, "A reactance theorem", Bell Syst. Tech. J., 1924
"""

import numpy as np
import torch
from typing import Dict, List, Tuple, Callable
from dataclasses import dataclass


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def safe_power(base: torch.Tensor, exp: torch.Tensor) -> torch.Tensor:
    """Safe complex power for fractional exponents."""
    mag = torch.abs(base)
    phase = torch.angle(base)
    return (mag ** exp) * torch.exp(1j * exp * phase)


def safe_divide(num: torch.Tensor, den: torch.Tensor, eps: float = 1e-15) -> torch.Tensor:
    """Safe division avoiding divide by zero."""
    return num / (den + eps * torch.sign(den.real + 1e-30))


# =============================================================================
# CATEGORY 1: DEBYE FAMILY (RC Parallel Networks)
# =============================================================================
# Circuit: Each Debye term = R || C (RC parallel)
# Multiple Debye = Foster network (sum of RC parallel stages)

def debye(omega: torch.Tensor, tau: torch.Tensor) -> torch.Tensor:
    """
    Debye Relaxation
    Formula: 1 / (1 + j*omega*tau)

    Circuit: R || C where tau = R*C
    """
    return 1.0 / (1.0 + 1j * omega * tau)


def cole_cole(omega: torch.Tensor, tau: torch.Tensor, alpha: torch.Tensor) -> torch.Tensor:
    """
    Cole-Cole Relaxation
    Formula: 1 / (1 + (j*omega*tau)^alpha)

    Circuit: RC ladder (Foster network with distributed tau)
    Approximation: 5-7 stage RC ladder using Charef method
    Parameters: 0 < alpha <= 1 (alpha=1 -> Debye)
    """
    alpha = torch.clamp(alpha, 0.1, 1.0)
    return 1.0 / (1.0 + safe_power(1j * omega * tau, alpha))


def cole_davidson(omega: torch.Tensor, tau: torch.Tensor, beta: torch.Tensor) -> torch.Tensor:
    """
    Cole-Davidson Relaxation
    Formula: 1 / (1 + j*omega*tau)^beta

    Circuit: RC ladder (asymmetric distribution)
    Parameters: 0 < beta <= 1
    """
    beta = torch.clamp(beta, 0.1, 1.0)
    return 1.0 / (1.0 + 1j * omega * tau) ** beta


def havriliak_negami(omega: torch.Tensor, tau: torch.Tensor,
                     alpha: torch.Tensor, beta: torch.Tensor) -> torch.Tensor:
    """
    Havriliak-Negami Relaxation
    Formula: 1 / (1 + (j*omega*tau)^alpha)^beta

    Circuit: RC ladder (generalized Foster network)
    Special cases: alpha=1,beta=1 -> Debye; beta=1 -> Cole-Cole
    """
    alpha = torch.clamp(alpha, 0.1, 1.0)
    beta = torch.clamp(beta, 0.1, 1.0)
    return 1.0 / safe_power(1.0 + safe_power(1j * omega * tau, alpha), beta)


def debye_two_site(omega: torch.Tensor, tau1: torch.Tensor, tau2: torch.Tensor,
                   f: torch.Tensor) -> torch.Tensor:
    """
    Two-Site Debye Model
    Formula: f / (1 + j*omega*tau1) + (1-f) / (1 + j*omega*tau2)

    Circuit: Two RC parallel stages in series
    """
    f = torch.clamp(f, 0.0, 1.0)
    return f / (1.0 + 1j * omega * tau1) + (1.0 - f) / (1.0 + 1j * omega * tau2)


# =============================================================================
# CATEGORY 2: CONSTANT PHASE ELEMENTS (RC Ladder)
# =============================================================================
# Circuit: RC ladder with geometrically distributed R and C values
# Valsa method: R_k = R_0 * a^k, C_k = C_0 * b^k

def cpe(omega: torch.Tensor, Q: torch.Tensor, n: torch.Tensor) -> torch.Tensor:
    """
    Constant Phase Element (CPE)
    Formula: 1 / (Q * (j*omega)^n)

    Circuit: RC ladder (Valsa approximation)
    - n=1: capacitor (C = Q)
    - n=0.5: Warburg (RC ladder with R_k/C_k = const)
    - n=0: resistor (R = 1/Q)
    """
    n = torch.clamp(n, 0.01, 0.99)
    return 1.0 / (Q * safe_power(1j * omega, n))


def cpe_bounded(omega: torch.Tensor, R_0: torch.Tensor, tau: torch.Tensor,
                n: torch.Tensor) -> torch.Tensor:
    """
    Bounded CPE (finite DC resistance)
    Formula: R_0 / (1 + (j*omega*tau)^n)

    Circuit: R_0 in parallel with CPE (finite at DC)
    """
    n = torch.clamp(n, 0.1, 0.99)
    return R_0 / (1.0 + safe_power(1j * omega * tau, n))


# =============================================================================
# CATEGORY 3: DIFFUSION (RC Ladder - Warburg Type)
# =============================================================================
# Circuit: Semi-infinite = infinite RC ladder
# Finite = truncated RC ladder with termination

def warburg_infinite(omega: torch.Tensor, Aw: torch.Tensor) -> torch.Tensor:
    """
    Warburg Impedance (Semi-infinite diffusion)
    Formula: Aw / sqrt(j*omega) = Aw * (1-j) / sqrt(2*omega)

    Circuit: Infinite RC ladder with R_k = C_k (Warburg ladder)
    Approximation: 5-10 stage RC ladder
    """
    return Aw / torch.sqrt(1j * omega + 1e-15)


def warburg_finite(omega: torch.Tensor, R_d: torch.Tensor,
                   tau_d: torch.Tensor) -> torch.Tensor:
    """
    Finite-Length Warburg (tanh form)
    Formula: R_d * tanh(sqrt(j*omega*tau_d)) / sqrt(j*omega*tau_d)

    Circuit: RC ladder with open termination (blocking electrode)
    """
    s = torch.sqrt(1j * omega * tau_d + 1e-15)
    return R_d * torch.tanh(s) / s


def gerischer(omega: torch.Tensor, R_g: torch.Tensor, tau_g: torch.Tensor) -> torch.Tensor:
    """
    Gerischer Element
    Formula: R_g / sqrt(1 + j*omega*tau_g)

    Circuit: RC ladder with parallel R termination (reaction)
    """
    return R_g / torch.sqrt(1.0 + 1j * omega * tau_g)


# =============================================================================
# CATEGORY 4: TRANSMISSION LINE (Distributed RC)
# =============================================================================
# Circuit: Cascaded T or Pi sections of R and C

def transmission_line_open(omega: torch.Tensor, R: torch.Tensor,
                           C: torch.Tensor, length: torch.Tensor) -> torch.Tensor:
    """
    RC Transmission Line (Open termination)
    Formula: sqrt(R/(j*omega*C)) * coth(length * sqrt(R*j*omega*C))

    Circuit: N-section RC ladder with open end
    """
    gamma = torch.sqrt(R * 1j * omega * C + 1e-15)
    Z0 = torch.sqrt(R / (1j * omega * C + 1e-15))
    return Z0 / torch.tanh(gamma * length + 1e-10)


def transmission_line_short(omega: torch.Tensor, R: torch.Tensor,
                            C: torch.Tensor, length: torch.Tensor) -> torch.Tensor:
    """
    RC Transmission Line (Short termination)
    Formula: sqrt(R/(j*omega*C)) * tanh(length * sqrt(R*j*omega*C))

    Circuit: N-section RC ladder with short end
    """
    gamma = torch.sqrt(R * 1j * omega * C + 1e-15)
    Z0 = torch.sqrt(R / (1j * omega * C + 1e-15))
    return Z0 * torch.tanh(gamma * length)


# =============================================================================
# CATEGORY 5: SKIN EFFECT (RL Ladder - Dowell)
# =============================================================================
# Circuit: RL ladder with Dowell coefficients [3, 5, 7, 9, 11, ...]
# Each stage: R_k in series with L_k to ground

def skin_effect_dowell(omega: torch.Tensor, R_dc: torch.Tensor,
                       delta_ref: torch.Tensor) -> torch.Tensor:
    """
    Skin Effect - Dowell Formula
    Formula: R_dc * z * coth(z) where z = (1+j) * d/(2*delta)

    Circuit: RL ladder with Dowell coefficients
    R_k = R_dc/a_k, L_k = R_dc*tau/a_k where a_k = [3, 5, 7, 9, 11]
    """
    tau = delta_ref**2 / 2
    z = (1 + 1j) * torch.sqrt(omega * tau + 1e-15)
    return R_dc * z / torch.tanh(z + 1e-10)


def skin_effect_cylindrical(omega: torch.Tensor, R_dc: torch.Tensor,
                            radius: torch.Tensor, delta: torch.Tensor) -> torch.Tensor:
    """
    Skin Effect - Cylindrical Conductor (approximation)
    Formula: R_dc * (kr/2) * I0(kr)/I1(kr)

    Circuit: RL ladder (modified Dowell for round wire)
    """
    kr = (1 + 1j) * radius / delta
    ratio = 1.0 + 1.0 / (kr + 1e-10)
    return R_dc * (kr / 2.0) * ratio


def multilayer_winding(omega: torch.Tensor, R_dc: torch.Tensor,
                       delta_ref: torch.Tensor, n_layers: torch.Tensor) -> torch.Tensor:
    """
    Multi-Layer Winding (Dowell extended)
    Formula: R_dc * z * coth(z) * (1 + (2/3)*(n^2-1)*(z*tanh(z)-1))

    Circuit: RL ladder with proximity effect correction
    """
    tau = delta_ref**2 / 2
    z = (1 + 1j) * torch.sqrt(omega * tau + 1e-15)
    skin = z / torch.tanh(z + 1e-10)
    multilayer = (2.0/3.0) * (n_layers**2 - 1.0) * (z * torch.tanh(z) - 1.0)
    return R_dc * skin * (1.0 + multilayer)


# =============================================================================
# CATEGORY 6: MAGNETIC MATERIALS (RC/RL Networks)
# =============================================================================
# Magnetic relaxation maps to electrical via:
# mu(omega) -> L(omega) or Z_magnetic(omega)
#
# Physical mechanisms in magnetic materials:
# 1. Domain wall motion: Low frequency, high permeability
# 2. Spin rotation: Higher frequency, lower contribution
# 3. Ferromagnetic resonance (FMR): GHz range for thin films/nanoparticles
# 4. Eddy current losses: Depends on conductivity and geometry

def magnetic_debye(omega: torch.Tensor, mu_s: torch.Tensor, mu_inf: torch.Tensor,
                   tau: torch.Tensor) -> torch.Tensor:
    """
    Magnetic Debye Relaxation
    Formula: mu_inf + (mu_s - mu_inf) / (1 + j*omega*tau)

    Circuit equivalent: L_inf in series with (L_s-L_inf) || R
    where tau = (L_s - L_inf) / R
    """
    return mu_inf + (mu_s - mu_inf) / (1.0 + 1j * omega * tau)


def magnetic_cole_cole(omega: torch.Tensor, mu_s: torch.Tensor, mu_inf: torch.Tensor,
                       tau: torch.Tensor, alpha: torch.Tensor) -> torch.Tensor:
    """
    Magnetic Cole-Cole Relaxation
    Formula: mu_inf + (mu_s - mu_inf) / (1 + (j*omega*tau)^alpha)

    Circuit: RL ladder (Foster network for permeability)
    """
    alpha = torch.clamp(alpha, 0.1, 1.0)
    return mu_inf + (mu_s - mu_inf) / (1.0 + safe_power(1j * omega * tau, alpha))


def two_relaxation_permeability(omega: torch.Tensor, mu_dc: torch.Tensor,
                                 omega_1: torch.Tensor, omega_2: torch.Tensor) -> torch.Tensor:
    """
    Two-Relaxation Permeability Model
    Formula: mu_dc / ((1 + j*omega/omega_1) * (1 + j*omega/omega_2))

    Circuit: Two RL stages in series (domain wall + spin rotation)
    """
    return mu_dc / ((1.0 + 1j * omega / omega_1) * (1.0 + 1j * omega / omega_2))


def domain_wall_relaxation(omega: torch.Tensor, mu_dw: torch.Tensor,
                            omega_dw: torch.Tensor, beta: torch.Tensor) -> torch.Tensor:
    """
    Domain Wall Relaxation (Cole-Cole type)
    Formula: mu_dw / (1 + (j*omega/omega_dw)^beta)

    Physical basis: Domain wall motion in polycrystalline ferrites.
    Domain walls move in response to applied field, but are pinned by
    grain boundaries and defects, leading to distribution of relaxation times.

    Circuit: RL ladder (distributed time constants)

    Parameters:
        mu_dw: Domain wall contribution to permeability
        omega_dw: Characteristic angular frequency of domain wall motion
        beta: Distribution parameter (0.5-1.0, beta=1 -> single relaxation)

    Applications:
        - MnZn ferrites (power frequencies)
        - NiZn ferrites (RF frequencies)
        - Soft magnetic composites
    """
    beta = torch.clamp(beta, 0.3, 1.0)
    return mu_dw / (1.0 + safe_power(1j * omega / omega_dw, beta))


def spin_rotation_relaxation(omega: torch.Tensor, mu_spin: torch.Tensor,
                              omega_spin: torch.Tensor) -> torch.Tensor:
    """
    Spin Rotation Relaxation (Debye type)
    Formula: mu_spin / (1 + j*omega/omega_spin)

    Physical basis: Rotation of magnetic moments within domains.
    Higher frequency process than domain wall motion.
    Damping from spin-orbit coupling and magnon-phonon interactions.

    Circuit: Single RL parallel (L_spin || R_spin)

    Parameters:
        mu_spin: Spin rotation contribution to permeability
        omega_spin: Characteristic angular frequency of spin rotation

    Applications:
        - High-frequency ferrites
        - Thin magnetic films
        - RF absorbers
    """
    return mu_spin / (1.0 + 1j * omega / omega_spin)


def ferromagnetic_resonance(omega: torch.Tensor, mu_0_r: torch.Tensor,
                            omega_0: torch.Tensor, alpha_G: torch.Tensor) -> torch.Tensor:
    """
    Ferromagnetic Resonance (FMR) - Landau-Lifshitz-Gilbert form
    Formula: mu_0_r * omega_0^2 / (omega_0^2 - omega^2 + j*alpha_G*omega*omega_0)

    Physical basis: Precession of magnetization around effective field.
    Based on Landau-Lifshitz-Gilbert equation:
        dM/dt = -gamma * M x H_eff + (alpha/M_s) * M x dM/dt

    Resonance occurs at omega = omega_0 = gamma * H_eff (Kittel formula)

    Circuit: RLC parallel resonator
        L = 1/(mu_0_r * omega_0^2)
        C = mu_0_r
        R = 1/(alpha_G * omega_0 * mu_0_r)

    Parameters:
        mu_0_r: DC relative permeability (susceptibility contribution)
        omega_0: Resonance angular frequency (= gamma * B_eff)
        alpha_G: Gilbert damping parameter (typically 0.001 - 0.1)

    Applications:
        - Thin magnetic films (GHz range)
        - Magnetic nanoparticles
        - Ferrite resonators
        - Microwave absorbers

    Reference:
        Kittel, "On the Theory of Ferromagnetic Resonance Absorption", Phys. Rev. 73, 155 (1948)
    """
    alpha_G = torch.clamp(alpha_G, 0.001, 0.5)
    denominator = omega_0**2 - omega**2 + 1j * alpha_G * omega * omega_0
    return mu_0_r * omega_0**2 / (denominator + 1e-15)


def domain_wall_resonance(omega: torch.Tensor, chi_dw: torch.Tensor,
                           omega_dw: torch.Tensor, gamma_dw: torch.Tensor) -> torch.Tensor:
    """
    Domain Wall Resonance (oscillatory domain wall motion)
    Formula: chi_dw * omega_dw^2 / (omega_dw^2 - omega^2 + j*gamma_dw*omega)

    Physical basis: Domain walls can resonate like mechanical oscillators.
    The wall has effective mass (from eddy currents and spin inertia) and
    restoring force (from magnetostatic energy and pinning).

    Typically lower frequency than FMR (kHz to MHz for bulk materials).

    Circuit: RLC series resonator
        L = 1/chi_dw (effective mass)
        C = chi_dw/omega_dw^2 (stiffness)
        R = gamma_dw/chi_dw (damping)

    Parameters:
        chi_dw: Domain wall susceptibility contribution
        omega_dw: Domain wall resonance frequency
        gamma_dw: Damping coefficient (from eddy currents)

    Applications:
        - Power ferrites with domain wall oscillations
        - Magnetic after-effect studies
        - Low-frequency absorbers
    """
    gamma_dw = torch.clamp(gamma_dw, 1e-3, 10.0)
    denominator = omega_dw**2 - omega**2 + 1j * gamma_dw * omega
    return chi_dw * omega_dw**2 / (denominator + 1e-15)


def snoek_limit_model(omega: torch.Tensor, mu_s: torch.Tensor,
                      f_cutoff: torch.Tensor) -> torch.Tensor:
    """
    Snoek's Limit Permeability Model
    Formula: mu_s / (1 + j*omega/(2*pi*f_cutoff)) with constraint (mu_s-1)*f_cutoff = const

    Physical basis: Snoek's law states that for polycrystalline ferrites:
        (mu_s - 1) * f_r = (2/3) * gamma * M_s / (2*pi)

    This fundamental limit arises from the gyromagnetic ratio and saturation
    magnetization, independent of material details.

    Circuit: Single RL parallel with L proportional to mu_s

    Parameters:
        mu_s: Static permeability (limited by Snoek's law)
        f_cutoff: Cutoff frequency (Hz)

    Applications:
        - Ferrite material selection
        - High-frequency inductor design
        - Estimating permeability rolloff

    Reference:
        Snoek, "Dispersion and absorption in magnetic ferrites at frequencies above one Mc/s",
        Physica 14, 207 (1948)
    """
    omega_cutoff = 2 * np.pi * f_cutoff
    return 1.0 + (mu_s - 1.0) / (1.0 + 1j * omega / omega_cutoff)


# =============================================================================
# CATEGORY 7: RLC RESONANCE (Direct Circuit Elements)
# =============================================================================
# These are direct circuit representations

def rlc_series(omega: torch.Tensor, R: torch.Tensor, L: torch.Tensor,
               C: torch.Tensor) -> torch.Tensor:
    """
    RLC Series Resonance
    Formula: R + j*omega*L + 1/(j*omega*C)

    Circuit: R, L, C in series
    """
    return R + 1j * omega * L + 1.0 / (1j * omega * C + 1e-15)


def rlc_parallel(omega: torch.Tensor, R: torch.Tensor, L: torch.Tensor,
                 C: torch.Tensor) -> torch.Tensor:
    """
    RLC Parallel Resonance
    Formula: 1 / (1/R + j*omega*C + 1/(j*omega*L))

    Circuit: R, L, C in parallel
    """
    Y = 1.0/R + 1j * omega * C + 1.0 / (1j * omega * L + 1e-15)
    return 1.0 / (Y + 1e-15)


def piezoelectric_bvd(omega: torch.Tensor, C_0: torch.Tensor, R_m: torch.Tensor,
                      L_m: torch.Tensor, C_m: torch.Tensor) -> torch.Tensor:
    """
    Piezoelectric (Butterworth-Van Dyke) Model
    Formula: 1/(j*omega*C_0) || (R_m + j*omega*L_m + 1/(j*omega*C_m))

    Circuit: Static capacitor C_0 in parallel with motional RLC
    """
    Z_motional = R_m + 1j * omega * L_m + 1.0 / (1j * omega * C_m + 1e-15)
    Z_static = 1.0 / (1j * omega * C_0 + 1e-15)
    return 1.0 / (1.0/Z_motional + 1.0/Z_static)


# =============================================================================
# CATEGORY 8: VISCOELASTIC (Mechanical-Electrical Analog)
# =============================================================================
# Mechanical elements map to electrical:
# Spring (stiffness K) -> Capacitor (C = 1/K)
# Dashpot (viscosity eta) -> Resistor (R = eta)
# Mass (M) -> Inductor (L = M)

def maxwell_element(omega: torch.Tensor, R: torch.Tensor, C: torch.Tensor) -> torch.Tensor:
    """
    Maxwell Model (Electrical Analog)
    Mechanical: spring-dashpot in series
    Electrical: R in series with C

    Formula: R + 1/(j*omega*C) ... but for impedance form:
    Z = R * 1/(j*omega*C) / (R + 1/(j*omega*C)) = R/(1 + j*omega*R*C)
    """
    return R / (1.0 + 1j * omega * R * C)


def voigt_element(omega: torch.Tensor, R: torch.Tensor, C: torch.Tensor) -> torch.Tensor:
    """
    Voigt/Kelvin Model (Electrical Analog)
    Mechanical: spring-dashpot in parallel
    Electrical: R in parallel with C

    Formula: R || C = R / (1 + j*omega*R*C)
    """
    return R / (1.0 + 1j * omega * R * C)


def standard_linear_solid(omega: torch.Tensor, R_1: torch.Tensor, R_2: torch.Tensor,
                          C: torch.Tensor) -> torch.Tensor:
    """
    Standard Linear Solid (Zener Model - Electrical Analog)
    Mechanical: Maxwell element in parallel with spring
    Electrical: R_1 in parallel with (R_2 series C)

    Formula: R_1 || (R_2 + 1/(j*omega*C))
    """
    Z_rc = R_2 + 1.0 / (1j * omega * C + 1e-15)
    return 1.0 / (1.0/R_1 + 1.0/Z_rc)


# =============================================================================
# BASIS FUNCTION REGISTRY (Circuit-Compatible Only)
# =============================================================================

BASIS_FUNCTIONS: Dict[str, Dict] = {
    # CATEGORY 1: Debye Family (RC Parallel)
    'debye': {
        'func': debye,
        'params': ['tau'],
        'category': 'Debye family',
        'circuit': 'R || C (tau = RC)'
    },
    'cole_cole': {
        'func': cole_cole,
        'params': ['tau', 'alpha'],
        'category': 'Debye family',
        'circuit': 'RC ladder (5-7 stages, Charef method)'
    },
    'cole_davidson': {
        'func': cole_davidson,
        'params': ['tau', 'beta'],
        'category': 'Debye family',
        'circuit': 'RC ladder (asymmetric distribution)'
    },
    'havriliak_negami': {
        'func': havriliak_negami,
        'params': ['tau', 'alpha', 'beta'],
        'category': 'Debye family',
        'circuit': 'RC ladder (generalized Foster)'
    },
    'debye_two_site': {
        'func': debye_two_site,
        'params': ['tau1', 'tau2', 'f'],
        'category': 'Debye family',
        'circuit': 'Two RC parallel stages'
    },

    # CATEGORY 2: CPE (RC Ladder)
    'cpe': {
        'func': cpe,
        'params': ['Q', 'n'],
        'category': 'CPE',
        'circuit': 'RC ladder (Valsa approximation, 5-10 stages)'
    },
    'cpe_bounded': {
        'func': cpe_bounded,
        'params': ['R_0', 'tau', 'n'],
        'category': 'CPE',
        'circuit': 'R_0 || CPE ladder'
    },

    # CATEGORY 3: Diffusion (RC Ladder)
    'warburg_infinite': {
        'func': warburg_infinite,
        'params': ['Aw'],
        'category': 'Diffusion',
        'circuit': 'RC ladder (n=0.5 CPE, 5-10 stages)'
    },
    'warburg_finite': {
        'func': warburg_finite,
        'params': ['R_d', 'tau_d'],
        'category': 'Diffusion',
        'circuit': 'RC ladder with open termination'
    },
    'gerischer': {
        'func': gerischer,
        'params': ['R_g', 'tau_g'],
        'category': 'Diffusion',
        'circuit': 'RC ladder with R termination'
    },

    # CATEGORY 4: Transmission Line
    'tl_open': {
        'func': transmission_line_open,
        'params': ['R', 'C', 'length'],
        'category': 'Transmission line',
        'circuit': 'N-section RC ladder (open end)'
    },
    'tl_short': {
        'func': transmission_line_short,
        'params': ['R', 'C', 'length'],
        'category': 'Transmission line',
        'circuit': 'N-section RC ladder (shorted end)'
    },

    # CATEGORY 5: Skin Effect (RL Ladder)
    'skin_dowell': {
        'func': skin_effect_dowell,
        'params': ['R_dc', 'delta_ref'],
        'category': 'Skin effect',
        'circuit': 'RL ladder (Dowell: a_k = 3,5,7,9,11)'
    },
    'skin_cylindrical': {
        'func': skin_effect_cylindrical,
        'params': ['R_dc', 'radius', 'delta'],
        'category': 'Skin effect',
        'circuit': 'RL ladder (modified for round wire)'
    },
    'multilayer_winding': {
        'func': multilayer_winding,
        'params': ['R_dc', 'delta_ref', 'n_layers'],
        'category': 'Skin effect',
        'circuit': 'RL ladder with proximity correction'
    },

    # CATEGORY 6: Magnetic (RL Network)
    'magnetic_debye': {
        'func': magnetic_debye,
        'params': ['mu_s', 'mu_inf', 'tau'],
        'category': 'Magnetic',
        'circuit': 'L_inf + (L_delta || R)'
    },
    'magnetic_cole_cole': {
        'func': magnetic_cole_cole,
        'params': ['mu_s', 'mu_inf', 'tau', 'alpha'],
        'category': 'Magnetic',
        'circuit': 'RL ladder (Foster for permeability)'
    },
    'two_relaxation_mu': {
        'func': two_relaxation_permeability,
        'params': ['mu_dc', 'omega_1', 'omega_2'],
        'category': 'Magnetic',
        'circuit': 'Two RL stages (domain wall + spin)'
    },
    'domain_wall_relaxation': {
        'func': domain_wall_relaxation,
        'params': ['mu_dw', 'omega_dw', 'beta'],
        'category': 'Magnetic',
        'circuit': 'RL ladder (distributed domain wall motion)'
    },
    'spin_rotation': {
        'func': spin_rotation_relaxation,
        'params': ['mu_spin', 'omega_spin'],
        'category': 'Magnetic',
        'circuit': 'L_spin || R_spin (single Debye)'
    },
    'fmr': {
        'func': ferromagnetic_resonance,
        'params': ['mu_0_r', 'omega_0', 'alpha_G'],
        'category': 'Magnetic',
        'circuit': 'RLC parallel (Landau-Lifshitz-Gilbert)'
    },
    'domain_wall_resonance': {
        'func': domain_wall_resonance,
        'params': ['chi_dw', 'omega_dw', 'gamma_dw'],
        'category': 'Magnetic',
        'circuit': 'RLC series (oscillatory domain wall)'
    },
    'snoek_limit': {
        'func': snoek_limit_model,
        'params': ['mu_s', 'f_cutoff'],
        'category': 'Magnetic',
        'circuit': 'RL parallel (Snoek-limited)'
    },

    # CATEGORY 7: RLC Resonance
    'rlc_series': {
        'func': rlc_series,
        'params': ['R', 'L', 'C'],
        'category': 'Resonance',
        'circuit': 'R-L-C series'
    },
    'rlc_parallel': {
        'func': rlc_parallel,
        'params': ['R', 'L', 'C'],
        'category': 'Resonance',
        'circuit': 'R || L || C'
    },
    'piezo_bvd': {
        'func': piezoelectric_bvd,
        'params': ['C_0', 'R_m', 'L_m', 'C_m'],
        'category': 'Resonance',
        'circuit': 'C_0 || (R_m + L_m + C_m)'
    },

    # CATEGORY 8: Viscoelastic (Mechanical-Electrical Analog)
    'maxwell': {
        'func': maxwell_element,
        'params': ['R', 'C'],
        'category': 'Viscoelastic',
        'circuit': 'R in series with C (parallel form)'
    },
    'voigt': {
        'func': voigt_element,
        'params': ['R', 'C'],
        'category': 'Viscoelastic',
        'circuit': 'R || C'
    },
    'sls': {
        'func': standard_linear_solid,
        'params': ['R_1', 'R_2', 'C'],
        'category': 'Viscoelastic',
        'circuit': 'R_1 || (R_2 + C)'
    },
}


def list_basis_functions():
    """Print available basis functions with circuit equivalents."""
    categories = {}
    for name, info in BASIS_FUNCTIONS.items():
        cat = info['category']
        if cat not in categories:
            categories[cat] = []
        categories[cat].append((name, info['params'], info['circuit']))

    print("=" * 80)
    print("Relaxation Basis Function Library (Circuit-Compatible)")
    print("=" * 80)
    print("All functions have direct circuit equivalents (RLC ladders)")
    print("=" * 80)

    total = 0
    for cat in ['Debye family', 'CPE', 'Diffusion', 'Transmission line',
                'Skin effect', 'Magnetic', 'Resonance', 'Viscoelastic']:
        if cat in categories:
            print(f"\n{cat} ({len(categories[cat])} functions):")
            for name, params, circuit in categories[cat]:
                params_str = ', '.join(params)
                print(f"  {name:20s} ({params_str})")
                print(f"      Circuit: {circuit}")
            total += len(categories[cat])

    print("\n" + "=" * 80)
    print(f"Total: {total} circuit-compatible basis functions")
    print("=" * 80)


if __name__ == '__main__':
    list_basis_functions()
