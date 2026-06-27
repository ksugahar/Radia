"""
URN Benchmark Suite for Paper

Comprehensive benchmark tests for KAN-inspired Universal Relaxation Network.
Collects data for publication-quality results.

Test Categories:
1. Ferrite permeability (various materials)
2. Conductor skin effect (flat, round, multilayer)
3. Electrochemical impedance (Randles, porous electrode)
4. Transformer/inductor windings
5. Dielectric relaxation (polymers, ceramics)
6. Comparison with Vector Fitting

Output: JSON files with detailed results for paper figures.
"""

import numpy as np
import torch
import json
import time
from datetime import datetime
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple
import os

# Import URN
from universal_relaxation_network import URNConfig, train_urn, UniversalRelaxationNetwork


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class BenchmarkResult:
    """Result from a single benchmark test."""
    test_name: str
    category: str
    true_model: str
    true_params: Dict

    # Fit quality
    max_rel_error: float
    mean_rel_error: float
    rms_error: float

    # Discovered model
    active_components: Dict
    n_active: int

    # Timing
    train_time_sec: float

    # Raw data for plotting
    freqs: List[float]
    Z_true_real: List[float]
    Z_true_imag: List[float]
    Z_pred_real: List[float]
    Z_pred_imag: List[float]


def save_results(results: List[BenchmarkResult], filename: str):
    """Save benchmark results to JSON."""
    data = {
        'timestamp': datetime.now().isoformat(),
        'n_tests': len(results),
        'results': [asdict(r) for r in results]
    }
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"Results saved to {filename}")


# =============================================================================
# CATEGORY 1: FERRITE PERMEABILITY
# =============================================================================

def generate_ferrite_data():
    """Generate ferrite permeability test cases."""
    test_cases = []

    # 1.1 MnZn Ferrite (power applications, ~100 kHz)
    test_cases.append({
        'name': 'MnZn_power_ferrite',
        'mu_s': 2500,
        'mu_inf': 1,
        'f0': 500e3,  # 500 kHz
        'alpha': 1.0,  # Pure Debye
        'description': 'MnZn ferrite for power electronics (100-500 kHz)'
    })

    # 1.2 NiZn Ferrite (EMI suppression, ~10 MHz)
    test_cases.append({
        'name': 'NiZn_EMI_ferrite',
        'mu_s': 500,
        'mu_inf': 1,
        'f0': 10e6,
        'alpha': 0.85,  # Cole-Cole
        'description': 'NiZn ferrite for EMI suppression (1-100 MHz)'
    })

    # 1.3 High-frequency ferrite (~100 MHz)
    test_cases.append({
        'name': 'HF_ferrite',
        'mu_s': 120,
        'mu_inf': 1,
        'f0': 100e6,
        'alpha': 0.9,
        'description': 'High-frequency ferrite for RF applications'
    })

    # 1.4 Lossy ferrite with complex permeability
    test_cases.append({
        'name': 'Lossy_ferrite',
        'mu_s': 1000,
        'mu_inf': 1,
        'f0': 2e6,
        'alpha': 0.7,  # Strong dispersion
        'description': 'Lossy ferrite with broad relaxation spectrum'
    })

    # 1.5 Two-relaxation ferrite (domain wall + spin)
    test_cases.append({
        'name': 'Two_relaxation_ferrite',
        'mu_dc': 3000,
        'f1': 100e3,  # Domain wall
        'f2': 5e6,    # Spin rotation
        'description': 'Ferrite with domain wall and spin relaxation'
    })

    return test_cases


def run_ferrite_benchmarks(verbose: bool = True) -> List[BenchmarkResult]:
    """Run ferrite permeability benchmarks."""
    results = []
    test_cases = generate_ferrite_data()

    print("\n" + "="*70)
    print("CATEGORY 1: FERRITE PERMEABILITY")
    print("="*70)

    for tc in test_cases:
        print(f"\nTest: {tc['name']}")
        print(f"  {tc['description']}")

        # Generate frequency sweep
        freqs = np.logspace(3, 9, 100)  # 1 kHz to 1 GHz
        omega = 2 * np.pi * freqs

        # Generate true data
        if 'alpha' in tc:
            # Cole-Cole model
            tau = 1 / (2 * np.pi * tc['f0'])
            mu_true = tc['mu_inf'] + (tc['mu_s'] - tc['mu_inf']) / (1 + (1j * omega * tau)**tc['alpha'])
            true_model = f"Cole-Cole (alpha={tc['alpha']})"
            true_params = {'mu_s': tc['mu_s'], 'mu_inf': tc['mu_inf'], 'f0': tc['f0'], 'alpha': tc['alpha']}
        else:
            # Two-relaxation model
            omega_1 = 2 * np.pi * tc['f1']
            omega_2 = 2 * np.pi * tc['f2']
            mu_true = tc['mu_dc'] / ((1 + 1j * omega / omega_1) * (1 + 1j * omega / omega_2))
            true_model = "Two-relaxation"
            true_params = {'mu_dc': tc['mu_dc'], 'f1': tc['f1'], 'f2': tc['f2']}

        # Configure URN
        config = URNConfig(
            n_debye=3,
            n_cole_cole=2,
            n_cole_davidson=1,
            n_havriliak_negami=1,
            n_cpe=1,
            n_warburg=0,
            n_gerischer=0,
            n_rlc=0,
            n_skin_effect=0,
            sparsity_weight=0.005,
            n_epochs=4000,
            n_restarts=3
        )

        # Train
        start_time = time.time()
        model = train_urn(freqs, mu_true, config, verbose=False)
        train_time = time.time() - start_time

        # Evaluate
        omega_torch = torch.tensor(omega, dtype=torch.float64)
        with torch.no_grad():
            mu_pred = model(omega_torch).numpy()

        rel_err = np.abs(mu_pred - mu_true) / np.abs(mu_true)

        # Get active components
        active = model.get_active_components()
        n_active = sum(len(v) for v in active.values())

        # Convert active components to serializable format
        active_serializable = {}
        for key, comps in active.items():
            active_serializable[key] = []
            for c in comps:
                comp_dict = {}
                for k, v in c.items():
                    if hasattr(v, 'item'):
                        comp_dict[k] = float(v.item())
                    elif isinstance(v, (np.floating, np.integer)):
                        comp_dict[k] = float(v)
                    else:
                        comp_dict[k] = v
                active_serializable[key].append(comp_dict)

        result = BenchmarkResult(
            test_name=tc['name'],
            category='ferrite',
            true_model=true_model,
            true_params=true_params,
            max_rel_error=float(rel_err.max()),
            mean_rel_error=float(rel_err.mean()),
            rms_error=float(np.sqrt(np.mean(rel_err**2))),
            active_components=active_serializable,
            n_active=n_active,
            train_time_sec=train_time,
            freqs=freqs.tolist(),
            Z_true_real=mu_true.real.tolist(),
            Z_true_imag=mu_true.imag.tolist(),
            Z_pred_real=mu_pred.real.tolist(),
            Z_pred_imag=mu_pred.imag.tolist()
        )
        results.append(result)

        print(f"  Max error: {rel_err.max()*100:.2f}%, Mean: {rel_err.mean()*100:.2f}%")
        print(f"  Active components: {n_active} ({list(active.keys())})")
        print(f"  Train time: {train_time:.1f} s")

    return results


# =============================================================================
# CATEGORY 2: CONDUCTOR SKIN EFFECT
# =============================================================================

def generate_skin_effect_data():
    """Generate skin effect test cases."""
    test_cases = []

    # 2.1 Flat conductor (transformer foil winding)
    test_cases.append({
        'name': 'flat_foil_winding',
        'R_dc': 0.01,  # 10 mOhm
        'thickness': 0.1e-3,  # 0.1 mm foil
        'sigma': 5.8e7,  # Copper
        'geometry': 'flat',
        'description': 'Copper foil winding (0.1 mm thick)'
    })

    # 2.2 Round wire
    test_cases.append({
        'name': 'round_wire_AWG20',
        'R_dc': 0.033,  # 33 mOhm/m for AWG20
        'radius': 0.4e-3,  # 0.4 mm radius
        'sigma': 5.8e7,
        'geometry': 'cylindrical',
        'description': 'AWG20 copper wire (0.8 mm diameter)'
    })

    # 2.3 Multilayer winding (5 layers)
    test_cases.append({
        'name': 'multilayer_5layers',
        'R_dc': 0.01,
        'thickness': 0.2e-3,
        'n_layers': 5,
        'sigma': 5.8e7,
        'geometry': 'multilayer',
        'description': '5-layer transformer winding with proximity effect'
    })

    # 2.4 Aluminum conductor
    test_cases.append({
        'name': 'aluminum_busbar',
        'R_dc': 0.005,  # 5 mOhm
        'thickness': 2e-3,  # 2 mm
        'sigma': 3.5e7,  # Aluminum
        'geometry': 'flat',
        'description': 'Aluminum busbar (2 mm thick)'
    })

    # 2.5 Litz wire approximation
    test_cases.append({
        'name': 'litz_wire',
        'R_dc': 0.02,
        'strand_radius': 0.05e-3,  # 50 um strand
        'n_strands': 100,
        'sigma': 5.8e7,
        'geometry': 'litz',
        'description': 'Litz wire (100 x 50um strands)'
    })

    return test_cases


def compute_skin_effect_impedance(tc: Dict, freqs: np.ndarray) -> np.ndarray:
    """Compute skin effect impedance for various geometries."""
    omega = 2 * np.pi * freqs
    mu_0 = 4 * np.pi * 1e-7

    if tc['geometry'] == 'flat':
        # Flat conductor: z*coth(z)
        delta = np.sqrt(2 / (omega * mu_0 * tc['sigma'] + 1e-20))
        d = tc['thickness']
        z = (1 + 1j) * d / (2 * delta)
        Z = tc['R_dc'] * z / np.tanh(z + 1e-10)

    elif tc['geometry'] == 'cylindrical':
        # Cylindrical conductor: simplified approximation
        delta = np.sqrt(2 / (omega * mu_0 * tc['sigma'] + 1e-20))
        r = tc['radius']
        kr = (1 + 1j) * r / delta
        # Approximate I0/I1 ratio
        Z = tc['R_dc'] * (kr / 2) * (1 + 1 / (kr + 0.1))

    elif tc['geometry'] == 'multilayer':
        # Dowell model with proximity effect
        delta = np.sqrt(2 / (omega * mu_0 * tc['sigma'] + 1e-20))
        d = tc['thickness']
        n = tc['n_layers']
        z = (1 + 1j) * d / (2 * delta)
        skin = z / np.tanh(z + 1e-10)
        proximity = (2/3) * (n**2 - 1) * (z * np.tanh(z) - 1)
        Z = tc['R_dc'] * skin * (1 + proximity)

    elif tc['geometry'] == 'litz':
        # Litz wire: skin effect in each strand
        delta = np.sqrt(2 / (omega * mu_0 * tc['sigma'] + 1e-20))
        r = tc['strand_radius']
        z = (1 + 1j) * r / delta
        # Individual strand impedance
        Z_strand = z / np.tanh(z + 1e-10)
        Z = tc['R_dc'] * Z_strand  # Parallel strands already in R_dc

    return Z


def run_skin_effect_benchmarks(verbose: bool = True) -> List[BenchmarkResult]:
    """Run skin effect benchmarks."""
    results = []
    test_cases = generate_skin_effect_data()

    print("\n" + "="*70)
    print("CATEGORY 2: CONDUCTOR SKIN EFFECT")
    print("="*70)

    for tc in test_cases:
        print(f"\nTest: {tc['name']}")
        print(f"  {tc['description']}")

        # Generate frequency sweep
        freqs = np.logspace(1, 7, 100)  # 10 Hz to 10 MHz

        # Generate true data
        Z_true = compute_skin_effect_impedance(tc, freqs)

        # Configure URN (emphasize skin effect)
        config = URNConfig(
            n_debye=2,
            n_cole_cole=1,
            n_cole_davidson=0,
            n_havriliak_negami=0,
            n_cpe=2,
            n_warburg=1,
            n_gerischer=0,
            n_rlc=0,
            n_skin_effect=3,
            sparsity_weight=0.01,
            n_epochs=5000,
            n_restarts=3
        )

        # Train
        start_time = time.time()
        model = train_urn(freqs, Z_true, config, verbose=False)
        train_time = time.time() - start_time

        # Evaluate
        omega = 2 * np.pi * freqs
        omega_torch = torch.tensor(omega, dtype=torch.float64)
        with torch.no_grad():
            Z_pred = model(omega_torch).numpy()

        rel_err = np.abs(Z_pred - Z_true) / np.abs(Z_true)

        # Get active components
        active = model.get_active_components()
        n_active = sum(len(v) for v in active.values())

        # Convert to serializable
        active_serializable = {}
        for key, comps in active.items():
            active_serializable[key] = []
            for c in comps:
                comp_dict = {}
                for k, v in c.items():
                    if hasattr(v, 'item'):
                        comp_dict[k] = float(v.item())
                    elif isinstance(v, (np.floating, np.integer)):
                        comp_dict[k] = float(v)
                    else:
                        comp_dict[k] = v
                active_serializable[key].append(comp_dict)

        result = BenchmarkResult(
            test_name=tc['name'],
            category='skin_effect',
            true_model=tc['geometry'],
            true_params={k: v for k, v in tc.items() if k not in ['name', 'description', 'geometry']},
            max_rel_error=float(rel_err.max()),
            mean_rel_error=float(rel_err.mean()),
            rms_error=float(np.sqrt(np.mean(rel_err**2))),
            active_components=active_serializable,
            n_active=n_active,
            train_time_sec=train_time,
            freqs=freqs.tolist(),
            Z_true_real=Z_true.real.tolist(),
            Z_true_imag=Z_true.imag.tolist(),
            Z_pred_real=Z_pred.real.tolist(),
            Z_pred_imag=Z_pred.imag.tolist()
        )
        results.append(result)

        # Check if skin effect was identified
        skin_detected = 'skin_effect' in active

        print(f"  Max error: {rel_err.max()*100:.2f}%, Mean: {rel_err.mean()*100:.2f}%")
        print(f"  Active: {n_active} ({list(active.keys())})")
        print(f"  Skin effect detected: {'YES' if skin_detected else 'NO'}")
        print(f"  Train time: {train_time:.1f} s")

    return results


# =============================================================================
# CATEGORY 3: ELECTROCHEMICAL IMPEDANCE
# =============================================================================

def generate_eis_data():
    """Generate electrochemical impedance test cases."""
    test_cases = []

    # 3.1 Simple Randles circuit
    test_cases.append({
        'name': 'randles_simple',
        'Rs': 10,      # Solution resistance
        'Cdl': 1e-6,   # Double layer capacitance
        'Rct': 100,    # Charge transfer resistance
        'Aw': 50,      # Warburg coefficient
        'model': 'randles',
        'description': 'Simple Randles circuit: Rs + (Cdl || (Rct + W))'
    })

    # 3.2 Randles with CPE (non-ideal capacitance)
    test_cases.append({
        'name': 'randles_cpe',
        'Rs': 5,
        'Q': 1e-5,     # CPE coefficient
        'n': 0.85,     # CPE exponent
        'Rct': 200,
        'Aw': 30,
        'model': 'randles_cpe',
        'description': 'Randles with CPE: Rs + (CPE || (Rct + W))'
    })

    # 3.3 Porous electrode (transmission line)
    test_cases.append({
        'name': 'porous_electrode',
        'Rs': 1,
        'R_pore': 50,   # Pore resistance
        'C_dl': 1e-4,   # Distributed capacitance
        'L': 1e-3,      # Pore length
        'model': 'porous',
        'description': 'Porous electrode (de Levie model)'
    })

    # 3.4 Li-ion battery (two RC + Warburg)
    test_cases.append({
        'name': 'liion_battery',
        'Rs': 0.05,
        'R1': 0.02,    # SEI resistance
        'C1': 1e-3,    # SEI capacitance
        'R2': 0.1,     # Charge transfer
        'C2': 0.1,     # Double layer
        'Aw': 0.5,
        'model': 'battery',
        'description': 'Li-ion battery: Rs + (R1||C1) + (R2||C2) + W'
    })

    # 3.5 Fuel cell (finite Warburg)
    test_cases.append({
        'name': 'fuel_cell',
        'Rs': 0.1,
        'Rct': 0.5,
        'Cdl': 0.01,
        'R_d': 0.3,    # Diffusion resistance
        'tau_d': 1.0,  # Diffusion time constant
        'model': 'fuel_cell',
        'description': 'PEMFC: Rs + (Cdl || (Rct + W_finite))'
    })

    return test_cases


def compute_eis_impedance(tc: Dict, freqs: np.ndarray) -> np.ndarray:
    """Compute EIS impedance for various models."""
    omega = 2 * np.pi * freqs

    if tc['model'] == 'randles':
        Z_warburg = tc['Aw'] / np.sqrt(1j * omega + 1e-15)
        Z_faradaic = tc['Rct'] + Z_warburg
        Y_parallel = 1j * omega * tc['Cdl'] + 1 / Z_faradaic
        Z = tc['Rs'] + 1 / Y_parallel

    elif tc['model'] == 'randles_cpe':
        Z_warburg = tc['Aw'] / np.sqrt(1j * omega + 1e-15)
        Z_faradaic = tc['Rct'] + Z_warburg
        Y_cpe = tc['Q'] * (1j * omega) ** tc['n']
        Y_parallel = Y_cpe + 1 / Z_faradaic
        Z = tc['Rs'] + 1 / Y_parallel

    elif tc['model'] == 'porous':
        # de Levie transmission line
        gamma = np.sqrt(tc['R_pore'] * 1j * omega * tc['C_dl'] + 1e-15)
        Z0 = np.sqrt(tc['R_pore'] / (1j * omega * tc['C_dl'] + 1e-15))
        Z_pore = Z0 / np.tanh(gamma * tc['L'] + 1e-10)
        Z = tc['Rs'] + Z_pore

    elif tc['model'] == 'battery':
        Z_sei = tc['R1'] / (1 + 1j * omega * tc['R1'] * tc['C1'])
        Z_ct = tc['R2'] / (1 + 1j * omega * tc['R2'] * tc['C2'])
        Z_warburg = tc['Aw'] / np.sqrt(1j * omega + 1e-15)
        Z = tc['Rs'] + Z_sei + Z_ct + Z_warburg

    elif tc['model'] == 'fuel_cell':
        # Finite Warburg
        s = np.sqrt(1j * omega * tc['tau_d'] + 1e-15)
        Z_warburg = tc['R_d'] * np.tanh(s) / s
        Z_faradaic = tc['Rct'] + Z_warburg
        Y_parallel = 1j * omega * tc['Cdl'] + 1 / Z_faradaic
        Z = tc['Rs'] + 1 / Y_parallel

    return Z


def run_eis_benchmarks(verbose: bool = True) -> List[BenchmarkResult]:
    """Run electrochemical impedance benchmarks."""
    results = []
    test_cases = generate_eis_data()

    print("\n" + "="*70)
    print("CATEGORY 3: ELECTROCHEMICAL IMPEDANCE SPECTROSCOPY")
    print("="*70)

    for tc in test_cases:
        print(f"\nTest: {tc['name']}")
        print(f"  {tc['description']}")

        # Generate frequency sweep (typical EIS range)
        freqs = np.logspace(-2, 6, 100)  # 0.01 Hz to 1 MHz

        # Generate true data
        Z_true = compute_eis_impedance(tc, freqs)

        # Configure URN (emphasize diffusion)
        config = URNConfig(
            n_debye=3,
            n_cole_cole=2,
            n_cole_davidson=1,
            n_havriliak_negami=0,
            n_cpe=2,
            n_warburg=2,
            n_gerischer=1,
            n_rlc=0,
            n_skin_effect=0,
            sparsity_weight=0.01,
            n_epochs=5000,
            n_restarts=3
        )

        # Train
        start_time = time.time()
        model = train_urn(freqs, Z_true, config, verbose=False)
        train_time = time.time() - start_time

        # Evaluate
        omega = 2 * np.pi * freqs
        omega_torch = torch.tensor(omega, dtype=torch.float64)
        with torch.no_grad():
            Z_pred = model(omega_torch).numpy()

        rel_err = np.abs(Z_pred - Z_true) / np.abs(Z_true)

        # Get active components
        active = model.get_active_components()
        n_active = sum(len(v) for v in active.values())

        # Convert to serializable
        active_serializable = {}
        for key, comps in active.items():
            active_serializable[key] = []
            for c in comps:
                comp_dict = {}
                for k, v in c.items():
                    if hasattr(v, 'item'):
                        comp_dict[k] = float(v.item())
                    elif isinstance(v, (np.floating, np.integer)):
                        comp_dict[k] = float(v)
                    else:
                        comp_dict[k] = v
                active_serializable[key].append(comp_dict)

        result = BenchmarkResult(
            test_name=tc['name'],
            category='eis',
            true_model=tc['model'],
            true_params={k: v for k, v in tc.items() if k not in ['name', 'description', 'model']},
            max_rel_error=float(rel_err.max()),
            mean_rel_error=float(rel_err.mean()),
            rms_error=float(np.sqrt(np.mean(rel_err**2))),
            active_components=active_serializable,
            n_active=n_active,
            train_time_sec=train_time,
            freqs=freqs.tolist(),
            Z_true_real=Z_true.real.tolist(),
            Z_true_imag=Z_true.imag.tolist(),
            Z_pred_real=Z_pred.real.tolist(),
            Z_pred_imag=Z_pred.imag.tolist()
        )
        results.append(result)

        # Check if Warburg was identified
        warburg_detected = 'warburg' in active

        print(f"  Max error: {rel_err.max()*100:.2f}%, Mean: {rel_err.mean()*100:.2f}%")
        print(f"  Active: {n_active} ({list(active.keys())})")
        print(f"  Warburg detected: {'YES' if warburg_detected else 'NO'}")
        print(f"  Train time: {train_time:.1f} s")

    return results


# =============================================================================
# CATEGORY 4: DIELECTRIC RELAXATION
# =============================================================================

def generate_dielectric_data():
    """Generate dielectric relaxation test cases."""
    test_cases = []

    # 4.1 Polymer (Havriliak-Negami)
    test_cases.append({
        'name': 'polymer_hn',
        'eps_s': 10,
        'eps_inf': 2.5,
        'tau': 1e-6,
        'alpha': 0.8,
        'beta': 0.6,
        'model': 'havriliak_negami',
        'description': 'Polymer dielectric (Havriliak-Negami)'
    })

    # 4.2 Glass (Cole-Cole)
    test_cases.append({
        'name': 'glass_cole_cole',
        'eps_s': 8,
        'eps_inf': 2.2,
        'tau': 1e-9,
        'alpha': 0.7,
        'model': 'cole_cole',
        'description': 'Glass dielectric (Cole-Cole)'
    })

    # 4.3 Ceramic (Debye + conductivity)
    test_cases.append({
        'name': 'ceramic_debye',
        'eps_s': 1000,
        'eps_inf': 100,
        'tau': 1e-7,
        'sigma_dc': 1e-8,  # DC conductivity
        'model': 'debye_conductive',
        'description': 'Ceramic dielectric with DC conductivity'
    })

    # 4.4 Water (two Debye)
    test_cases.append({
        'name': 'water_two_debye',
        'eps_s': 80,
        'eps_1': 6,
        'eps_inf': 4,
        'tau_1': 8e-12,  # 8 ps
        'tau_2': 1e-13,  # 0.1 ps
        'model': 'two_debye',
        'description': 'Water (two Debye relaxations)'
    })

    return test_cases


def compute_dielectric_permittivity(tc: Dict, freqs: np.ndarray) -> np.ndarray:
    """Compute complex permittivity for various models."""
    omega = 2 * np.pi * freqs
    eps_0 = 8.854e-12

    if tc['model'] == 'havriliak_negami':
        eps = tc['eps_inf'] + (tc['eps_s'] - tc['eps_inf']) / (1 + (1j * omega * tc['tau'])**tc['alpha'])**tc['beta']

    elif tc['model'] == 'cole_cole':
        eps = tc['eps_inf'] + (tc['eps_s'] - tc['eps_inf']) / (1 + (1j * omega * tc['tau'])**tc['alpha'])

    elif tc['model'] == 'debye_conductive':
        eps_relax = tc['eps_inf'] + (tc['eps_s'] - tc['eps_inf']) / (1 + 1j * omega * tc['tau'])
        eps = eps_relax - 1j * tc['sigma_dc'] / (omega * eps_0 + 1e-20)

    elif tc['model'] == 'two_debye':
        eps = tc['eps_inf'] + (tc['eps_1'] - tc['eps_inf']) / (1 + 1j * omega * tc['tau_2'])
        eps += (tc['eps_s'] - tc['eps_1']) / (1 + 1j * omega * tc['tau_1'])

    return eps


def run_dielectric_benchmarks(verbose: bool = True) -> List[BenchmarkResult]:
    """Run dielectric relaxation benchmarks."""
    results = []
    test_cases = generate_dielectric_data()

    print("\n" + "="*70)
    print("CATEGORY 4: DIELECTRIC RELAXATION")
    print("="*70)

    for tc in test_cases:
        print(f"\nTest: {tc['name']}")
        print(f"  {tc['description']}")

        # Generate frequency sweep
        freqs = np.logspace(3, 12, 100)  # 1 kHz to 1 THz

        # Generate true data
        eps_true = compute_dielectric_permittivity(tc, freqs)

        # Configure URN (emphasize Debye family)
        config = URNConfig(
            n_debye=3,
            n_cole_cole=2,
            n_cole_davidson=1,
            n_havriliak_negami=2,
            n_cpe=1,
            n_warburg=0,
            n_gerischer=0,
            n_rlc=0,
            n_skin_effect=0,
            sparsity_weight=0.005,
            n_epochs=4000,
            n_restarts=3
        )

        # Train
        start_time = time.time()
        model = train_urn(freqs, eps_true, config, verbose=False)
        train_time = time.time() - start_time

        # Evaluate
        omega = 2 * np.pi * freqs
        omega_torch = torch.tensor(omega, dtype=torch.float64)
        with torch.no_grad():
            eps_pred = model(omega_torch).numpy()

        rel_err = np.abs(eps_pred - eps_true) / np.abs(eps_true)

        # Get active components
        active = model.get_active_components()
        n_active = sum(len(v) for v in active.values())

        # Convert to serializable
        active_serializable = {}
        for key, comps in active.items():
            active_serializable[key] = []
            for c in comps:
                comp_dict = {}
                for k, v in c.items():
                    if hasattr(v, 'item'):
                        comp_dict[k] = float(v.item())
                    elif isinstance(v, (np.floating, np.integer)):
                        comp_dict[k] = float(v)
                    else:
                        comp_dict[k] = v
                active_serializable[key].append(comp_dict)

        result = BenchmarkResult(
            test_name=tc['name'],
            category='dielectric',
            true_model=tc['model'],
            true_params={k: v for k, v in tc.items() if k not in ['name', 'description', 'model']},
            max_rel_error=float(rel_err.max()),
            mean_rel_error=float(rel_err.mean()),
            rms_error=float(np.sqrt(np.mean(rel_err**2))),
            active_components=active_serializable,
            n_active=n_active,
            train_time_sec=train_time,
            freqs=freqs.tolist(),
            Z_true_real=eps_true.real.tolist(),
            Z_true_imag=eps_true.imag.tolist(),
            Z_pred_real=eps_pred.real.tolist(),
            Z_pred_imag=eps_pred.imag.tolist()
        )
        results.append(result)

        print(f"  Max error: {rel_err.max()*100:.2f}%, Mean: {rel_err.mean()*100:.2f}%")
        print(f"  Active: {n_active} ({list(active.keys())})")
        print(f"  Train time: {train_time:.1f} s")

    return results


# =============================================================================
# MAIN: RUN ALL BENCHMARKS
# =============================================================================

def run_all_benchmarks():
    """Run complete benchmark suite."""
    print("="*70)
    print("KAN-inspired URN Benchmark Suite")
    print("="*70)
    print(f"Start time: {datetime.now().isoformat()}")

    all_results = []

    # Category 1: Ferrite
    ferrite_results = run_ferrite_benchmarks()
    all_results.extend(ferrite_results)

    # Category 2: Skin effect
    skin_results = run_skin_effect_benchmarks()
    all_results.extend(skin_results)

    # Category 3: EIS
    eis_results = run_eis_benchmarks()
    all_results.extend(eis_results)

    # Category 4: Dielectric
    dielectric_results = run_dielectric_benchmarks()
    all_results.extend(dielectric_results)

    # Save results
    output_dir = os.path.dirname(os.path.abspath(__file__))
    output_file = os.path.join(output_dir, 'urn_benchmark_results.json')
    save_results(all_results, output_file)

    # Print summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)

    categories = {}
    for r in all_results:
        if r.category not in categories:
            categories[r.category] = []
        categories[r.category].append(r)

    print(f"\n{'Category':<20} {'Tests':<8} {'Avg Max Err':<12} {'Avg Mean Err':<12} {'Avg Time':<10}")
    print("-"*70)

    for cat, results in categories.items():
        n_tests = len(results)
        avg_max = np.mean([r.max_rel_error for r in results]) * 100
        avg_mean = np.mean([r.mean_rel_error for r in results]) * 100
        avg_time = np.mean([r.train_time_sec for r in results])
        print(f"{cat:<20} {n_tests:<8} {avg_max:<12.2f}% {avg_mean:<12.2f}% {avg_time:<10.1f}s")

    print("-"*70)
    print(f"{'TOTAL':<20} {len(all_results):<8}")

    # Mechanism detection summary
    print("\n\nMECHANISM DETECTION:")
    print("-"*70)

    mechanism_counts = {}
    for r in all_results:
        for mech in r.active_components.keys():
            if mech not in mechanism_counts:
                mechanism_counts[mech] = 0
            mechanism_counts[mech] += 1

    for mech, count in sorted(mechanism_counts.items(), key=lambda x: -x[1]):
        print(f"  {mech:<20}: detected in {count}/{len(all_results)} tests")

    print(f"\nResults saved to: {output_file}")
    print(f"End time: {datetime.now().isoformat()}")

    return all_results


if __name__ == '__main__':
    run_all_benchmarks()
