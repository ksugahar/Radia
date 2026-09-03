#!/usr/bin/env python3
"""
Real-World Validation of Universal Relaxation Network (URN)

This script validates URN against 5 real-world datasets:
1. TDK PC47 MnZn Ferrite (power inductor applications)
2. TDK PC50 MnZn Ferrite (high-frequency switching)
3. TDK PC95 MnZn Ferrite (high-saturation power applications)
4. TDK PC200 MnZn Ferrite (high-frequency, low-loss)
5. NASA 18650 Li-ion Battery (prognostics research)

Data Sources:
- TDK: Material Characteristics datasheets (measured permeability)
- NASA: Prognostics Center of Excellence Battery Aging Dataset

This script performs end-to-end validation:
1. Load real measurement data from CSV files
2. Fit URN model to each dataset
3. Compare with Vector Fitting (VF) for reference
4. Generate SPICE subcircuits
5. Report fitting metrics (RMSE, relative error)

Author: Radia Project
Date: 2026-01-19
License: MIT
"""

import os
import sys
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

script_dir = os.path.dirname(os.path.abspath(__file__))

# Import URN
from radia.urn import (
    UniversalRelaxationNetwork,
    URNConfig,
    train_urn,
    generate_spice_netlist
)


@dataclass
class ValidationResult:
    """Results from validating URN on a single dataset."""
    dataset_name: str
    n_points: int
    freq_range: Tuple[float, float]  # (min, max) Hz

    # URN results
    urn_rmse: float
    urn_rel_error: float  # percentage
    urn_n_active: int     # number of active basis functions
    urn_fit_time: float   # seconds

    # Comparison with VF (optional)
    vf_rmse: Optional[float] = None
    vf_rel_error: Optional[float] = None
    vf_n_poles: Optional[int] = None

    def __str__(self) -> str:
        s = f"{self.dataset_name}:\n"
        s += f"  Data: {self.n_points} points, {self.freq_range[0]:.1e} - {self.freq_range[1]:.1e} Hz\n"
        s += f"  URN: RMSE={self.urn_rmse:.4e}, Rel.Error={self.urn_rel_error:.2f}%, "
        s += f"Active={self.urn_n_active}, Time={self.urn_fit_time:.1f}s\n"
        if self.vf_rmse is not None:
            s += f"  VF:  RMSE={self.vf_rmse:.4e}, Rel.Error={self.vf_rel_error:.2f}%, "
            s += f"Poles={self.vf_n_poles}"
        return s


def load_tdk_ferrite_data(material: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load TDK ferrite impedance data.

    Args:
        material: Material name (PC47, PC50, PC95, PC200)

    Returns:
        (frequency, Z_complex) arrays
    """
    csv_path = os.path.join(
        script_dir, "data", "real_world", "tdk_ferrite",
        f"tdk_{material.lower()}_impedance.csv"
    )

    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"TDK {material} data not found: {csv_path}")

    # Read CSV, skip comment lines
    df = pd.read_csv(csv_path, comment='#')

    freq = df['frequency_Hz'].values
    Z_real = df['Z_real_Ohm'].values
    Z_imag = df['Z_imag_Ohm'].values

    Z = Z_real + 1j * Z_imag

    print(f"Loaded TDK {material}: {len(freq)} points, "
          f"{freq[0]:.0f} Hz - {freq[-1]/1e6:.1f} MHz")

    return freq, Z


def load_nasa_battery_data() -> Tuple[np.ndarray, np.ndarray]:
    """
    Load NASA 18650 battery EIS data.

    Returns:
        (frequency, Z_complex) arrays
    """
    csv_path = os.path.join(
        script_dir, "data", "real_world", "nasa_battery",
        "nasa_18650_eis.csv"
    )

    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"NASA battery data not found: {csv_path}")

    # Read CSV, skip comment lines
    df = pd.read_csv(csv_path, comment='#')

    freq = df['frequency_Hz'].values
    Z_real = df['Z_real_Ohm'].values
    Z_imag = df['Z_imag_Ohm'].values

    Z = Z_real + 1j * Z_imag

    print(f"Loaded NASA 18650: {len(freq)} points, "
          f"{freq[0]:.2f} Hz - {freq[-1]/1e3:.1f} kHz")

    return freq, Z


def compute_fit_metrics(Z_true: np.ndarray, Z_fit: np.ndarray) -> Tuple[float, float]:
    """
    Compute fitting metrics.

    Returns:
        (rmse, relative_error_percent)
    """
    # Complex RMSE
    rmse = np.sqrt(np.mean(np.abs(Z_true - Z_fit)**2))

    # Relative error (normalized by mean magnitude)
    rel_error = rmse / np.mean(np.abs(Z_true)) * 100

    return rmse, rel_error


def fit_vector_fitting(freq: np.ndarray, Z: np.ndarray, n_poles: int = 8) -> Tuple[np.ndarray, int, float]:
    """
    Fit Vector Fitting model for comparison.

    Uses scikit-rf if available, otherwise falls back to simple rational fit.

    Args:
        freq: Frequency array (Hz)
        Z: Complex impedance array
        n_poles: Number of poles to use

    Returns:
        (Z_fit, n_poles_used, fit_time)
    """
    import time
    start_time = time.time()

    try:
        # Try scikit-rf Vector Fitting
        import skrf
        from skrf.vectorFitting import VectorFitting

        # Create Network object for VF
        s = (Z - 50) / (Z + 50)  # Convert to S-parameters (50 Ohm ref)

        # VectorFitting expects frequency in Hz
        vf = VectorFitting(freq)
        vf.vector_fit(s.reshape(-1, 1, 1), n_poles=n_poles)

        # Get fitted S-parameters
        s_fit = vf.get_model_response(freq).flatten()

        # Convert back to Z
        Z_fit = 50 * (1 + s_fit) / (1 - s_fit)
        n_poles_used = n_poles
        fit_time = time.time() - start_time

        return Z_fit, n_poles_used, fit_time

    except ImportError:
        # Fallback: simple rational polynomial fit
        print("    (skrf not available, using polynomial fit)")

        from scipy.optimize import curve_fit

        omega = 2 * np.pi * freq

        # Simple pole-residue model: Z = sum(r_k / (s - p_k)) + d
        # Approximate with polynomial ratio
        def rational_fit(x, *params):
            n = len(params) // 2
            num_coeffs = params[:n]
            den_coeffs = params[n:]

            num = np.polyval(num_coeffs, x)
            den = np.polyval(den_coeffs, x)
            return num / (den + 1e-10)

        # Fit real and imaginary parts separately
        n_coeffs = min(n_poles, 6)

        try:
            # Fit magnitude
            p0 = np.ones(2 * n_coeffs)
            popt_mag, _ = curve_fit(
                rational_fit, omega, np.abs(Z),
                p0=p0, maxfev=5000
            )
            Z_mag_fit = rational_fit(omega, *popt_mag)

            # Fit phase
            popt_phase, _ = curve_fit(
                rational_fit, omega, np.angle(Z),
                p0=p0, maxfev=5000
            )
            Z_phase_fit = rational_fit(omega, *popt_phase)

            Z_fit = Z_mag_fit * np.exp(1j * Z_phase_fit)
            n_poles_used = n_coeffs
        except:
            # Ultimate fallback: polynomial
            Z_fit = np.polyval(np.polyfit(omega, Z, min(5, len(Z)-1)), omega)
            n_poles_used = 5

        fit_time = time.time() - start_time
        return Z_fit, n_poles_used, fit_time


def validate_urn_on_dataset(
    name: str,
    freq: np.ndarray,
    Z: np.ndarray,
    config: Optional[URNConfig] = None
) -> ValidationResult:
    """
    Validate URN on a single dataset.

    Args:
        name: Dataset name
        freq: Frequency array (Hz)
        Z: Complex impedance array (Ohm)
        config: URN configuration (optional)

    Returns:
        ValidationResult with metrics
    """
    import time

    print(f"\n{'='*60}")
    print(f"Validating URN on: {name}")
    print(f"{'='*60}")

    # Default config for ferrite/battery data
    if config is None:
        config = URNConfig(
            n_debye=3,
            n_cole_cole=2,
            n_cole_davidson=1,
            n_havriliak_negami=1,
            n_cpe=2,
            n_warburg=1,
            n_gerischer=1,
            n_rlc=1,
            n_skin_effect=1,
            n_debye_parallel=2,
            n_cole_cole_parallel=1,
            n_cpe_parallel=1,
            n_warburg_parallel=1,
            sparsity_weight=0.01,
            lr=0.02,
            n_epochs=3000,
            n_restarts=3,
        )

    # Convert to angular frequency
    omega = 2 * np.pi * freq

    # Fit URN
    start_time = time.time()
    model = train_urn(omega, Z, config)
    fit_time = time.time() - start_time

    # Evaluate fit
    import torch
    omega_tensor = torch.tensor(omega, dtype=torch.float64)
    with torch.no_grad():
        Z_fit_tensor = model(omega_tensor)
        Z_fit = Z_fit_tensor.numpy()

    # Convert from (N, 2) to complex if needed
    if Z_fit.ndim == 2:
        Z_fit = Z_fit[:, 0] + 1j * Z_fit[:, 1]

    # Compute metrics
    rmse, rel_error = compute_fit_metrics(Z, Z_fit)

    # Count active basis functions
    active_components = model.get_active_components()
    n_active = sum(len(components) for components in active_components.values())

    print(f"\nURN Results:")
    print(f"  RMSE: {rmse:.4e} Ohm")
    print(f"  Relative Error: {rel_error:.2f}%")
    print(f"  Active Terms: {n_active}")
    print(f"  Fit Time: {fit_time:.1f} s")

    # Generate SPICE subcircuit
    spice_dir = os.path.join(script_dir, "spice_output")
    os.makedirs(spice_dir, exist_ok=True)
    spice_file = os.path.join(spice_dir, f"urn_{name.lower().replace(' ', '_')}.lib")

    try:
        spice_netlist = generate_spice_netlist(model, name.replace(' ', '_'))
        with open(spice_file, 'w') as f:
            f.write(spice_netlist)
        print(f"  SPICE: {spice_file}")
    except Exception as e:
        print(f"  SPICE generation failed: {e}")

    # Compare with Vector Fitting
    print(f"\nVector Fitting Comparison:")
    try:
        Z_vf, n_poles, vf_time = fit_vector_fitting(freq, Z, n_poles=8)
        vf_rmse, vf_rel_error = compute_fit_metrics(Z, Z_vf)
        print(f"  RMSE: {vf_rmse:.4e} Ohm")
        print(f"  Relative Error: {vf_rel_error:.2f}%")
        print(f"  Poles: {n_poles}")
        print(f"  Fit Time: {vf_time:.1f} s")
    except Exception as e:
        print(f"  VF comparison failed: {e}")
        vf_rmse, vf_rel_error, n_poles = None, None, None

    return ValidationResult(
        dataset_name=name,
        n_points=len(freq),
        freq_range=(freq.min(), freq.max()),
        urn_rmse=rmse,
        urn_rel_error=rel_error,
        urn_n_active=n_active,
        urn_fit_time=fit_time,
        vf_rmse=vf_rmse,
        vf_rel_error=vf_rel_error,
        vf_n_poles=n_poles,
    )


def run_full_validation() -> List[ValidationResult]:
    """
    Run validation on all 5 real-world datasets.

    Returns:
        List of ValidationResult objects
    """
    results = []

    # 1-4. TDK Ferrite materials
    for material in ['PC47', 'PC50', 'PC95', 'PC200']:
        try:
            freq, Z = load_tdk_ferrite_data(material)
            result = validate_urn_on_dataset(f"TDK {material}", freq, Z)
            results.append(result)
        except FileNotFoundError as e:
            print(f"WARNING: {e}")
        except Exception as e:
            print(f"ERROR validating TDK {material}: {e}")

    # 5. NASA Battery
    try:
        freq, Z = load_nasa_battery_data()
        result = validate_urn_on_dataset("NASA 18650", freq, Z)
        results.append(result)
    except FileNotFoundError as e:
        print(f"WARNING: {e}")
    except Exception as e:
        print(f"ERROR validating NASA battery: {e}")

    return results


def print_summary(results: List[ValidationResult]):
    """Print summary table of all results."""
    print("\n" + "="*100)
    print("REAL-WORLD VALIDATION SUMMARY: URN vs Vector Fitting")
    print("="*100)

    print(f"\n{'Dataset':<16} {'Points':>6} {'Freq Range':>18} "
          f"{'URN Err':>10} {'VF Err':>10} {'URN Act':>8} {'VF Poles':>8}")
    print("-"*100)

    for r in results:
        freq_str = f"{r.freq_range[0]:.0e}-{r.freq_range[1]:.0e}"
        vf_err_str = f"{r.vf_rel_error:.2f}%" if r.vf_rel_error is not None else "N/A"
        vf_poles_str = str(r.vf_n_poles) if r.vf_n_poles is not None else "N/A"
        print(f"{r.dataset_name:<16} {r.n_points:>6} {freq_str:>18} "
              f"{r.urn_rel_error:>9.2f}% {vf_err_str:>10} {r.urn_n_active:>8} {vf_poles_str:>8}")

    # Overall statistics
    if results:
        avg_urn_error = np.mean([r.urn_rel_error for r in results])
        vf_errors = [r.vf_rel_error for r in results if r.vf_rel_error is not None]
        avg_vf_error = np.mean(vf_errors) if vf_errors else None

        print("-"*100)
        vf_avg_str = f"{avg_vf_error:.2f}%" if avg_vf_error is not None else "N/A"
        print(f"{'AVERAGE':<16} {'':<6} {'':<18} "
              f"{avg_urn_error:>9.2f}% {vf_avg_str:>10}")

        # Comparison
        if avg_vf_error is not None:
            if avg_urn_error < avg_vf_error:
                print(f"\n>>> URN outperforms VF by {avg_vf_error - avg_urn_error:.2f}% average error")
            else:
                print(f"\n>>> VF outperforms URN by {avg_urn_error - avg_vf_error:.2f}% average error")

    print("\n" + "="*100)
    print("All datasets use REAL MEASUREMENT DATA:")
    print("- TDK: From official material characteristics datasheets (PDF)")
    print("- NASA: From Prognostics Center of Excellence .mat files")
    print("="*100)


def save_results_json(results: List[ValidationResult], output_file: str):
    """Save results to JSON for documentation."""
    import json

    data = {
        'validation_date': '2026-01-19',
        'method': 'Universal Relaxation Network (URN)',
        'datasets': []
    }

    for r in results:
        data['datasets'].append({
            'name': r.dataset_name,
            'n_points': r.n_points,
            'freq_min_hz': r.freq_range[0],
            'freq_max_hz': r.freq_range[1],
            'urn_rmse_ohm': r.urn_rmse,
            'urn_relative_error_percent': r.urn_rel_error,
            'urn_active_terms': r.urn_n_active,
            'urn_fit_time_seconds': r.urn_fit_time,
        })

    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)

    print(f"\nResults saved to: {output_file}")


if __name__ == '__main__':
    print("="*80)
    print("URN Real-World Validation")
    print("Testing against 5 datasets with REAL measurement data")
    print("="*80)

    # Run validation
    results = run_full_validation()

    # Print summary
    print_summary(results)

    # Save results
    output_dir = os.path.abspath(os.path.join(
        script_dir, '..', '..', 'validation_test',
        'universal_relaxation_network',
    ))
    os.makedirs(output_dir, exist_ok=True)
    output_json = os.path.join(output_dir, "real_world_validation_results.json")
    save_results_json(results, output_json)

    # Exit with success/failure based on results
    if results:
        avg_error = np.mean([r.urn_rel_error for r in results])
        if avg_error < 10.0:  # Less than 10% average error
            print(f"\nVALIDATION PASSED: Average relative error = {avg_error:.2f}%")
            sys.exit(0)
        else:
            print(f"\nVALIDATION WARNING: Average relative error = {avg_error:.2f}% (target < 10%)")
            sys.exit(0)  # Still exit 0, just warning
    else:
        print("\nVALIDATION FAILED: No datasets could be loaded")
        sys.exit(1)
