#!/usr/bin/env python
"""
Comprehensive URN vs Vector Fitting Validation on All Real-World Datasets

This script validates URN on all available real-world datasets:
1. NASA 18650 Battery EIS
2. TDK PC47 Ferrite
3. TDK PC50 Ferrite
4. TDK PC95 Ferrite
5. TDK PC200 Ferrite

Output:
- Comprehensive comparison table
- Individual plots for each dataset
- JSON results file

Usage:
    python validation_test/universal_relaxation_network/validate_all_datasets.py
"""

import os
import json
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

from radia.urn import (
    UniversalRelaxationNetwork, URNConfig, train_urn, generate_spice_netlist
)
import torch

DOCS_URN_DIR = (
    Path(__file__).resolve().parents[2] / 'docs' / 'universal_relaxation_network'
)


def load_data(data_path):
    """Load impedance data from CSV file."""
    data = pd.read_csv(data_path, comment='#')
    freq = data['frequency_Hz'].to_numpy()
    Z_real = data['Z_real_Ohm'].to_numpy()
    Z_imag = data['Z_imag_Ohm'].to_numpy()
    return freq, Z_real + 1j * Z_imag


def compute_nrmse(Z_true, Z_pred):
    """Compute Normalized Root Mean Square Error."""
    error = np.abs(Z_true - Z_pred)
    nrmse = np.sqrt(np.mean(error**2)) / np.mean(np.abs(Z_true))
    return nrmse


def vector_fitting(freq, Z_data, n_poles=10, n_iter=5):
    """Simple Vector Fitting implementation."""
    omega = 2 * np.pi * freq
    s = 1j * omega

    # Initialize poles (logarithmically spaced)
    f_min, f_max = freq.min(), freq.max()
    pole_freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_poles // 2)

    # Real poles
    poles_real = -2 * np.pi * pole_freqs[:n_poles//4]

    # Complex conjugate poles
    poles_complex = []
    for f in pole_freqs[n_poles//4:]:
        Q = 2.0
        omega_p = 2 * np.pi * f
        sigma_p = -omega_p / (2 * Q)
        omega_d = omega_p * np.sqrt(1 - 1/(4*Q**2))
        poles_complex.append(sigma_p + 1j * omega_d)
        poles_complex.append(sigma_p - 1j * omega_d)

    poles = np.concatenate([poles_real, np.array(poles_complex)])

    # Iterative pole relocation
    for iteration in range(n_iter):
        # Build basis matrix
        A = np.zeros((len(s), len(poles) + 1), dtype=complex)
        A[:, 0] = 1  # d term
        for i, p in enumerate(poles):
            A[:, i + 1] = 1 / (s - p)

        # Solve least squares
        coeffs, residuals, rank, sv = np.linalg.lstsq(A, Z_data, rcond=None)

        # Compute fitted response
        Z_fit = A @ coeffs

    return Z_fit


def run_urn_fitting(freq, Z_data, n_restarts=3):
    """Run URN fitting on impedance data."""
    config = URNConfig(
        n_debye=3,
        n_cole_cole=2,
        n_warburg=2,
        sparsity_weight=0.01,
        lr=0.02,
        n_epochs=4000,
        n_restarts=n_restarts
    )

    start_time = time.time()
    model = train_urn(freq, Z_data, config)
    train_time = time.time() - start_time

    # Get fitted response
    omega = 2 * np.pi * freq
    omega_tensor = torch.tensor(omega, dtype=torch.float64)
    Z_fit = model(omega_tensor).detach().numpy()

    # Count active parameters
    active_components = model.get_active_components()
    active_params = sum(len(v) for v in active_components.values())

    return Z_fit, model, active_params, train_time


def main():
    raise RuntimeError('Legacy bundled NASA mode is retired: this all-dataset aggregate depended on private measurements with an unverified frequency axis. Use the documented private-input consumers in docs/universal_relaxation_network.')


if __name__ == '__main__':
    main()
