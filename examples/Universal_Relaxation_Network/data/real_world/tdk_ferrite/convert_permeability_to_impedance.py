#!/usr/bin/env python3
"""
Convert TDK ferrite complex permeability data to impedance Z(f).

This script converts the complex permeability (mu' - j*mu") from TDK datasheets
to complex impedance Z(f) for a toroidal inductor geometry.

Source: TDK Corporation datasheet "Mn-Zn Ferrite Material characteristics" (May 2022)
Document: ferrite_mn-zn_material_characteristics_en.pdf

Physics:
    Complex permeability: mu = mu' - j*mu"
    Complex inductance: L = L0 * (mu' - j*mu") where L0 = mu0 * N^2 * Ae / le
    Complex impedance: Z = j*omega*L = j*omega*L0*(mu' - j*mu")
                        = omega*L0*mu" + j*omega*L0*mu'
                        = R(f) + j*X(f)

    Where:
        R(f) = omega * L0 * mu"  (resistance from magnetic losses)
        X(f) = omega * L0 * mu'  (reactance from energy storage)

Reference toroidal core geometry (from TDK datasheet test conditions):
    OD = 31 mm, ID = 19 mm, TH = 8 mm
    Ae = TH * (OD - ID) / 2 = 8 * 6 = 48 mm^2 = 48e-6 m^2
    le = pi * (OD + ID) / 2 = pi * 25 = 78.54 mm = 0.07854 m
    N = 10 turns (assumed for demonstration)
    L0 = mu0 * N^2 * Ae / le = 4*pi*1e-7 * 100 * 48e-6 / 0.07854 = 7.68e-8 H

Author: Radia Project
Date: 2026-01-19
"""

import numpy as np
import pandas as pd
import os

# Physical constants
MU0 = 4 * np.pi * 1e-7  # H/m

# Toroidal core geometry (TDK test core: OD=31mm, ID=19mm, TH=8mm)
OD = 31e-3  # m
ID = 19e-3  # m
TH = 8e-3   # m
N_TURNS = 10  # Number of turns (assumed)

# Calculate core parameters
Ae = TH * (OD - ID) / 2  # Effective cross-sectional area [m^2]
le = np.pi * (OD + ID) / 2  # Effective magnetic path length [m]
L0 = MU0 * N_TURNS**2 * Ae / le  # Air-core inductance [H]

print(f"Toroidal core geometry:")
print(f"  OD = {OD*1000:.1f} mm, ID = {ID*1000:.1f} mm, TH = {TH*1000:.1f} mm")
print(f"  Ae = {Ae*1e6:.2f} mm^2")
print(f"  le = {le*1000:.2f} mm")
print(f"  N = {N_TURNS} turns")
print(f"  L0 = {L0*1e9:.3f} nH")


def convert_permeability_to_impedance(input_file: str, output_file: str, material_name: str):
    """
    Convert permeability CSV to impedance CSV.

    Args:
        input_file: Path to permeability CSV (frequency_Hz, mu_real, mu_imag)
        output_file: Path to output impedance CSV
        material_name: Material name for header
    """
    # Read permeability data
    df = pd.read_csv(input_file, comment='#')

    freq = df['frequency_Hz'].values
    mu_real = df['mu_real'].values
    mu_imag = df['mu_imag'].values

    # Calculate angular frequency
    omega = 2 * np.pi * freq

    # Calculate impedance components
    # Z = j*omega*L = j*omega*L0*(mu' - j*mu")
    #   = omega*L0*mu" + j*omega*L0*mu'
    R = omega * L0 * mu_imag  # Resistance (from losses)
    X = omega * L0 * mu_real  # Reactance (from inductance)

    # Complex impedance
    Z_real = R
    Z_imag = X
    Z_mag = np.sqrt(Z_real**2 + Z_imag**2)
    Z_phase_deg = np.degrees(np.arctan2(Z_imag, Z_real))

    # Create output dataframe
    output_df = pd.DataFrame({
        'frequency_Hz': freq,
        'Z_real_Ohm': Z_real,
        'Z_imag_Ohm': Z_imag,
        'Z_mag_Ohm': Z_mag,
        'Z_phase_deg': Z_phase_deg
    })

    # Write with header
    header = f"""# TDK MnZn Ferrite {material_name} Impedance Data
# Source: TDK Corporation datasheet "Mn-Zn Ferrite Material characteristics" (May 2022)
# Document: ferrite_mn-zn_material_characteristics_en.pdf
# REAL MANUFACTURER DATA - Converted from complex permeability
#
# Conversion parameters:
#   Toroidal core: OD={OD*1000:.0f}mm, ID={ID*1000:.0f}mm, TH={TH*1000:.0f}mm
#   Number of turns: N={N_TURNS}
#   Air-core inductance: L0={L0*1e9:.3f} nH
#
# Z = j*omega*L0*(mu' - j*mu") = omega*L0*mu" + j*omega*L0*mu'
#
# frequency_Hz: Frequency in Hz
# Z_real_Ohm: Real part of impedance (resistance from magnetic losses)
# Z_imag_Ohm: Imaginary part of impedance (reactance from inductance)
# Z_mag_Ohm: Magnitude of impedance
# Z_phase_deg: Phase angle in degrees
#
"""

    with open(output_file, 'w') as f:
        f.write(header)
        output_df.to_csv(f, index=False)

    print(f"\nConverted {material_name}:")
    print(f"  Input: {input_file}")
    print(f"  Output: {output_file}")
    print(f"  Frequency range: {freq[0]/1e3:.1f} kHz - {freq[-1]/1e6:.1f} MHz")
    print(f"  |Z| range: {Z_mag[0]:.3f} - {Z_mag.max():.1f} Ohm")

    return output_df


if __name__ == '__main__':
    # Get script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # List of materials to convert
    materials = ['PC47', 'PC50', 'PC95', 'PC200']

    for material in materials:
        input_file = os.path.join(script_dir, f'tdk_{material.lower()}_permeability.csv')
        output_file = os.path.join(script_dir, f'tdk_{material.lower()}_impedance.csv')
        if os.path.exists(input_file):
            convert_permeability_to_impedance(input_file, output_file, material)

    print("\n" + "="*60)
    print("Conversion complete!")
    print("These impedance files can be used with URN for validation.")
    print("="*60)
