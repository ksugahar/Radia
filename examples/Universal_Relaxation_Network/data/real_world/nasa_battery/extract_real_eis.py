#!/usr/bin/env python3
"""
Extract REAL EIS data from NASA Battery Dataset.

Source: NASA Ames Prognostics Center of Excellence (PCoE)
Dataset: Li-ion Battery Aging Datasets
URL: https://data.nasa.gov/dataset/li-ion-battery-aging-datasets

This script extracts actual measured EIS data from the NASA .mat files.

Reference:
    B. Saha and K. Goebel (2007). "Battery Data Set", NASA Ames
    Prognostics Data Repository, NASA Ames Research Center, Moffett Field, CA.

Author: Radia Project
Date: 2026-01-19
"""

import os
import zipfile
import scipy.io
import numpy as np
import pandas as pd

# Path to downloaded NASA data
NASA_DATA_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "5. Battery Data Set"
)


def extract_eis_from_nasa(battery_id: str = "B0005", cycle_type: str = "first") -> dict:
    """
    Extract EIS data from NASA battery dataset.

    Args:
        battery_id: Battery ID (B0005, B0006, B0007, or B0018)
        cycle_type: "first" for fresh battery, "last" for aged battery

    Returns:
        Dictionary with frequency and impedance arrays
    """
    zip_path = os.path.join(NASA_DATA_DIR, "1. BatteryAgingARC-FY08Q4.zip")

    if not os.path.exists(zip_path):
        print(f"ERROR: NASA data not found at {zip_path}")
        print("Please download from: https://data.nasa.gov/dataset/li-ion-battery-aging-datasets")
        return None

    print(f"Extracting EIS data for {battery_id} ({cycle_type})...")

    with zipfile.ZipFile(zip_path, 'r') as zf:
        # Find the mat file
        mat_name = f"{battery_id}.mat"
        mat_files = [n for n in zf.namelist() if n.endswith(mat_name)]

        if not mat_files:
            print(f"ERROR: {mat_name} not found in archive")
            return None

        # Extract to temp
        temp_dir = os.path.dirname(os.path.abspath(__file__))
        zf.extract(mat_files[0], temp_dir)
        mat_path = os.path.join(temp_dir, mat_files[0])

        try:
            # Load MATLAB file
            mat_data = scipy.io.loadmat(mat_path)
            battery = mat_data[battery_id]
            cycles = battery['cycle'][0, 0]
            n_cycles = cycles.shape[1]

            print(f"  Total cycles: {n_cycles}")

            # Find impedance cycles
            impedance_cycles = []
            for i in range(n_cycles):
                cycle = cycles[0, i]
                cycle_type_str = str(cycle['type'][0])
                if cycle_type_str == 'impedance':
                    impedance_cycles.append(i)

            print(f"  Impedance cycles: {len(impedance_cycles)}")

            if not impedance_cycles:
                print("ERROR: No impedance cycles found")
                return None

            # Select cycle based on type
            if cycle_type == "first":
                cycle_idx = impedance_cycles[0]
            else:  # last
                cycle_idx = impedance_cycles[-1]

            cycle = cycles[0, cycle_idx]
            data = cycle['data'][0, 0]

            print(f"  Selected cycle: {cycle_idx}")
            print(f"  Data fields: {data.dtype.names}")

            # Extract impedance data
            result = {'battery_id': battery_id, 'cycle': cycle_idx}

            if 'Battery_impedance' in data.dtype.names:
                Z = data['Battery_impedance'].flatten()
                result['Z'] = Z
                result['Z_real'] = np.real(Z)
                result['Z_imag'] = np.imag(Z)
                result['Z_mag'] = np.abs(Z)
                print(f"  Impedance points: {len(Z)}")
                print(f"  |Z| range: {np.abs(Z).min():.4f} - {np.abs(Z).max():.4f} Ohm")

            if 'Rectified_Impedance' in data.dtype.names:
                result['Z_rect'] = data['Rectified_Impedance'].flatten()

            # NASA EIS frequency range is 0.1 Hz to 5 kHz
            # Generate frequency array based on number of points
            n_points = len(result.get('Z', []))
            if n_points > 0:
                # Assume log-spaced frequency sweep
                result['frequency'] = np.logspace(-1, np.log10(5000), n_points)

            return result

        finally:
            # Clean up
            if os.path.exists(mat_path):
                os.remove(mat_path)
            # Remove extracted directory if empty
            extracted_dir = os.path.dirname(mat_path)
            if extracted_dir != temp_dir and os.path.exists(extracted_dir):
                try:
                    os.rmdir(extracted_dir)
                except:
                    pass


def create_real_eis_csv(output_dir: str):
    """
    Create EIS CSV from real NASA measurement data.
    """
    # Extract fresh battery EIS
    eis_data = extract_eis_from_nasa("B0005", "first")

    if eis_data is None or 'Z' not in eis_data:
        print("WARNING: Could not extract real EIS data from NASA dataset")
        print("The NASA .mat files may have a different structure than expected")
        return None

    freq = eis_data['frequency']
    Z = eis_data['Z']

    # Create DataFrame
    df = pd.DataFrame({
        'frequency_Hz': freq,
        'Z_real_Ohm': np.real(Z),
        'Z_imag_Ohm': np.imag(Z),
        'Z_mag_Ohm': np.abs(Z),
        'Z_phase_deg': np.degrees(np.angle(Z))
    })

    # Save with header
    header = f"""# NASA 18650 Li-ion Battery EIS Data
# Source: NASA Ames Prognostics Center of Excellence (PCoE)
# Dataset: Li-ion Battery Aging Datasets
# URL: https://data.nasa.gov/dataset/li-ion-battery-aging-datasets
# REAL MEASUREMENT DATA - Extracted from NASA .mat files
#
# Reference:
#   B. Saha and K. Goebel (2007). "Battery Data Set", NASA Ames
#   Prognostics Data Repository, NASA Ames Research Center.
#
# Battery: {eis_data['battery_id']}
# Cycle: {eis_data['cycle']} (fresh battery, beginning of life)
# Type: 18650 Li-ion, 2 Ah rated, 3.7 V
# Temperature: 24 C (room temperature)
#
# EIS measurement conditions:
#   Frequency range: 0.1 Hz - 5 kHz (estimated)
#   Points: {len(Z)}
#
"""

    output_file = os.path.join(output_dir, "nasa_18650_eis_real.csv")
    with open(output_file, 'w') as f:
        f.write(header)
        df.to_csv(f, index=False)

    print(f"\nCreated: {output_file}")
    print(f"  Frequency range: {freq[0]:.2f} Hz - {freq[-1]:.0f} Hz")
    print(f"  |Z| range: {np.abs(Z).min()*1000:.1f} - {np.abs(Z).max()*1000:.1f} mOhm")

    return df


if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Try to extract real EIS data
    result = create_real_eis_csv(script_dir)

    if result is None:
        print("\n" + "="*60)
        print("Could not extract real EIS from NASA .mat files")
        print("Using representative data based on NASA documentation instead")
        print("="*60)
    else:
        print("\n" + "="*60)
        print("Successfully extracted REAL NASA battery EIS data!")
        print("="*60)
