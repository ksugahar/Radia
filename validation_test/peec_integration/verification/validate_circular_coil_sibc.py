"""
PEEC SIBC Validation: Circular Coil vs Analytical Solution

Compares Radia PEEC SIBC results with analytical formulas for circular coils.

References:
  - F.W. Grover, "Inductance Calculations", Dover, 1946
  - E.B. Rosa, "The Self-Inductance of Circular Coils of Rectangular Cross-Section", 1906

Author: Radia Development Team
Date: 2026-02-13
"""

import sys
import os
HERE = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(HERE, '../../..'))
sys.path.insert(0, os.path.join(REPO_ROOT, 'src'))
sys.path.insert(0, os.path.join(REPO_ROOT, 'src/radia'))

import numpy as np
from scipy.special import iv, ellipk, ellipe  # Modified Bessel and elliptic integrals
import matplotlib.pyplot as plt

from radia.peec_mesh_import import GMSHCenterlineReader

print("=" * 70)
print("PEEC SIBC Validation: Circular Coil vs Analytical")
print("=" * 70)

# Import PEEC module
try:
    from peec_matrices import PEECBuilder
    print("\nUsing C++ PEEC implementation with SIBC")
except ImportError as e:
    print(f"\nERROR: PEEC module not available: {e}")
    print("Build with: pwsh -NoProfile -ExecutionPolicy Bypass -File .\\Build.ps1")
    sys.exit(1)

# Physical constants
SIGMA_COPPER = 5.8e7  # S/m
MU_0 = 4 * np.pi * 1e-7  # H/m

# ============================================================================
# Read Coil Parameters and Mesh from Cubit
# ============================================================================

print("\n[Step 1] Reading Cubit-generated mesh...")

# Read coil parameters
param_file = os.path.join(os.path.dirname(__file__),
                          "cubit_mesh_generation/circular_coil_params.txt")
mesh_file = os.path.join(os.path.dirname(__file__),
                         "cubit_mesh_generation/circular_coil_centerline.msh")

if not os.path.exists(param_file):
    print(f"\nERROR: Parameter file not found: {param_file}")
    print("Run generate_circular_coil.py in Cubit first!")
    sys.exit(1)

if not os.path.exists(mesh_file):
    print(f"\nERROR: Mesh file not found: {mesh_file}")
    print("Run generate_circular_coil.py in Cubit first!")
    sys.exit(1)

# Parse parameters
coil_params = {}
with open(param_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('#') or not line:
            continue
        parts = line.split()
        if len(parts) >= 2:
            key = parts[0]
            value = float(parts[1])
            coil_params[key] = value

coil_radius = coil_params['coil_radius']  # Mean radius
wire_radius = coil_params['wire_radius']
sigma = coil_params.get('sigma', SIGMA_COPPER)
n_segments = int(coil_params.get('n_segments', 36))

# Read GMSH centerline mesh
reader = GMSHCenterlineReader(mesh_file)
nodes, edges = reader.read()
segments = reader.get_segments()

print(f"\nCoil parameters (from Cubit):")
print(f"  Coil radius: {coil_radius * 1e3:.1f} mm")
print(f"  Wire radius: {wire_radius * 1e3:.2f} mm")
print(f"  Conductivity: {sigma:.2e} S/m (copper)")
print(f"  Segments: {len(segments)}")

# Verify circumference
circumference_expected = 2.0 * np.pi * coil_radius
circumference_actual = reader.get_total_length()
circ_error = abs(circumference_actual - circumference_expected) / circumference_expected * 100

print(f"\nMesh validation:")
print(f"  Expected circumference: {circumference_expected * 1e3:.2f} mm")
print(f"  Actual circumference: {circumference_actual * 1e3:.2f} mm")
print(f"  Error: {circ_error:.4f}%")

if circ_error > 0.5:
    print(f"\nWARNING: Mesh circumference error exceeds 0.5%!")
    sys.exit(1)

# ============================================================================
# Analytical Solution: Grover's Formula for Circular Coil Inductance
# ============================================================================

def grover_inductance(radius, wire_radius):
    """
    Grover's formula for self-inductance of circular coil (single turn).

    L = mu_0 * radius * [log(8*radius/wire_radius) - 2] (approximate)

    More accurate formula uses elliptic integrals (Neumann's formula).

    Args:
        radius: Coil mean radius [m]
        wire_radius: Wire radius [m]

    Returns:
        Inductance [H]
    """
    # Approximate formula (valid for radius >> wire_radius)
    L = MU_0 * radius * (np.log(8.0 * radius / wire_radius) - 2.0)
    return L

def grover_dc_resistance(radius, wire_radius, sigma):
    """
    DC resistance of circular coil.

    R_dc = circumference / (sigma * cross_section_area)

    Args:
        radius: Coil mean radius [m]
        wire_radius: Wire radius [m]
        sigma: Conductivity [S/m]

    Returns:
        DC resistance [Ohm]
    """
    circumference = 2.0 * np.pi * radius
    cross_section_area = np.pi * wire_radius**2
    R_dc = circumference / (sigma * cross_section_area)
    return R_dc

# Compute analytical values
L_analytical = grover_inductance(coil_radius, wire_radius)
R_dc_analytical = grover_dc_resistance(coil_radius, wire_radius, sigma)

print(f"\nAnalytical solution (Grover):")
print(f"  DC inductance: {L_analytical * 1e6:.4f} uH")
print(f"  DC resistance: {R_dc_analytical * 1e3:.4f} mOhm")

# ============================================================================
# PEEC Calculation (DC and AC)
# ============================================================================

print(f"\nEquivalent square cross-section:")
equiv_side = np.sqrt(np.pi) * wire_radius
print(f"  Side length: {equiv_side * 1e3:.3f} mm x {equiv_side * 1e3:.3f} mm")
print(f"  (Area matches circular: {np.pi * wire_radius**2 * 1e6:.4f} mm^2)")

# DC analysis
print("\n[Step 2] PEEC DC analysis...")

builder_dc = PEECBuilder()

for seg in segments:
    start, end = seg
    start_arr = np.array(start)
    end_arr = np.array(end)

    center = (start_arr + end_arr) / 2.0
    direction_vec = end_arr - start_arr
    length = np.linalg.norm(direction_vec)
    direction = direction_vec / length

    builder_dc.add_segment(
        center=center.tolist(),
        direction=direction.tolist(),
        length=length,
        width=equiv_side,
        height=equiv_side,
        sigma=sigma,
        cross_section_type=1  # 1 = circular cross-section
    )

builder_dc.set_frequency(0.0)  # DC
L, R_dc, _, _ = builder_dc.build(include_star=False)

# Extract total L and R (single loop)
L_peec_dc = np.sum(L)  # Total self-inductance
R_peec_dc = np.sum(R_dc)  # Total DC resistance

print(f"\nPEEC DC results:")
print(f"  Inductance: {L_peec_dc * 1e6:.4f} uH")
print(f"  DC resistance: {R_peec_dc * 1e3:.4f} mOhm")

print(f"\nDC comparison:")
print(f"  L error: {abs(L_peec_dc - L_analytical) / L_analytical * 100:.2f}%")
print(f"  R_dc error: {abs(R_peec_dc - R_dc_analytical) / R_dc_analytical * 100:.2f}%")

# ============================================================================
# AC Frequency Sweep
# ============================================================================

print("\n" + "=" * 70)
print("Frequency Sweep: PEEC SIBC vs DC Baseline")
print("=" * 70)

frequencies = np.logspace(1, 6, 20)  # 10 Hz to 1 MHz
R_peec_list = []
L_peec_list = []

print(f"{'Frequency':>12} {'R (mOhm)':>12} {'L (uH)':>12} {'R/R_dc':>12}")
print("-" * 60)

for freq in frequencies:
    builder = PEECBuilder()

    for seg in segments:
        start, end = seg
        start_arr = np.array(start)
        end_arr = np.array(end)

        center = (start_arr + end_arr) / 2.0
        direction_vec = end_arr - start_arr
        length = np.linalg.norm(direction_vec)
        direction = direction_vec / length

        builder.add_segment(
            center=center.tolist(),
            direction=direction.tolist(),
            length=length,
            width=equiv_side,
            height=equiv_side,
            sigma=sigma,
            cross_section_type=1  # 1 = circular cross-section
        )

    builder.set_frequency(freq)
    L, R, _, _ = builder.build(include_star=False)

    R_peec = np.sum(R)
    L_peec = np.sum(L)

    R_peec_list.append(R_peec)
    L_peec_list.append(L_peec)

    print(f"{freq:>12.0f} {R_peec*1e3:>12.4f} {L_peec*1e6:>12.4f} {R_peec/R_peec_dc:>12.2f}")

# ============================================================================
# Visualization
# ============================================================================

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Resistance vs Frequency
ax1.loglog(frequencies, np.array(R_peec_list) * 1e3, 'b-o', label='PEEC with SIBC', linewidth=2, markersize=4)
ax1.axhline(R_peec_dc * 1e3, color='r', linestyle='--', label=f'DC resistance ({R_peec_dc*1e3:.4f} mOhm)')
ax1.set_xlabel('Frequency [Hz]')
ax1.set_ylabel('Resistance [mOhm]')
ax1.set_title('AC Resistance vs Frequency')
ax1.grid(True, which='both', alpha=0.3)
ax1.legend()

# Plot 2: Inductance vs Frequency
ax2.semilogx(frequencies, np.array(L_peec_list) * 1e6, 'g-o', label='PEEC', linewidth=2, markersize=4)
ax2.axhline(L_analytical * 1e6, color='r', linestyle='--', label=f'Analytical ({L_analytical*1e6:.4f} uH)')
ax2.set_xlabel('Frequency [Hz]')
ax2.set_ylabel('Inductance [uH]')
ax2.set_title('Inductance vs Frequency')
ax2.grid(True, which='both', alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.savefig('validation_circular_coil_sibc.png', dpi=300, bbox_inches='tight')
print(f"\nPlot saved: validation_circular_coil_sibc.png")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 70)
print("Validation Summary")
print("=" * 70)

print(f"\nDC validation:")
print(f"  Inductance error: {abs(L_peec_dc - L_analytical) / L_analytical * 100:.2f}%")
print(f"  Resistance error: {abs(R_peec_dc - R_dc_analytical) / R_dc_analytical * 100:.2f}%")

print(f"\nAC resistance (at 1 MHz):")
print(f"  R_dc: {R_peec_dc * 1e3:.4f} mOhm")
print(f"  R_ac: {R_peec_list[-1] * 1e3:.4f} mOhm")
print(f"  Ratio: {R_peec_list[-1] / R_peec_dc:.2f}x")

if abs(L_peec_dc - L_analytical) / L_analytical * 100 < 5.0:
    print("\nOK Inductance validation passed (error < 5%)")
else:
    print("\nWARNING: Inductance error exceeds 5%")

if abs(R_peec_dc - R_dc_analytical) / R_dc_analytical * 100 < 5.0:
    print("OK DC resistance validation passed (error < 5%)")
else:
    print("WARNING: DC resistance error exceeds 5%")

print("\n" + "=" * 70)
print("Validation complete!")
print("=" * 70)
