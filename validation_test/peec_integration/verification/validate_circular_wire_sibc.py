"""
PEEC SIBC Validation: Circular Wire vs Bessel Function Analytical Solution

Compares Radia PEEC SIBC results with exact Bessel function solution
for circular cross-section conductors.

Reference:
  F.W. Grover, "Inductance Calculations", Dover, 1946
  E.B. Rosa, F.W. Grover, "Formulas and Tables for Calculation of
  Mutual and Self-Inductance", NBS Bulletin, 1916

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
from scipy.special import iv  # Modified Bessel functions
import matplotlib.pyplot as plt

from radia.peec_mesh_import import GMSHCenterlineReader

print("=" * 70)
print("PEEC SIBC Validation: Circular Wire vs Bessel Analytical")
print("=" * 70)

# Import PEEC module
try:
    from peec_matrices import PEECBuilder
    print("\nUsing C++ PEEC implementation with SIBC")
except ImportError as e:
    print(f"\nERROR: PEEC module not available: {e}")
    print("Build with: powershell.exe -ExecutionPolicy Bypass -File Build.ps1")
    sys.exit(1)

# Physical constants
SIGMA_COPPER = 5.8e7  # S/m
MU_0 = 4 * np.pi * 1e-7  # H/m

# ============================================================================
# Read Wire Parameters and Mesh from Cubit
# ============================================================================

print("\n[Step 1] Reading Cubit-generated mesh...")

# Read wire parameters
param_file = os.path.join(os.path.dirname(__file__),
                          "cubit_mesh_generation/circular_wire_params.txt")
mesh_file = os.path.join(os.path.dirname(__file__),
                         "cubit_mesh_generation/circular_wire_centerline.msh")

if not os.path.exists(param_file):
    print(f"\nERROR: Parameter file not found: {param_file}")
    print("Run generate_circular_wire.py in Cubit first!")
    sys.exit(1)

if not os.path.exists(mesh_file):
    print(f"\nERROR: Mesh file not found: {mesh_file}")
    print("Run generate_circular_wire.py in Cubit first!")
    sys.exit(1)

# Parse parameters
wire_params = {}
with open(param_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('#') or not line:
            continue
        parts = line.split()
        if len(parts) >= 2:
            key = parts[0]
            value = float(parts[1])
            wire_params[key] = value

wire_radius = wire_params['radius']
wire_length = wire_params['length']
sigma = wire_params.get('sigma', SIGMA_COPPER)
n_segments = int(wire_params.get('n_segments', 10))
cross_section_area = np.pi * wire_radius**2

# Read GMSH centerline mesh
reader = GMSHCenterlineReader(mesh_file)
nodes, edges = reader.read()
segments = reader.get_segments()

print(f"\nWire parameters (from Cubit):")
print(f"  Length: {wire_length * 1e3:.0f} mm")
print(f"  Radius: {wire_radius * 1e3:.2f} mm")
print(f"  Area: {cross_section_area * 1e6:.4f} mm^2")
print(f"  Conductivity: {sigma:.2e} S/m (copper)")
print(f"  Segments: {len(segments)}")

# Verify total length
total_length = reader.get_total_length()
length_error = abs(total_length - wire_length) / wire_length * 100
print(f"\nMesh validation:")
print(f"  Expected length: {wire_length:.6f} m")
print(f"  Actual length: {total_length:.6f} m")
print(f"  Error: {length_error:.4f}%")

if length_error > 0.1:
    print(f"\nWARNING: Mesh length error exceeds 0.1%!")
    sys.exit(1)

# DC resistance (theoretical)
R_dc_theory = wire_length / (sigma * cross_section_area)
print(f"  DC resistance (theory): {R_dc_theory * 1e3:.4f} mOhm")

print(f"\nEquivalent square cross-section:")
equiv_side = np.sqrt(np.pi) * wire_radius
print(f"  Side length: {equiv_side * 1e3:.3f} mm x {equiv_side * 1e3:.3f} mm")
print(f"  (Area matches circular: {cross_section_area * 1e6:.4f} mm^2)")

# ============================================================================
# Analytical Solution: Bessel Function (Exact for Circular Wire)
# ============================================================================

def bessel_impedance_circular(frequency, radius, length, sigma, mu_r=1.0):
    """
    Exact impedance for circular wire using Bessel functions.

    Formula: Z = (k*length)/(2*pi*r*sigma) * I0(kr)/I1(kr)

    where:
      k = sqrt(j*omega*mu*sigma)
      I0, I1 = Modified Bessel functions of first kind, order 0 and 1

    Args:
      frequency: Frequency [Hz]
      radius: Wire radius [m]
      length: Wire length [m]
      sigma: Conductivity [S/m]
      mu_r: Relative permeability

    Returns:
      Complex impedance Z = R + jX [Ohm]
    """
    if frequency == 0:
        # DC resistance
        area = np.pi * radius**2
        return length / (sigma * area)

    omega = 2 * np.pi * frequency
    mu = mu_r * MU_0

    # Propagation constant
    k = np.sqrt(1j * omega * mu * sigma)
    kr = k * radius

    # Modified Bessel functions I0(kr), I1(kr)
    I0 = iv(0, kr)
    I1 = iv(1, kr)

    # Impedance [Ohm]
    # Z = (k*length)/(2*pi*r*sigma) * I0(kr)/I1(kr)
    Z = (k * length) / (2 * np.pi * radius * sigma) * (I0 / I1)

    return Z

# ============================================================================
# PEEC SIBC Calculation
# ============================================================================

frequencies = np.logspace(1, 6, 30)  # 10 Hz to 1 MHz
R_peec = []
X_peec = []
R_analytical = []
X_analytical = []

print("\n" + "=" * 70)
print("Frequency Sweep: PEEC SIBC vs Bessel Analytical")
print("=" * 70)
print(f"{'Frequency':>12} {'R_PEEC':>12} {'R_Bessel':>12} {'Error (%)':>12}")
print("-" * 60)

for freq in frequencies:
    # PEEC with SIBC (approximate circular as square)
    builder = PEECBuilder()

    # Note: PEEC uses rectangular SIBC, so approximate circle as square
    # Equivalent square: side = sqrt(pi) * radius ~ 1.77 * radius
    equiv_side = np.sqrt(np.pi) * wire_radius

    # Add each segment from Cubit mesh
    for seg in segments:
        start, end = seg
        start_arr = np.array(start)
        end_arr = np.array(end)

        # Compute center, direction, length for add_segment API
        center = (start_arr + end_arr) / 2.0
        direction_vec = end_arr - start_arr
        length = np.linalg.norm(direction_vec)
        direction = direction_vec / length  # Normalize

        builder.add_segment(
            center=center.tolist(),
            direction=direction.tolist(),
            length=length,
            width=equiv_side,
            height=equiv_side,
            sigma=sigma
        )

    builder.set_frequency(freq)
    L, R, _, _ = builder.build(include_star=False)

    # Extract R and L (sum all segments - series connection)
    R_peec_val = np.sum(R)  # Total resistance = sum of all segments
    L_peec_val = np.sum(L)  # Total inductance = sum of all self/mutual inductances
    X_peec_val = 2 * np.pi * freq * L_peec_val if freq > 0 else 0

    # Analytical (Bessel)
    Z_analytical = bessel_impedance_circular(freq, wire_radius, wire_length, sigma)
    R_analytical_val = Z_analytical.real
    X_analytical_val = Z_analytical.imag

    R_peec.append(R_peec_val)
    X_peec.append(X_peec_val)
    R_analytical.append(R_analytical_val)
    X_analytical.append(X_analytical_val)

    error_pct = abs(R_peec_val - R_analytical_val) / abs(R_analytical_val) * 100 if R_analytical_val != 0 else 0

    print(f"{freq:>12.0f} {R_peec_val:>12.6f} {R_analytical_val:>12.6f} {error_pct:>11.2f}%")

R_peec = np.array(R_peec)
X_peec = np.array(X_peec)
R_analytical = np.array(R_analytical)
X_analytical = np.array(X_analytical)

# ============================================================================
# Error Analysis
# ============================================================================

print("\n" + "=" * 70)
print("Error Analysis")
print("=" * 70)

# Resistance error
R_error_pct = np.abs(R_peec - R_analytical) / np.abs(R_analytical) * 100
print(f"\nResistance Error:")
print(f"  Mean error: {np.mean(R_error_pct):.2f}%")
print(f"  Max error: {np.max(R_error_pct):.2f}%")
print(f"  RMS error: {np.sqrt(np.mean(R_error_pct**2)):.2f}%")

# Reactance error
X_error_pct = np.abs(X_peec - X_analytical) / np.abs(X_analytical) * 100
X_error_pct = X_error_pct[frequencies > 100]  # Ignore low freq where X is small
print(f"\nReactance Error (f > 100 Hz):")
print(f"  Mean error: {np.mean(X_error_pct):.2f}%")
print(f"  Max error: {np.max(X_error_pct):.2f}%")

# ============================================================================
# Visualization
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Resistance vs Frequency
ax1 = axes[0, 0]
ax1.loglog(frequencies, R_peec, 'b-o', label='PEEC SIBC', linewidth=2, markersize=4)
ax1.loglog(frequencies, R_analytical, 'r--', label='Bessel Analytical', linewidth=2)
ax1.set_xlabel('Frequency [Hz]')
ax1.set_ylabel('Resistance [Ohm]')
ax1.set_title('Resistance vs Frequency')
ax1.grid(True, which='both', alpha=0.3)
ax1.legend()

# Plot 2: Reactance vs Frequency
ax2 = axes[0, 1]
ax2.loglog(frequencies, np.abs(X_peec), 'b-o', label='PEEC SIBC', linewidth=2, markersize=4)
ax2.loglog(frequencies, np.abs(X_analytical), 'r--', label='Bessel Analytical', linewidth=2)
ax2.set_xlabel('Frequency [Hz]')
ax2.set_ylabel('Reactance [Ohm]')
ax2.set_title('Reactance vs Frequency')
ax2.grid(True, which='both', alpha=0.3)
ax2.legend()

# Plot 3: Resistance Error
ax3 = axes[1, 0]
ax3.semilogx(frequencies, R_error_pct, 'g-o', linewidth=2, markersize=4)
ax3.set_xlabel('Frequency [Hz]')
ax3.set_ylabel('Error [%]')
ax3.set_title('Resistance Error vs Analytical')
ax3.grid(True, which='both', alpha=0.3)
ax3.axhline(5, color='r', linestyle='--', label='5% threshold')
ax3.axhline(10, color='orange', linestyle='--', label='10% threshold')
ax3.legend()

# Plot 4: Impedance Magnitude
ax4 = axes[1, 1]
Z_peec = np.sqrt(R_peec**2 + X_peec**2)
Z_analytical_mag = np.sqrt(R_analytical**2 + X_analytical**2)
ax4.loglog(frequencies, Z_peec, 'b-o', label='PEEC SIBC', linewidth=2, markersize=4)
ax4.loglog(frequencies, Z_analytical_mag, 'r--', label='Bessel Analytical', linewidth=2)
ax4.set_xlabel('Frequency [Hz]')
ax4.set_ylabel('Impedance Magnitude |Z| [Ohm]')
ax4.set_title('Impedance Magnitude vs Frequency')
ax4.grid(True, which='both', alpha=0.3)
ax4.legend()

plt.tight_layout()
plt.savefig('validation_circular_wire_sibc.png', dpi=300, bbox_inches='tight')
print(f"\nPlot saved: validation_circular_wire_sibc.png")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 70)
print("Validation Summary")
print("=" * 70)

print("\nPEEC SIBC Implementation:")
print("  - Rectangular cross-section approximation")
print(f"  - Equivalent square: {equiv_side * 1e3:.3f} mm x {equiv_side * 1e3:.3f} mm")
print(f"  - (Area matches circular: {cross_section_area * 1e6:.4f} mm^2)")

print("\nBessel Analytical Solution:")
print("  - Exact for circular cross-section")
print("  - I0(kr)/I1(kr) impedance formula")

print("\nError Summary:")
print(f"  - Resistance mean error: {np.mean(R_error_pct):.2f}%")
print(f"  - Resistance max error: {np.max(R_error_pct):.2f}%")

if np.max(R_error_pct) < 10:
    print("\nOK SIBC implementation validated (error < 10%)")
    print("   Rectangular approximation is acceptable for circular wires")
else:
    print("\nWARNING: Error exceeds 10% at some frequencies")
    print("         Consider implementing circular SIBC (ComputeCircular)")

print("\nNote:")
print("  - Error comes from rectangular approximation of circular cross-section")
print("  - For exact circular wire results, use ComputeCircular() method")
print("  - (Currently implemented but not exposed to Python API)")

print("\n" + "=" * 70)
print("Validation complete!")
print("=" * 70)
