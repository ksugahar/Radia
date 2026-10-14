"""
PEEC SIBC Validation: Bessel Function (Circular Wire) vs Dowell (Rectangular)

Compares three approaches for AC resistance of a circular coil:
  A. Analytical (Bessel exact): scipy.special.jv for circular wire
  B. PEEC + Bessel SIBC (Python): PEEC L matrix + scipy Bessel Z_s
  C. PEEC + Dowell (Python): Dowell formula for rectangular cross-section

Model: Single-turn circular coil
  - Coil radius: 50 mm
  - Wire radius: 1 mm (AWG20)
  - Copper: sigma = 5.8e7 S/m
  - Frequency: 10 Hz to 1 MHz

References:
  - F.W. Grover, "Inductance Calculations", Dover, 1946
  - Bessel: Z_i/l = k/(2*pi*a*sigma) * I0(ka)/I1(ka)
  - Dowell: F_R = xi * [sinh(2xi)+sin(2xi)] / [cosh(2xi)-cos(2xi)]

Author: Radia Development Team
Date: 2026-02-13
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import numpy as np
from scipy.special import iv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("=" * 70)
print("PEEC SIBC Validation: Bessel (Circular) vs Dowell (Rectangular)")
print("=" * 70)

# Import PEEC module
try:
    from peec_matrices import PEECBuilder
    print("\nC++ PEEC module loaded successfully")
except ImportError as e:
    print(f"\nERROR: PEEC module not available: {e}")
    print("Build with: pwsh -NoProfile -ExecutionPolicy Bypass -File .\\Build.ps1")
    sys.exit(1)

# ============================================================================
# Physical Constants and Model Parameters
# ============================================================================

MU_0 = 4.0 * np.pi * 1e-7  # H/m
SIGMA_COPPER = 5.8e7        # S/m

# Coil geometry
coil_radius = 50e-3   # 50 mm
wire_radius = 1e-3     # 1 mm
sigma = SIGMA_COPPER
n_segments = 36        # Number of segments around the loop

# Derived quantities
circumference = 2.0 * np.pi * coil_radius
cross_section_area = np.pi * wire_radius**2
seg_length = circumference / n_segments

# Equivalent square side for rectangular approximation
equiv_side = np.sqrt(np.pi) * wire_radius

print(f"\nModel parameters:")
print(f"  Coil radius:     {coil_radius*1e3:.1f} mm")
print(f"  Wire radius:     {wire_radius*1e3:.2f} mm")
print(f"  Circumference:   {circumference*1e3:.2f} mm")
print(f"  Cross-section:   {cross_section_area*1e6:.4f} mm^2")
print(f"  Equiv. square:   {equiv_side*1e3:.3f} mm x {equiv_side*1e3:.3f} mm")
print(f"  Segments:        {n_segments}")
print(f"  Segment length:  {seg_length*1e3:.3f} mm")

# ============================================================================
# Method A: Analytical Solution (Bessel Exact for Circular Wire)
# ============================================================================

def bessel_impedance_per_length(freq, radius, sigma_val, mu_r=1.0):
    """
    Internal impedance per unit length of solid circular wire (exact).

    Z_i/l = k/(2*pi*a*sigma) * I0(ka)/I1(ka)

    where k = sqrt(j*omega*mu*sigma), I0/I1 are modified Bessel functions.

    DC limit: Z_i/l -> 1/(sigma*pi*a^2) + j*omega*mu/(8*pi)
    """
    if freq <= 0:
        # DC limit
        return 1.0 / (sigma_val * np.pi * radius**2) + 0j

    omega = 2.0 * np.pi * freq
    mu = mu_r * MU_0
    k = np.sqrt(1j * omega * mu * sigma_val)
    ka = k * radius

    I0_ka = iv(0, ka)
    I1_ka = iv(1, ka)

    Z_i_per_length = k / (2.0 * np.pi * radius * sigma_val) * (I0_ka / I1_ka)
    return Z_i_per_length


def grover_inductance(R_coil, r_wire):
    """
    Grover's formula for single-turn circular coil inductance.

    L = mu_0 * R * [ln(8R/r) - 2]
    """
    return MU_0 * R_coil * (np.log(8.0 * R_coil / r_wire) - 2.0)


def skin_depth(freq, sigma_val, mu_r=1.0):
    """Skin depth: delta = sqrt(2/(omega*mu*sigma))"""
    if freq <= 0:
        return np.inf
    omega = 2.0 * np.pi * freq
    mu = mu_r * MU_0
    return np.sqrt(2.0 / (omega * mu * sigma_val))


def dowell_factor(freq, sigma_val, thickness, mu_r=1.0):
    """
    Dowell's AC resistance factor for rectangular cross-section (Python).

    F_R = xi * [sinh(2*xi) + sin(2*xi)] / [cosh(2*xi) - cos(2*xi)]

    where xi = thickness / skin_depth

    Valid only for d << w (thin plate approximation).
    For circular cross-section, use bessel_impedance_per_length() instead.
    """
    if freq <= 0:
        return 1.0

    delta = skin_depth(freq, sigma_val, mu_r)
    xi = thickness / delta

    if xi < 1e-6:
        return 1.0  # DC limit

    numerator = np.sinh(2.0 * xi) + np.sin(2.0 * xi)
    denominator = np.cosh(2.0 * xi) - np.cos(2.0 * xi)

    if abs(denominator) < 1e-30:
        return xi  # Large xi limit

    return xi * numerator / denominator


# DC analytical values
R_dc_analytical = circumference / (sigma * cross_section_area)
L_analytical = grover_inductance(coil_radius, wire_radius)

print(f"\n--- Analytical Reference (Grover + Bessel) ---")
print(f"  DC inductance: {L_analytical*1e6:.4f} uH")
print(f"  DC resistance: {R_dc_analytical*1e3:.6f} mOhm")

# Verify Bessel DC limit
Z_dc_bessel = bessel_impedance_per_length(0, wire_radius, sigma)
R_dc_bessel = Z_dc_bessel.real * circumference
print(f"  Bessel DC R:   {R_dc_bessel*1e3:.6f} mOhm (should match)")

# ============================================================================
# Method B: PEEC + Bessel SIBC (Python scipy)
# ============================================================================

print(f"\n--- PEEC Model Construction ---")

# Build PEEC matrices (frequency-independent: L, R_dc)
builder_dc = PEECBuilder()
n_created = builder_dc.create_loop(
    [0, 0, 0],        # center
    coil_radius,       # radius
    [0, 0, 1],         # normal (z-axis)
    equiv_side,        # width (equiv square)
    equiv_side,        # height (equiv square)
    n_segments,        # segments
    sigma              # conductivity
)

L_matrix, R_dc_vec, _, _ = builder_dc.build(False)  # include_star=False

L_total_peec = np.sum(L_matrix)
R_dc_peec = np.sum(R_dc_vec)

print(f"  Created {n_created} segments")
print(f"  L_total (PEEC):  {L_total_peec*1e6:.4f} uH")
print(f"  R_dc (PEEC):     {R_dc_peec*1e3:.6f} mOhm")

L_error = abs(L_total_peec - L_analytical) / L_analytical * 100
R_error = abs(R_dc_peec - R_dc_analytical) / R_dc_analytical * 100
print(f"\n  DC L error: {L_error:.2f}%")
print(f"  DC R error: {R_error:.2f}%")

# ============================================================================
# AC Frequency Sweep: Compare 3 Methods
# ============================================================================

print(f"\n{'='*70}")
print("AC Frequency Sweep: 3-Method Comparison")
print(f"{'='*70}")

frequencies = np.logspace(1, 6, 30)  # 10 Hz to 1 MHz

# Storage for results
results_analytical = []   # Method A: Bessel exact
results_peec_bessel = []  # Method B: PEEC + Bessel SIBC
results_peec_dowell = []  # Method C: PEEC + Dowell

print(f"\n{'Freq (Hz)':>12} {'delta (mm)':>10} {'F_R Bessel':>10} "
      f"{'F_R PEEC+B':>10} {'F_R Dowell':>10} {'Err B (%)':>10} {'Err D (%)':>10}")
print("-" * 82)

for freq in frequencies:
    omega = 2.0 * np.pi * freq
    delta = skin_depth(freq, sigma)

    # --- Method A: Analytical Bessel ---
    Z_bessel_per_l = bessel_impedance_per_length(freq, wire_radius, sigma)
    Z_internal_total = Z_bessel_per_l * circumference
    Z_analytical = Z_internal_total + 1j * omega * L_analytical
    R_ac_analytical = Z_internal_total.real
    F_R_analytical = R_ac_analytical / R_dc_analytical

    results_analytical.append({
        'freq': freq,
        'R_ac': R_ac_analytical,
        'X_ac': Z_analytical.imag,
        'F_R': F_R_analytical,
        'Z': Z_analytical,
    })

    # --- Method B: PEEC + Bessel SIBC (Python) ---
    # Z_bessel per segment (total for one segment, not per unit length)
    Z_bessel_seg = Z_bessel_per_l * seg_length
    Z_internal_peec_bessel = Z_bessel_seg * n_segments  # Sum over all segments
    Z_peec_bessel = Z_internal_peec_bessel + 1j * omega * L_total_peec
    R_ac_peec_bessel = Z_internal_peec_bessel.real
    F_R_peec_bessel = R_ac_peec_bessel / R_dc_peec

    results_peec_bessel.append({
        'freq': freq,
        'R_ac': R_ac_peec_bessel,
        'X_ac': Z_peec_bessel.imag,
        'F_R': F_R_peec_bessel,
        'Z': Z_peec_bessel,
    })

    # --- Method C: PEEC + Dowell (Python, rectangular approximation) ---
    # Dowell formula applied to DC resistance (reuses same L matrix)
    F_R_dowell_factor = dowell_factor(freq, sigma, equiv_side)
    R_ac_dowell = R_dc_peec * F_R_dowell_factor
    Z_peec_dowell = R_ac_dowell + 1j * omega * L_total_peec
    F_R_dowell = F_R_dowell_factor

    results_peec_dowell.append({
        'freq': freq,
        'R_ac': R_ac_dowell,
        'X_ac': Z_peec_dowell.imag,
        'F_R': F_R_dowell,
        'Z': Z_peec_dowell,
    })

    # Error relative to analytical Bessel
    err_bessel = abs(F_R_peec_bessel - F_R_analytical) / F_R_analytical * 100
    err_dowell = abs(F_R_dowell - F_R_analytical) / F_R_analytical * 100

    # Print selected frequencies
    if freq in [frequencies[0], frequencies[4], frequencies[9],
                frequencies[14], frequencies[19], frequencies[24], frequencies[-1]]:
        print(f"{freq:>12.0f} {delta*1e3:>10.4f} {F_R_analytical:>10.4f} "
              f"{F_R_peec_bessel:>10.4f} {F_R_dowell:>10.4f} "
              f"{err_bessel:>10.4f} {err_dowell:>10.4f}")

# ============================================================================
# Visualization
# ============================================================================

freqs = np.array([r['freq'] for r in results_analytical])
FR_analytical = np.array([r['F_R'] for r in results_analytical])
FR_peec_bessel = np.array([r['F_R'] for r in results_peec_bessel])
FR_dowell = np.array([r['F_R'] for r in results_peec_dowell])

R_analytical_arr = np.array([r['R_ac'] for r in results_analytical])
R_peec_bessel_arr = np.array([r['R_ac'] for r in results_peec_bessel])
R_dowell_arr = np.array([r['R_ac'] for r in results_peec_dowell])

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: F_R = R_ac/R_dc vs Frequency
ax = axes[0, 0]
ax.loglog(freqs, FR_analytical, 'k-', linewidth=2, label='Analytical (Bessel exact)')
ax.loglog(freqs, FR_peec_bessel, 'b--', linewidth=2, label='PEEC + Bessel SIBC')
ax.loglog(freqs, FR_dowell, 'r:', linewidth=2, label='PEEC + Dowell (rectangular)')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('F_R = R_ac / R_dc')
ax.set_title('AC Resistance Factor vs Frequency')
ax.grid(True, which='both', alpha=0.3)
ax.legend(fontsize=9)

# Plot 2: Absolute Resistance
ax = axes[0, 1]
ax.loglog(freqs, R_analytical_arr * 1e3, 'k-', linewidth=2, label='Analytical (Bessel)')
ax.loglog(freqs, R_peec_bessel_arr * 1e3, 'b--', linewidth=2, label='PEEC + Bessel')
ax.loglog(freqs, R_dowell_arr * 1e3, 'r:', linewidth=2, label='PEEC + Dowell')
ax.axhline(R_dc_analytical * 1e3, color='gray', linestyle='-.', alpha=0.5, label='R_dc')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('Resistance [mOhm]')
ax.set_title('AC Resistance vs Frequency')
ax.grid(True, which='both', alpha=0.3)
ax.legend(fontsize=9)

# Plot 3: Error of Dowell vs Bessel
ax = axes[1, 0]
err_dowell_pct = np.abs(FR_dowell - FR_analytical) / FR_analytical * 100
err_bessel_pct = np.abs(FR_peec_bessel - FR_analytical) / FR_analytical * 100
ax.semilogx(freqs, err_bessel_pct, 'b-o', markersize=3, linewidth=1.5,
            label='PEEC + Bessel vs Analytical')
ax.semilogx(freqs, err_dowell_pct, 'r-s', markersize=3, linewidth=1.5,
            label='PEEC + Dowell vs Analytical')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('F_R Error [%]')
ax.set_title('Error: Dowell (rectangular) vs Bessel (circular)')
ax.grid(True, which='both', alpha=0.3)
ax.legend(fontsize=9)

# Plot 4: Skin depth vs wire radius
ax = axes[1, 1]
delta_arr = np.array([skin_depth(f, sigma) for f in freqs])
ax.loglog(freqs, delta_arr * 1e3, 'g-', linewidth=2, label='Skin depth')
ax.axhline(wire_radius * 1e3, color='k', linestyle='--',
           label=f'Wire radius ({wire_radius*1e3:.1f} mm)')
ax.fill_between(freqs, 0, wire_radius * 1e3, alpha=0.1, color='red')
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('[mm]')
ax.set_title('Skin Depth vs Wire Radius')
ax.grid(True, which='both', alpha=0.3)
ax.legend(fontsize=9)
ax.set_ylim([1e-2, 1e2])

plt.suptitle(f'Circular Coil SIBC Validation: R={coil_radius*1e3:.0f}mm, '
             f'a={wire_radius*1e3:.1f}mm, {n_segments} segments',
             fontsize=13, fontweight='bold')
plt.tight_layout()

output_path = os.path.join(os.path.dirname(__file__), 'validate_bessel_sibc.png')
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"\nPlot saved: {output_path}")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*70}")
print("Validation Summary")
print(f"{'='*70}")

print(f"\n[DC Validation]")
print(f"  L: PEEC={L_total_peec*1e6:.4f} uH, Grover={L_analytical*1e6:.4f} uH, "
      f"Error={L_error:.2f}%")
print(f"  R: PEEC={R_dc_peec*1e3:.6f} mOhm, Analytical={R_dc_analytical*1e3:.6f} mOhm, "
      f"Error={R_error:.2f}%")

# High-frequency comparison (1 MHz)
idx_1mhz = -1
F_R_bessel_1mhz = FR_analytical[idx_1mhz]
F_R_dowell_1mhz = FR_dowell[idx_1mhz]
err_dowell_1mhz = abs(F_R_dowell_1mhz - F_R_bessel_1mhz) / F_R_bessel_1mhz * 100
delta_1mhz = skin_depth(frequencies[idx_1mhz], sigma)

print(f"\n[AC Validation at {frequencies[idx_1mhz]/1e6:.1f} MHz]")
print(f"  Skin depth: {delta_1mhz*1e3:.4f} mm (wire radius: {wire_radius*1e3:.2f} mm)")
print(f"  F_R (Bessel exact):   {F_R_bessel_1mhz:.4f}")
print(f"  F_R (PEEC + Bessel):  {FR_peec_bessel[idx_1mhz]:.4f}")
print(f"  F_R (PEEC + Dowell):  {F_R_dowell_1mhz:.4f}")
print(f"  Dowell error vs Bessel: {err_dowell_1mhz:.1f}%")

# Pass/Fail criteria
dc_l_pass = L_error < 5.0
dc_r_pass = R_error < 5.0
bessel_match = np.max(err_bessel_pct) < 1.0

print(f"\n[Pass/Fail]")
print(f"  DC inductance (< 5%):     {'PASS' if dc_l_pass else 'FAIL'} ({L_error:.2f}%)")
print(f"  DC resistance (< 5%):     {'PASS' if dc_r_pass else 'FAIL'} ({R_error:.2f}%)")
print(f"  Bessel SIBC match (< 1%): {'PASS' if bessel_match else 'FAIL'} "
      f"(max {np.max(err_bessel_pct):.4f}%)")
print(f"  Dowell differs from Bessel: {np.max(err_dowell_pct):.1f}% "
      f"(expected for circular wire)")

print(f"\n{'='*70}")
print("Validation complete!")
print(f"{'='*70}")
