"""
NMR Magnet Homogeneity Analysis with Spherical Harmonics

This demo shows NMR magnet field analysis:
1. Field homogeneity metrics (ppm)
2. Spherical harmonic expansion
3. Shim coil field patterns
4. Homogeneity optimization

Physical setup:
- Main magnet: Cylindrical permanent magnet assembly
- Target region: Spherical volume (DSV - Diameter of Spherical Volume)

Author: Radia Development Team
Date: 2026-01-16
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lpmv  # Associated Legendre polynomials


# Physical constants
MU_0 = 4 * np.pi * 1e-7  # [H/m]


def spherical_to_cartesian(r, theta, phi):
    """Convert spherical (r, theta, phi) to Cartesian (x, y, z)."""
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z


def cartesian_to_spherical(x, y, z):
    """Convert Cartesian (x, y, z) to spherical (r, theta, phi)."""
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r) if r > 0 else 0
    phi = np.arctan2(y, x)
    return r, theta, phi


def spherical_harmonic_real(l, m, theta, phi):
    """
    Real spherical harmonic Y_l^m(theta, phi).

    Uses convention:
    - Y_l^0 = P_l^0(cos(theta))
    - Y_l^m = P_l^m(cos(theta)) * cos(m*phi)  for m > 0
    - Y_l^(-m) = P_l^m(cos(theta)) * sin(m*phi)  for m > 0
    """
    cos_theta = np.cos(theta)
    abs_m = abs(m)

    # Associated Legendre polynomial
    P_lm = lpmv(abs_m, l, cos_theta)

    if m == 0:
        return P_lm
    elif m > 0:
        return P_lm * np.cos(m * phi)
    else:
        return P_lm * np.sin(abs_m * phi)


def expand_field_spherical_harmonics(Bz_samples, r_samples, theta_samples, phi_samples,
                                       l_max=6):
    """
    Expand Bz field in spherical harmonics.

    Bz(r, theta, phi) = sum_l sum_m C_lm * (r/R)^l * Y_l^m(theta, phi)

    Parameters:
        Bz_samples: Bz values at sample points [T]
        r_samples, theta_samples, phi_samples: Spherical coordinates of samples
        l_max: Maximum order of expansion

    Returns:
        coefficients: Dict of C_lm values [T]
    """
    n_samples = len(Bz_samples)

    # Reference radius (maximum sampling radius)
    R = np.max(r_samples)

    # Build design matrix
    n_coeffs = (l_max + 1)**2
    A = np.zeros((n_samples, n_coeffs))

    col = 0
    for l in range(l_max + 1):
        for m in range(-l, l + 1):
            for i in range(n_samples):
                r_norm = r_samples[i] / R if R > 0 else 0
                Y_lm = spherical_harmonic_real(l, m, theta_samples[i], phi_samples[i])
                A[i, col] = (r_norm ** l) * Y_lm
            col += 1

    # Least squares fit
    coeffs, residuals, rank, s = np.linalg.lstsq(A, Bz_samples, rcond=None)

    # Organize coefficients
    result = {}
    col = 0
    for l in range(l_max + 1):
        for m in range(-l, l + 1):
            result[(l, m)] = coeffs[col]
            col += 1

    return result


def calc_homogeneity_ppm(Bz_samples, B0=None):
    """
    Calculate field homogeneity in ppm.

    Parameters:
        Bz_samples: Bz values at sample points [T]
        B0: Reference field (default: mean of samples) [T]

    Returns:
        dict with peak-to-peak and RMS homogeneity in ppm
    """
    if B0 is None:
        B0 = np.mean(Bz_samples)

    dB = Bz_samples - B0

    # Peak-to-peak homogeneity
    pp_ppm = (np.max(dB) - np.min(dB)) / np.abs(B0) * 1e6

    # RMS homogeneity
    rms_ppm = np.std(dB) / np.abs(B0) * 1e6

    return {
        'B0': B0,
        'peak_to_peak_ppm': pp_ppm,
        'rms_ppm': rms_ppm,
        'max_deviation_ppm': np.max(np.abs(dB)) / np.abs(B0) * 1e6
    }


def generate_dsv_samples(R_dsv, n_radial=5, n_theta=10, n_phi=20):
    """
    Generate sample points within a spherical volume.

    Parameters:
        R_dsv: Radius of DSV [m]
        n_radial: Number of radial points
        n_theta: Number of theta points
        n_phi: Number of phi points

    Returns:
        r, theta, phi arrays and x, y, z arrays
    """
    r_vals = np.linspace(0, R_dsv, n_radial)
    theta_vals = np.linspace(0.1, np.pi - 0.1, n_theta)  # Avoid poles
    phi_vals = np.linspace(0, 2*np.pi, n_phi, endpoint=False)

    # Generate all combinations
    r_list, theta_list, phi_list = [], [], []
    x_list, y_list, z_list = [], [], []

    for r in r_vals[1:]:  # Skip r=0
        for theta in theta_vals:
            for phi in phi_vals:
                r_list.append(r)
                theta_list.append(theta)
                phi_list.append(phi)
                x, y, z = spherical_to_cartesian(r, theta, phi)
                x_list.append(x)
                y_list.append(y)
                z_list.append(z)

    return (np.array(r_list), np.array(theta_list), np.array(phi_list),
            np.array(x_list), np.array(y_list), np.array(z_list))


def ideal_dipole_field(x, y, z, M, direction='z'):
    """
    Ideal magnetic dipole field (far-field approximation).

    Parameters:
        x, y, z: Observation points [m]
        M: Dipole moment [A*m^2]
        direction: Dipole direction ('x', 'y', or 'z')

    Returns:
        Bx, By, Bz components [T]
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    r3 = r**3
    r5 = r**5

    # Dipole moment vector
    if direction == 'z':
        mx, my, mz = 0, 0, M
    elif direction == 'x':
        mx, my, mz = M, 0, 0
    elif direction == 'y':
        mx, my, mz = 0, M, 0

    # m dot r
    m_dot_r = mx * x + my * y + mz * z

    # B = (mu0/4pi) * (3r(m.r)/r^5 - m/r^3)
    factor = MU_0 / (4 * np.pi)

    Bx = factor * (3 * x * m_dot_r / r5 - mx / r3)
    By = factor * (3 * y * m_dot_r / r5 - my / r3)
    Bz = factor * (3 * z * m_dot_r / r5 - mz / r3)

    return Bx, By, Bz


def shim_coil_field_pattern(l, m, x, y, z, R_shim=0.1):
    """
    Idealized shim coil field pattern for harmonics cancellation.

    A shim coil designed to produce Y_l^m harmonic pattern.

    Parameters:
        l, m: Spherical harmonic order
        x, y, z: Observation points [m]
        R_shim: Shim coil radius [m]

    Returns:
        Bz contribution [T] (normalized)
    """
    r, theta, phi = [], [], []
    for i in range(len(x)):
        ri, ti, pi = cartesian_to_spherical(x[i], y[i], z[i])
        r.append(ri)
        theta.append(ti)
        phi.append(pi)

    r = np.array(r)
    theta = np.array(theta)
    phi = np.array(phi)

    # Normalized field pattern
    Bz = np.zeros(len(x))
    for i in range(len(x)):
        if r[i] > 0:
            Bz[i] = (r[i] / R_shim)**l * spherical_harmonic_real(l, m, theta[i], phi[i])

    return Bz


def main():
    """NMR magnet homogeneity analysis demo."""
    print("=" * 70)
    print("NMR Magnet Homogeneity Analysis")
    print("=" * 70)

    # === System Parameters ===
    # Target field and DSV
    B0_target = 1.5      # Target field [T] (1.5T MRI)
    R_dsv = 0.25         # DSV radius [m] (50 cm diameter)

    print(f"\n--- System Parameters ---")
    print(f"Target B0: {B0_target} T")
    print(f"DSV diameter: {R_dsv*2*100:.0f} cm")

    # === Generate Sample Points ===
    print(f"\n--- Generating DSV Sample Points ---")
    r_samp, theta_samp, phi_samp, x_samp, y_samp, z_samp = generate_dsv_samples(
        R_dsv, n_radial=5, n_theta=10, n_phi=20
    )
    n_points = len(r_samp)
    print(f"Number of sample points: {n_points}")

    # === Simulate Imperfect Magnet Field ===
    # Add realistic inhomogeneities (2nd and 4th order)
    print(f"\n--- Simulating Imperfect Magnet Field ---")

    # Base field plus inhomogeneities
    Bz = np.ones(n_points) * B0_target

    # Add Z2 (l=2, m=0) inhomogeneity: 100 ppm
    C20 = B0_target * 100e-6  # 100 ppm
    for i in range(n_points):
        Bz[i] += C20 * (r_samp[i] / R_dsv)**2 * spherical_harmonic_real(2, 0, theta_samp[i], phi_samp[i])

    # Add Z4 (l=4, m=0) inhomogeneity: 20 ppm
    C40 = B0_target * 20e-6  # 20 ppm
    for i in range(n_points):
        Bz[i] += C40 * (r_samp[i] / R_dsv)**4 * spherical_harmonic_real(4, 0, theta_samp[i], phi_samp[i])

    # Add X (l=1, m=1) gradient: 50 ppm
    C11 = B0_target * 50e-6
    for i in range(n_points):
        Bz[i] += C11 * (r_samp[i] / R_dsv) * spherical_harmonic_real(1, 1, theta_samp[i], phi_samp[i])

    # === Calculate Initial Homogeneity ===
    homo = calc_homogeneity_ppm(Bz, B0_target)
    print(f"\nInitial field homogeneity:")
    print(f"  B0 = {homo['B0']:.6f} T")
    print(f"  Peak-to-peak: {homo['peak_to_peak_ppm']:.1f} ppm")
    print(f"  RMS: {homo['rms_ppm']:.1f} ppm")
    print(f"  Max deviation: {homo['max_deviation_ppm']:.1f} ppm")

    # === Spherical Harmonic Expansion ===
    print(f"\n--- Spherical Harmonic Analysis ---")
    coeffs = expand_field_spherical_harmonics(Bz, r_samp, theta_samp, phi_samp, l_max=6)

    print(f"\nSignificant harmonics (> 1 ppm):")
    print(f"{'(l,m)':>10} {'C_lm [T]':>15} {'[ppm]':>10}")
    print("-" * 40)

    for (l, m), C in sorted(coeffs.items()):
        ppm = np.abs(C) / B0_target * 1e6
        if ppm > 1:
            print(f"{str((l,m)):>10} {C:>15.3e} {ppm:>10.1f}")

    # === Shim Coil Optimization ===
    print(f"\n--- Shim Coil Optimization ---")
    print("Applying shim currents to cancel harmonics...")

    # Calculate shim currents needed
    shim_currents = {}
    harmonics_to_shim = [(1, 1), (2, 0), (4, 0)]  # X, Z2, Z4

    for (l, m) in harmonics_to_shim:
        if (l, m) in coeffs:
            shim_currents[(l, m)] = -coeffs[(l, m)]  # Negative to cancel

    print(f"\nShim currents (arbitrary units):")
    for (l, m), current in shim_currents.items():
        name = ['Z', 'X', 'Y', 'Z2', 'ZX', 'ZY', 'X2-Y2', 'XY'][l + abs(m)] if l <= 2 else f"({l},{m})"
        print(f"  {name}: {current:.3e} T (to cancel {coeffs[(l,m)]/B0_target*1e6:.1f} ppm)")

    # Apply shimming
    Bz_shimmed = Bz.copy()
    for (l, m), current in shim_currents.items():
        for i in range(n_points):
            shim_field = current * (r_samp[i] / R_dsv)**l * spherical_harmonic_real(l, m, theta_samp[i], phi_samp[i])
            Bz_shimmed[i] += shim_field

    # === Calculate Shimmed Homogeneity ===
    homo_shimmed = calc_homogeneity_ppm(Bz_shimmed, B0_target)
    print(f"\nShimmed field homogeneity:")
    print(f"  Peak-to-peak: {homo_shimmed['peak_to_peak_ppm']:.1f} ppm")
    print(f"  RMS: {homo_shimmed['rms_ppm']:.1f} ppm")
    print(f"  Max deviation: {homo_shimmed['max_deviation_ppm']:.1f} ppm")

    improvement = homo['rms_ppm'] / homo_shimmed['rms_ppm']
    print(f"\n  Improvement factor: {improvement:.1f}x")

    # === Plot Results ===
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: Field distribution before shimming
    ax1 = axes[0, 0]
    scatter = ax1.scatter(z_samp*1e3, (Bz - B0_target)/B0_target*1e6,
                          c=np.sqrt(x_samp**2 + y_samp**2)*1e3,
                          cmap='viridis', alpha=0.7)
    ax1.set_xlabel('z [mm]')
    ax1.set_ylabel('dB/B0 [ppm]')
    ax1.set_title('Field Inhomogeneity vs Z (Before Shimming)')
    plt.colorbar(scatter, ax=ax1, label='Radial distance [mm]')
    ax1.grid(True, alpha=0.3)

    # Plot 2: Field distribution after shimming
    ax2 = axes[0, 1]
    scatter = ax2.scatter(z_samp*1e3, (Bz_shimmed - B0_target)/B0_target*1e6,
                          c=np.sqrt(x_samp**2 + y_samp**2)*1e3,
                          cmap='viridis', alpha=0.7)
    ax2.set_xlabel('z [mm]')
    ax2.set_ylabel('dB/B0 [ppm]')
    ax2.set_title('Field Inhomogeneity vs Z (After Shimming)')
    plt.colorbar(scatter, ax=ax2, label='Radial distance [mm]')
    ax2.grid(True, alpha=0.3)

    # Plot 3: Histogram comparison
    ax3 = axes[1, 0]
    bins = np.linspace(-150, 150, 50)
    ax3.hist((Bz - B0_target)/B0_target*1e6, bins=bins, alpha=0.5,
             label=f'Before: {homo["rms_ppm"]:.1f} ppm RMS', density=True)
    ax3.hist((Bz_shimmed - B0_target)/B0_target*1e6, bins=bins, alpha=0.5,
             label=f'After: {homo_shimmed["rms_ppm"]:.1f} ppm RMS', density=True)
    ax3.set_xlabel('dB/B0 [ppm]')
    ax3.set_ylabel('Probability Density')
    ax3.set_title('Field Homogeneity Distribution')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # Plot 4: Spherical harmonic spectrum
    ax4 = axes[1, 1]
    l_vals = []
    ppm_vals = []
    labels = []
    for (l, m), C in sorted(coeffs.items()):
        ppm = np.abs(C) / B0_target * 1e6
        if ppm > 0.1:
            l_vals.append(l + m * 0.1)  # Offset for visibility
            ppm_vals.append(ppm)
            labels.append(f"({l},{m})")

    bars = ax4.bar(range(len(l_vals)), ppm_vals, alpha=0.7)
    ax4.set_xticks(range(len(l_vals)))
    ax4.set_xticklabels(labels, rotation=45)
    ax4.set_xlabel('Harmonic (l, m)')
    ax4.set_ylabel('Amplitude [ppm]')
    ax4.set_title('Spherical Harmonic Spectrum')
    ax4.set_yscale('log')
    ax4.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig('demo_nmr_magnet.png', dpi=150)
    print(f"\nSaved: demo_nmr_magnet.png")
    plt.close()

    # === NMR Application Guidelines ===
    print("\n" + "=" * 70)
    print("NMR Magnet Homogeneity Guidelines")
    print("=" * 70)
    print("""
Application         B0 [T]    DSV [cm]   Homogeneity [ppm]
------------------------------------------------------------
Clinical MRI        1.5-3.0    40-50      < 5 ppm (peak-peak)
Research MRI        7.0-9.4    20-30      < 10 ppm
High-field NMR      11.7-23.5  5-10       < 0.01 ppm (10 ppb)
Portable NMR        0.5-1.0    5-15       < 100 ppm

Key Harmonics for Shimming:
  Z0  (0,0)  : Main field offset
  Z1  (1,0)  : Z gradient
  X   (1,1)  : X gradient
  Y   (1,-1) : Y gradient
  Z2  (2,0)  : Quadratic Z
  ZX  (2,1)  : ZX cross-term
  ZY  (2,-1) : ZY cross-term
  X2-Y2 (2,2): Quadratic XY
  XY  (2,-2) : XY cross-term
  Higher orders (l >= 3): Fine shimming

Shimming Methods:
  1. Passive: Iron shims positioned around bore
  2. Active: Shim coil currents adjusted
  3. Gradient coils: Linear corrections
  4. High-order: Additional shim coil sets
""")

    print("=" * 70)


if __name__ == '__main__':
    main()
