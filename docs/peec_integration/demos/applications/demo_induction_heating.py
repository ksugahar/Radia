"""
Induction Heating Analysis with PEEC and Surface Impedance

This demo shows induction heating coil analysis:
1. Coil-workpiece electromagnetic coupling
2. Surface impedance with skin effect (ESIM model)
3. Power loss distribution in workpiece
4. Frequency selection for optimal heating depth

Physical setup:
- Heating coil: Circular, water-cooled copper tube
- Workpiece: Cylindrical steel or aluminum

Author: Radia Development Team
Date: 2026-01-16
"""

import numpy as np
import matplotlib.pyplot as plt


# Physical constants
MU_0 = 4 * np.pi * 1e-7  # [H/m]


def calc_skin_depth(freq: float, sigma: float, mu_r: float = 1.0) -> float:
    """Calculate skin depth delta = sqrt(2 / (omega * mu * sigma))."""
    omega = 2 * np.pi * freq
    if omega < 1e-10:
        return np.inf
    return np.sqrt(2 / (omega * MU_0 * mu_r * sigma))


def calc_surface_impedance(freq: float, sigma: float,
                            mu_r: float = 1.0) -> complex:
    """
    Calculate surface impedance Zs = (1+j) / (sigma * delta).

    For high-frequency (strong skin effect):
        Zs = sqrt(j * omega * mu / sigma)
    """
    omega = 2 * np.pi * freq
    if omega < 1e-10:
        return 0.0
    delta = calc_skin_depth(freq, sigma, mu_r)
    Zs = (1 + 1j) / (sigma * delta)
    return Zs


def calc_power_density(H_surface: float, Zs: complex) -> dict:
    """
    Calculate power density from surface H-field and impedance.

    P = Re(Zs) * |H|^2 [W/m^2]

    Parameters:
        H_surface: Surface magnetic field [A/m]
        Zs: Surface impedance [Ohm]

    Returns:
        dict with active power, reactive power, power factor
    """
    # Active power (heating)
    P_active = np.real(Zs) * H_surface**2

    # Reactive power (magnetic energy)
    Q_reactive = np.imag(Zs) * H_surface**2

    # Apparent power
    S_apparent = np.abs(Zs) * H_surface**2

    # Power factor
    pf = P_active / S_apparent if S_apparent > 0 else 0

    return {
        'P_active': P_active,       # W/m^2
        'Q_reactive': Q_reactive,   # VAr/m^2
        'S_apparent': S_apparent,   # VA/m^2
        'power_factor': pf
    }


def calc_temperature_rise(P_density: float, time: float,
                           rho: float, c_p: float, delta: float) -> float:
    """
    Estimate temperature rise in skin depth layer.

    Delta_T = P * t / (rho * c_p * delta)

    Parameters:
        P_density: Power density [W/m^2]
        time: Heating time [s]
        rho: Material density [kg/m^3]
        c_p: Specific heat [J/(kg*K)]
        delta: Skin depth [m]

    Returns:
        Temperature rise [K]
    """
    # Energy deposited per unit area
    E = P_density * time

    # Volume of heated layer per unit area
    V = delta  # m^3/m^2 = m

    # Mass per unit area
    m = rho * V

    # Temperature rise
    if m > 0:
        delta_T = E / (m * c_p)
    else:
        delta_T = 0

    return delta_T


def calc_coil_workpiece_coupling(r_coil: float, r_wp: float,
                                   gap: float, n_turns: int) -> float:
    """
    Calculate coil-workpiece mutual inductance approximation.

    Uses simplified model for coaxial coil and cylinder.

    Parameters:
        r_coil: Coil radius [m]
        r_wp: Workpiece radius [m]
        gap: Coil-workpiece gap [m]
        n_turns: Number of coil turns

    Returns:
        Mutual inductance [H]
    """
    # Effective distance
    d = gap + (r_coil - r_wp) / 2 if r_coil > r_wp else gap

    # Simplified mutual inductance
    # M ~ mu_0 * n * pi * r_wp^2 / (2 * d) for d << r_coil
    M = MU_0 * n_turns * np.pi * min(r_coil, r_wp)**2 / (2 * d + r_wp)

    return M


# Material database
MATERIALS = {
    'carbon_steel': {
        'sigma': 5e6,        # S/m (varies with temperature)
        'mu_r': 100,         # At low field (decreases with temperature)
        'rho': 7850,         # kg/m^3
        'c_p': 480,          # J/(kg*K)
        'T_curie': 770,      # C
    },
    'stainless_steel': {
        'sigma': 1.4e6,      # S/m (austenitic 304)
        'mu_r': 1.0,         # Non-magnetic
        'rho': 8000,
        'c_p': 500,
        'T_curie': None,
    },
    'aluminum': {
        'sigma': 3.5e7,      # S/m
        'mu_r': 1.0,
        'rho': 2700,
        'c_p': 900,
        'T_curie': None,
    },
    'copper': {
        'sigma': 5.8e7,      # S/m
        'mu_r': 1.0,
        'rho': 8960,
        'c_p': 385,
        'T_curie': None,
    },
}


def main():
    """Induction heating analysis demo."""
    print("=" * 70)
    print("Induction Heating Analysis with PEEC")
    print("=" * 70)

    # === System Parameters ===
    # Heating coil
    r_coil = 0.05        # 50 mm radius
    n_turns = 5          # 5 turns
    tube_od = 10e-3      # 10 mm outer diameter
    tube_wall = 1e-3     # 1 mm wall thickness
    coil_length = 0.05   # 50 mm coil length

    # Workpiece
    r_wp = 0.03          # 30 mm radius
    gap = 5e-3           # 5 mm coil-workpiece gap

    # Operating current
    I_coil = 100         # 100 A (RMS)

    print("\n--- System Parameters ---")
    print(f"Coil: r = {r_coil*1e3:.1f} mm, N = {n_turns}, tube OD = {tube_od*1e3:.1f} mm")
    print(f"Workpiece: r = {r_wp*1e3:.1f} mm")
    print(f"Gap: {gap*1e3:.1f} mm")
    print(f"Coil current: {I_coil:.0f} A (RMS)")

    # === Material Comparison ===
    print("\n" + "=" * 70)
    print("Material Comparison at Different Frequencies")
    print("=" * 70)

    freqs = [1e3, 10e3, 100e3, 500e3]  # 1 kHz to 500 kHz

    # Surface H-field estimation (simplified)
    # H ~ n * I / (2 * (r_coil - r_wp + gap))
    H_surface = n_turns * I_coil / (2 * (gap + 0.01))

    print(f"\nEstimated surface H-field: {H_surface:.0f} A/m")

    for mat_name in ['carbon_steel', 'stainless_steel', 'aluminum']:
        mat = MATERIALS[mat_name]
        print(f"\n--- {mat_name.replace('_', ' ').title()} ---")
        print(f"sigma = {mat['sigma']:.1e} S/m, mu_r = {mat['mu_r']}")

        print(f"{'Freq':>10} {'delta':>12} {'Rs':>12} {'P':>12} {'eta_skin':>10}")
        print("-" * 60)

        for f in freqs:
            delta = calc_skin_depth(f, mat['sigma'], mat['mu_r'])
            Zs = calc_surface_impedance(f, mat['sigma'], mat['mu_r'])
            power = calc_power_density(H_surface, Zs)

            # Skin effect efficiency (how much of workpiece is heated)
            eta_skin = min(1.0, delta / r_wp)

            freq_str = f"{f/1e3:.0f} kHz"
            print(f"{freq_str:>10} {delta*1e3:>10.3f} mm {np.real(Zs)*1e6:>10.2f} uOhm {power['P_active']/1e6:>10.3f} MW/m^2 {eta_skin*100:>9.1f}%")

    # === Frequency Selection for Carbon Steel ===
    print("\n" + "=" * 70)
    print("Optimal Frequency Selection (Carbon Steel)")
    print("=" * 70)

    mat = MATERIALS['carbon_steel']

    # Target heating depths
    depths = [0.5e-3, 1e-3, 2e-3, 5e-3]  # 0.5 to 5 mm

    print("\nTarget depth -> Required frequency (delta = depth)")
    print("-" * 40)
    for depth in depths:
        # delta = sqrt(2 / (omega * mu * sigma))
        # omega = 2 / (delta^2 * mu * sigma)
        omega = 2 / (depth**2 * MU_0 * mat['mu_r'] * mat['sigma'])
        f_req = omega / (2 * np.pi)
        print(f"  {depth*1e3:.1f} mm -> {f_req/1e3:.1f} kHz")

    # === Power vs Frequency Analysis ===
    print("\n--- Power Distribution Analysis ---")

    freqs_fine = np.logspace(2, 6, 200)  # 100 Hz to 1 MHz
    P_surface = []
    deltas = []

    for f in freqs_fine:
        delta = calc_skin_depth(f, mat['sigma'], mat['mu_r'])
        Zs = calc_surface_impedance(f, mat['sigma'], mat['mu_r'])
        power = calc_power_density(H_surface, Zs)
        P_surface.append(power['P_active'])
        deltas.append(delta)

    P_surface = np.array(P_surface)
    deltas = np.array(deltas)

    # === Temperature Rise Estimation ===
    print("\n--- Temperature Rise Estimate (1 second heating) ---")

    t_heat = 1.0  # 1 second
    mat = MATERIALS['carbon_steel']

    print(f"{'Freq':>10} {'P [kW/m^2]':>14} {'delta [mm]':>12} {'Delta_T [K]':>12}")
    print("-" * 52)

    for f in [1e3, 10e3, 50e3, 100e3]:
        delta = calc_skin_depth(f, mat['sigma'], mat['mu_r'])
        Zs = calc_surface_impedance(f, mat['sigma'], mat['mu_r'])
        power = calc_power_density(H_surface, Zs)

        delta_T = calc_temperature_rise(
            power['P_active'], t_heat,
            mat['rho'], mat['c_p'], delta
        )

        freq_str = f"{f/1e3:.0f} kHz"
        print(f"{freq_str:>10} {power['P_active']/1e3:>12.1f} {delta*1e3:>12.3f} {delta_T:>12.1f}")

    # === Plot Results ===
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: Skin depth vs frequency for different materials
    ax1 = axes[0, 0]
    for mat_name, mat in MATERIALS.items():
        deltas_mat = [calc_skin_depth(f, mat['sigma'], mat['mu_r'])
                      for f in freqs_fine]
        ax1.loglog(freqs_fine/1e3, np.array(deltas_mat)*1e3,
                   linewidth=2, label=mat_name.replace('_', ' ').title())
    ax1.set_xlabel('Frequency [kHz]')
    ax1.set_ylabel('Skin Depth [mm]')
    ax1.set_title('Skin Depth vs Frequency')
    ax1.legend()
    ax1.grid(True, which='both', alpha=0.3)
    ax1.set_xlim([0.1, 1000])

    # Plot 2: Surface power density vs frequency
    ax2 = axes[0, 1]
    for mat_name in ['carbon_steel', 'stainless_steel', 'aluminum']:
        mat = MATERIALS[mat_name]
        P_mat = []
        for f in freqs_fine:
            Zs = calc_surface_impedance(f, mat['sigma'], mat['mu_r'])
            power = calc_power_density(H_surface, Zs)
            P_mat.append(power['P_active'])
        ax2.loglog(freqs_fine/1e3, np.array(P_mat)/1e6,
                   linewidth=2, label=mat_name.replace('_', ' ').title())
    ax2.set_xlabel('Frequency [kHz]')
    ax2.set_ylabel('Surface Power Density [MW/m^2]')
    ax2.set_title(f'Power Density vs Frequency (H = {H_surface:.0f} A/m)')
    ax2.legend()
    ax2.grid(True, which='both', alpha=0.3)

    # Plot 3: Surface resistance vs frequency
    ax3 = axes[1, 0]
    for mat_name in ['carbon_steel', 'stainless_steel', 'aluminum']:
        mat = MATERIALS[mat_name]
        Rs_mat = []
        for f in freqs_fine:
            Zs = calc_surface_impedance(f, mat['sigma'], mat['mu_r'])
            Rs_mat.append(np.real(Zs))
        ax3.loglog(freqs_fine/1e3, np.array(Rs_mat)*1e6,
                   linewidth=2, label=mat_name.replace('_', ' ').title())
    ax3.set_xlabel('Frequency [kHz]')
    ax3.set_ylabel('Surface Resistance [uOhm]')
    ax3.set_title('Surface Resistance vs Frequency')
    ax3.legend()
    ax3.grid(True, which='both', alpha=0.3)

    # Plot 4: Power factor vs frequency
    ax4 = axes[1, 1]
    for mat_name in ['carbon_steel', 'stainless_steel', 'aluminum']:
        mat = MATERIALS[mat_name]
        pf_mat = []
        for f in freqs_fine:
            Zs = calc_surface_impedance(f, mat['sigma'], mat['mu_r'])
            # Power factor for surface impedance is always cos(45) = 0.707
            # for ideal skin effect, but can vary near DC
            pf_mat.append(np.real(Zs) / np.abs(Zs) if np.abs(Zs) > 0 else 0)
        ax4.semilogx(freqs_fine/1e3, pf_mat,
                     linewidth=2, label=mat_name.replace('_', ' ').title())
    ax4.set_xlabel('Frequency [kHz]')
    ax4.set_ylabel('Power Factor')
    ax4.set_title('Power Factor vs Frequency')
    ax4.legend()
    ax4.grid(True, which='both', alpha=0.3)
    ax4.set_ylim([0, 1])

    plt.tight_layout()
    plt.savefig('demo_induction_heating.png', dpi=150)
    print(f"\nSaved: demo_induction_heating.png")
    plt.close()

    # === Application Guidelines ===
    print("\n" + "=" * 70)
    print("Application Guidelines")
    print("=" * 70)
    print("""
Frequency Selection by Application:

  Application           Material         Frequency    Skin Depth
  -------------------------------------------------------------------
  Through heating       Steel            1-10 kHz     3-10 mm
  Surface hardening     Steel            50-500 kHz   0.3-1 mm
  Brazing/Soldering     Various          10-100 kHz   0.5-3 mm
  Melting (small)       Metals           1-10 kHz     2-8 mm
  Sealing               Al foil          100-500 kHz  0.1-0.3 mm

Key Design Considerations:

  1. Higher mu_r (magnetic steel) -> Lower frequency needed
  2. Higher sigma (aluminum, copper) -> Higher frequency for same depth
  3. Power ~ sqrt(f) for constant skin depth heating
  4. Efficiency improves at higher frequency (higher Q)
  5. Above Curie point, mu_r -> 1 (steel becomes non-magnetic)
""")

    print("=" * 70)


if __name__ == '__main__':
    main()
