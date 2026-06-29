"""
Coil wound on magnetic core - Frequency-dependent characteristics analysis.

This example demonstrates the physics of a coil wound on a magnetic core,
showing how conductivity (sigma) and permeability (mu_r) affect:
- Inductance vs frequency
- Resistance vs frequency (eddy current losses)
- Skin depth in the core
- Q-factor vs frequency

Physical models used:
- DC: Pure magnetic circuit (Radia MSC)
- AC: Nonlocal SIBC for conductive magnetic materials

Reference:
- Bilicz et al., "Nonlocal SIBC for high-conductivity regions", ISEM 2023
"""

import numpy as np
import matplotlib.pyplot as plt

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # H/m
EPS_0 = 8.854187817e-12  # F/m


def skin_depth(freq, sigma, mu_r):
    """Calculate skin depth in meters."""
    if freq <= 0 or sigma <= 0:
        return np.inf
    omega = 2 * np.pi * freq
    mu = MU_0 * mu_r
    return np.sqrt(2 / (omega * mu * sigma))


def local_surface_impedance(freq, sigma, mu_r):
    """Local surface impedance Zs = (1+j) / (sigma * delta)."""
    delta = skin_depth(freq, sigma, mu_r)
    if np.isinf(delta):
        return 0 + 0j
    return (1 + 1j) / (sigma * delta)


class MagneticCoreMaterial:
    """Magnetic core material properties."""

    def __init__(self, name, mu_r, sigma, description=""):
        self.name = name
        self.mu_r = mu_r  # Relative permeability
        self.sigma = sigma  # Conductivity [S/m]
        self.description = description

    def skin_depth(self, freq):
        return skin_depth(freq, self.sigma, self.mu_r)

    def surface_impedance(self, freq):
        return local_surface_impedance(freq, self.sigma, self.mu_r)


class CoilOnCore:
    """Model of a coil wound on a magnetic core."""

    def __init__(self, core_material, N_turns, core_length, core_area,
                 wire_diameter, wire_sigma=5.8e7):
        """
        Parameters:
        -----------
        core_material : MagneticCoreMaterial
            Core material properties
        N_turns : int
            Number of turns
        core_length : float
            Magnetic path length [m]
        core_area : float
            Core cross-section area [m^2]
        wire_diameter : float
            Wire diameter [m]
        wire_sigma : float
            Wire conductivity [S/m] (default: copper)
        """
        self.core = core_material
        self.N = N_turns
        self.l_core = core_length
        self.A_core = core_area
        self.d_wire = wire_diameter
        self.sigma_wire = wire_sigma

        # Wire cross-section area
        self.A_wire = np.pi * (wire_diameter / 2) ** 2

        # Approximate wire length (single layer, close wound)
        self.l_wire = N_turns * np.pi * np.sqrt(core_area / np.pi) * 2

    def dc_inductance(self):
        """DC inductance assuming infinite permeability core."""
        # L = mu_0 * mu_r * N^2 * A / l
        return MU_0 * self.core.mu_r * self.N**2 * self.A_core / self.l_core

    def dc_resistance(self):
        """DC resistance of wire."""
        return self.l_wire / (self.sigma_wire * self.A_wire)

    def effective_permeability(self, freq):
        """
        Effective permeability considering eddy currents in core.

        At high frequencies, eddy currents shield the core interior,
        reducing effective permeability.
        """
        delta = self.core.skin_depth(freq)

        # Characteristic dimension of core (assume circular for simplicity)
        a = np.sqrt(self.A_core / np.pi)

        if delta > 10 * a:
            # Low frequency: full penetration
            return self.core.mu_r
        elif delta < 0.1 * a:
            # High frequency: surface layer only
            # Effective area ~ 2*pi*a*delta (annular ring)
            area_ratio = 2 * delta / a
            return self.core.mu_r * area_ratio
        else:
            # Transition region: approximate interpolation
            # Using Bessel function approximation for cylindrical core
            x = a / delta
            # For large x: mu_eff/mu_r ~ 2/x * (ber + j*bei) / (ber' + j*bei')
            # Simplified: mu_eff/mu_r ~ 2*delta/a for x >> 1
            return self.core.mu_r * (2 / x) * (1 - 1j / x)

    def inductance(self, freq):
        """Frequency-dependent inductance."""
        mu_eff = self.effective_permeability(freq)
        L_complex = MU_0 * mu_eff * self.N**2 * self.A_core / self.l_core
        return L_complex

    def wire_ac_resistance(self, freq):
        """AC resistance of wire due to skin effect."""
        delta_wire = skin_depth(freq, self.sigma_wire, 1.0)  # mu_r=1 for copper
        a_wire = self.d_wire / 2

        if delta_wire > a_wire:
            # DC regime
            return self.dc_resistance()
        else:
            # Skin effect regime
            # Approximate: R_ac/R_dc ~ a/(2*delta) for delta << a
            return self.dc_resistance() * a_wire / (2 * delta_wire)

    def core_loss_resistance(self, freq):
        """
        Equivalent resistance representing core losses.

        Core losses = eddy current + hysteresis
        P_eddy ~ f^2 * B^2
        P_hyst ~ f * B^n (n ~ 1.6-2.0)

        For simplicity, model as R_core ~ omega * L * tan(delta_core)
        """
        if freq <= 0:
            return 0

        omega = 2 * np.pi * freq
        delta = self.core.skin_depth(freq)
        a = np.sqrt(self.A_core / np.pi)

        # Eddy current loss factor (proportional to f^2)
        # For thin laminations or ferrite, this is reduced
        if self.core.sigma > 1e4:  # Conductive core (silicon steel, etc.)
            # Bulk eddy current loss
            if delta < a:
                # High frequency: loss in skin depth layer
                loss_factor = (a / delta) ** 2 * 0.01  # Empirical
            else:
                # Low frequency: uniform eddy currents
                loss_factor = (omega * MU_0 * self.core.mu_r * self.core.sigma * a**2) ** 2 * 1e-6
        else:
            # Low conductivity (ferrite)
            loss_factor = 0.001  # Small residual loss

        L_dc = self.dc_inductance()
        return omega * np.real(L_dc) * loss_factor

    def total_impedance(self, freq):
        """Total impedance Z = R + j*omega*L."""
        omega = 2 * np.pi * freq
        L = self.inductance(freq)
        R_wire = self.wire_ac_resistance(freq)
        R_core = self.core_loss_resistance(freq)

        Z = (R_wire + R_core) + 1j * omega * L
        return Z

    def quality_factor(self, freq):
        """Q = omega*L / R."""
        Z = self.total_impedance(freq)
        omega = 2 * np.pi * freq
        L = np.abs(self.inductance(freq))
        R = np.real(Z)

        if R < 1e-15:
            return np.inf
        return omega * L / R


def analyze_coil_characteristics():
    """Analyze and plot coil characteristics for different core materials."""

    # Define core materials
    materials = [
        MagneticCoreMaterial("Air", 1.0, 0, "No core (air)"),
        MagneticCoreMaterial("Ferrite (MnZn)", 2000, 0.1, "Low loss ferrite"),
        MagneticCoreMaterial("Ferrite (NiZn)", 200, 1e-4, "High frequency ferrite"),
        MagneticCoreMaterial("Silicon Steel", 4000, 2e6, "Laminated transformer core"),
        MagneticCoreMaterial("Pure Iron", 5000, 1e7, "Solid iron core"),
    ]

    # Coil parameters
    N_turns = 100
    core_length = 0.1  # 10 cm
    core_area = 1e-4   # 1 cm^2
    wire_diameter = 0.5e-3  # 0.5 mm

    # Frequency range
    freqs = np.logspace(1, 8, 100)  # 10 Hz to 100 MHz

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    for mat in materials:
        coil = CoilOnCore(mat, N_turns, core_length, core_area, wire_diameter)

        # Calculate characteristics
        L_values = np.array([np.abs(coil.inductance(f)) for f in freqs])
        R_values = np.array([np.real(coil.total_impedance(f)) for f in freqs])
        Q_values = np.array([coil.quality_factor(f) for f in freqs])
        delta_values = np.array([mat.skin_depth(f) for f in freqs])

        # Plot inductance vs frequency
        axes[0, 0].loglog(freqs, L_values * 1e3, label=mat.name)

        # Plot resistance vs frequency
        axes[0, 1].loglog(freqs, R_values, label=mat.name)

        # Plot Q factor vs frequency
        axes[1, 0].semilogx(freqs, Q_values, label=mat.name)

        # Plot skin depth vs frequency (only for conductive materials)
        if mat.sigma > 0:
            valid_mask = np.isfinite(delta_values) & (delta_values < 1)
            if np.any(valid_mask):
                axes[1, 1].loglog(freqs[valid_mask], delta_values[valid_mask] * 1e3,
                                  label=mat.name)

    # Formatting
    axes[0, 0].set_xlabel('Frequency [Hz]')
    axes[0, 0].set_ylabel('Inductance [mH]')
    axes[0, 0].set_title('Inductance vs Frequency')
    axes[0, 0].legend(loc='best', fontsize=8)
    axes[0, 0].grid(True, which='both', alpha=0.3)
    axes[0, 0].set_ylim([1e-3, 1e3])

    axes[0, 1].set_xlabel('Frequency [Hz]')
    axes[0, 1].set_ylabel('Resistance [Ohm]')
    axes[0, 1].set_title('Total Resistance vs Frequency')
    axes[0, 1].legend(loc='best', fontsize=8)
    axes[0, 1].grid(True, which='both', alpha=0.3)

    axes[1, 0].set_xlabel('Frequency [Hz]')
    axes[1, 0].set_ylabel('Q Factor')
    axes[1, 0].set_title('Quality Factor vs Frequency')
    axes[1, 0].legend(loc='best', fontsize=8)
    axes[1, 0].grid(True, which='both', alpha=0.3)
    axes[1, 0].set_ylim([0, 200])

    axes[1, 1].set_xlabel('Frequency [Hz]')
    axes[1, 1].set_ylabel('Skin Depth [mm]')
    axes[1, 1].set_title('Skin Depth in Core vs Frequency')
    axes[1, 1].legend(loc='best', fontsize=8)
    axes[1, 1].grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig('coil_on_magnetic_core_analysis.png', dpi=150)
    print("Saved: coil_on_magnetic_core_analysis.png")
    plt.close()

    # Print summary table
    print("\n" + "="*80)
    print("Coil on Magnetic Core - Characteristics Summary")
    print("="*80)
    print(f"Coil: {N_turns} turns, core length={core_length*100:.1f}cm, "
          f"core area={core_area*1e4:.2f}cm^2")
    print("-"*80)
    print(f"{'Material':<20} {'mu_r':<10} {'sigma [S/m]':<12} "
          f"{'L_DC [mH]':<12} {'R_DC [Ohm]':<12} {'Q@1kHz':<10}")
    print("-"*80)

    for mat in materials:
        coil = CoilOnCore(mat, N_turns, core_length, core_area, wire_diameter)
        L_dc = coil.dc_inductance()
        R_dc = coil.dc_resistance()
        Q_1k = coil.quality_factor(1000)

        print(f"{mat.name:<20} {mat.mu_r:<10.0f} {mat.sigma:<12.2e} "
              f"{L_dc*1e3:<12.4f} {R_dc:<12.4f} {Q_1k:<10.1f}")

    print("="*80)


def compare_sibc_requirements():
    """
    Compare when local vs nonlocal SIBC is needed.

    Key criterion: skin depth vs characteristic dimension
    - delta >> d: DC/quasi-static, no SIBC needed
    - delta ~ d: Nonlocal SIBC required (captures internal field distribution)
    - delta << d: Local SIBC sufficient (thin skin approximation)
    """

    print("\n" + "="*80)
    print("SIBC Requirement Analysis")
    print("="*80)

    # Core characteristic dimension (assume 1 cm diameter)
    d_core = 0.01  # 10 mm

    materials = [
        ("Ferrite (MnZn)", 2000, 0.1),
        ("Silicon Steel", 4000, 2e6),
        ("Pure Iron", 5000, 1e7),
        ("Copper", 1, 5.8e7),
    ]

    frequencies = [50, 1000, 10000, 100000, 1e6]

    print(f"\nCore diameter: {d_core*1000:.1f} mm")
    print("-"*80)
    print(f"{'Material':<20}", end="")
    for f in frequencies:
        if f >= 1e6:
            print(f"{'f='+str(int(f/1e6))+'MHz':<12}", end="")
        elif f >= 1000:
            print(f"{'f='+str(int(f/1000))+'kHz':<12}", end="")
        else:
            print(f"{'f='+str(int(f))+'Hz':<12}", end="")
    print()
    print("-"*80)

    for name, mu_r, sigma in materials:
        print(f"{name:<20}", end="")
        for f in frequencies:
            delta = skin_depth(f, sigma, mu_r)
            ratio = delta / d_core

            if ratio > 10:
                status = "DC"  # No SIBC needed
            elif ratio > 1:
                status = "Nonloc"  # Nonlocal SIBC
            elif ratio > 0.1:
                status = "Local"  # Local SIBC OK
            else:
                status = "Surf"  # Surface current only

            if np.isinf(delta):
                print(f"{'inf':<12}", end="")
            else:
                print(f"{delta*1000:.2f}mm({status})", end="")
                # Adjust spacing
                extra = 12 - len(f"{delta*1000:.2f}mm({status})")
                print(" " * max(0, extra), end="")
        print()

    print("-"*80)
    print("Legend: DC=quasi-static, Nonloc=Nonlocal SIBC, Local=Local SIBC, Surf=Surface only")
    print("="*80)


def main():
    """Main function."""
    print("Coil on Magnetic Core - Frequency Analysis")
    print("=" * 50)

    # Run analyses
    analyze_coil_characteristics()
    compare_sibc_requirements()

    print("\nAnalysis complete!")


if __name__ == '__main__':
    main()
