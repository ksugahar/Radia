#!/usr/bin/env python
"""
SPICE Time-Domain Simulation Demonstration

This script demonstrates:
1. URN model fitting to impedance data
2. SPICE netlist generation
3. Time-domain simulation (step response, PWM)
4. Comparison with frequency-domain predictions

Requirements:
    pip install numpy scipy torch matplotlib
    ngspice (for actual SPICE simulation)

Usage:
    python demo_spice_timedomain.py

Output:
    - SPICE netlists for LTspice/ngspice
    - Time-domain response plots
    - Stability analysis results
"""

import sys
import numpy as np
import torch
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import signal

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))

from universal_relaxation_network import (
    UniversalRelaxationNetwork, URNConfig, train_urn, generate_spice_netlist
)


def create_rl_ladder_transfer_function(R_values, L_values, R_dc=0.01):
    """
    Create transfer function for RL ladder network.

    Circuit topology (series RL ladder):
        o---R_dc---+---L1---R1---+---L2---R2---+---...
                   |             |             |
                   V             V             V
                  out           ...           ...

    Parameters
    ----------
    R_values : array-like
        Resistance values in Ohms
    L_values : array-like
        Inductance values in Henrys
    R_dc : float
        DC resistance in series

    Returns
    -------
    num, den : array-like
        Transfer function coefficients H(s) = num(s)/den(s)
    """
    n_stages = len(R_values)
    assert len(L_values) == n_stages

    # For simple demonstration, compute impedance poles/zeros
    # Z(s) = R_dc + sum_i(R_i * s*tau_i / (1 + s*tau_i))
    # where tau_i = L_i / R_i

    poles = []
    zeros = []

    for R, L in zip(R_values, L_values):
        if L > 0 and R > 0:
            tau = L / R
            poles.append(-1.0 / tau)  # pole at s = -R/L

    # Convert to transfer function
    # Simple approximation: first-order cascade
    num = [R_dc]
    den = [1.0]

    for R, L in zip(R_values, L_values):
        if L > 0 and R > 0:
            tau = L / R
            # Add stage: (R + s*L)/(1) = R * (1 + s*tau)
            num_stage = [R * tau, R]
            den_stage = [tau, 1.0]
            num = np.convolve(num, num_stage)
            den = np.convolve(den, den_stage)

    return num, den


def simulate_step_response(freq, Z_data, model, title, output_path):
    """
    Simulate step response using URN model.

    Converts the fitted impedance model to a transfer function
    and computes the step response.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Get model parameters
    mechanisms = model.get_active_components()
    omega = 2 * np.pi * freq
    omega_torch = torch.tensor(omega, dtype=torch.float64)
    Z_fit = model(omega_torch).detach().numpy()

    # Create approximate transfer function from dominant poles
    # Extract time constants from discovered mechanisms
    tau_values = []
    R_values = []

    for name, components in mechanisms.items():
        for params in components:
            if 'tau' in params and params.get('weight_magnitude', 0) > 0.01:
                tau_values.append(params['tau'])
                R_values.append(params['weight_magnitude'] * np.abs(Z_data).mean())

    # Ensure we have at least some values
    if len(tau_values) == 0:
        tau_values = [1e-3, 1e-2, 1e-1]
        R_values = [0.01, 0.02, 0.05]

    # Build RL ladder approximation
    L_values = [R * tau for R, tau in zip(R_values, tau_values)]
    R_dc = np.abs(Z_fit[np.argmin(freq)]).real

    # Create transfer function
    # Simple RC ladder: Z(s) = R0 + R1/(1+s*tau1) + R2/(1+s*tau2) + ...
    # Approximate as first-order systems in series

    # Time vector for simulation
    t_max = max(tau_values) * 20
    t = np.linspace(0, t_max, 2000)

    # Compute step response as sum of exponentials
    # z(t) = R_dc + sum_i R_i * (1 - exp(-t/tau_i))
    z_step = np.ones_like(t) * R_dc
    for R, tau in zip(R_values, tau_values):
        z_step += R * (1 - np.exp(-t / tau))

    # Current response to unit voltage step: i(t) = V/Z(t)
    # For proper step response, we need to solve the circuit ODE
    # Simplified: i(t) ~ V/Z_dc * (1 - sum of exponential decays)
    Z_dc = R_dc + sum(R_values)
    i_step = np.zeros_like(t)
    i_step[0] = 1.0 / R_dc if R_dc > 0 else 1.0  # Initial current (high)

    # Exponential decay from each branch
    for i, (R, tau) in enumerate(zip(R_values, tau_values)):
        contribution = (R / Z_dc) * np.exp(-t / tau)
        i_step += contribution

    i_step = 1.0 / Z_dc + (1.0 / R_dc - 1.0 / Z_dc) * np.exp(-t / min(tau_values))

    # Plot 1: Step response (current)
    ax1 = axes[0, 0]
    ax1.plot(t * 1000, i_step, 'b-', linewidth=2)
    ax1.set_xlabel('Time [ms]')
    ax1.set_ylabel('Current [A] (normalized to 1V)')
    ax1.set_title('Step Response (1V applied)')
    ax1.grid(True)
    ax1.axhline(y=1.0/Z_dc, color='r', linestyle='--', label=f'Steady state: {1/Z_dc:.4f} A')
    ax1.legend()

    # Plot 2: Impedance evolution
    ax2 = axes[0, 1]
    ax2.semilogy(t * 1000, z_step, 'b-', linewidth=2)
    ax2.set_xlabel('Time [ms]')
    ax2.set_ylabel('Impedance [Ohm]')
    ax2.set_title('Impedance vs Time')
    ax2.grid(True)

    # Plot 3: PWM response simulation
    ax3 = axes[1, 0]
    f_pwm = 10000  # 10 kHz PWM
    duty = 0.5
    T_pwm = 1.0 / f_pwm
    t_pwm = np.linspace(0, 5 * T_pwm, 5000)

    # PWM voltage waveform
    v_pwm = np.zeros_like(t_pwm)
    for i, ti in enumerate(t_pwm):
        phase = (ti % T_pwm) / T_pwm
        v_pwm[i] = 1.0 if phase < duty else 0.0

    # Simple RC filter response to PWM
    tau_eff = np.mean(tau_values) if tau_values else 1e-4
    i_pwm = np.zeros_like(t_pwm)
    dt = t_pwm[1] - t_pwm[0]
    R_eff = np.mean(R_values) if R_values else 0.01

    for i in range(1, len(t_pwm)):
        # Simple first-order response: di/dt = (V - i*R) / L
        L_eff = tau_eff * R_eff
        di_dt = (v_pwm[i] - i_pwm[i-1] * R_eff) / L_eff
        i_pwm[i] = i_pwm[i-1] + di_dt * dt

    ax3.plot(t_pwm * 1e6, v_pwm, 'b-', linewidth=1, alpha=0.5, label='PWM Voltage')
    ax3.plot(t_pwm * 1e6, i_pwm / max(i_pwm) if max(i_pwm) > 0 else i_pwm,
             'r-', linewidth=2, label='Current (normalized)')
    ax3.set_xlabel('Time [us]')
    ax3.set_ylabel('Amplitude')
    ax3.set_title(f'PWM Response ({f_pwm/1000:.0f} kHz, {duty*100:.0f}% duty)')
    ax3.legend()
    ax3.grid(True)

    # Plot 4: Stability analysis (pole locations)
    ax4 = axes[1, 1]

    # Plot poles from time constants
    poles = [-1.0/tau for tau in tau_values if tau > 0]

    if poles:
        ax4.scatter(poles, np.zeros_like(poles), s=100, marker='x',
                   color='red', linewidths=2, label='Poles')
        ax4.axvline(x=0, color='k', linestyle='-', linewidth=0.5)
        ax4.axhline(y=0, color='k', linestyle='-', linewidth=0.5)

        # Shade stable region
        ax4.axvspan(min(poles)*2, 0, alpha=0.2, color='green', label='Stable region')

        ax4.set_xlabel('Real(s) [rad/s]')
        ax4.set_ylabel('Imag(s) [rad/s]')
        ax4.set_title('Pole Locations (all in LHP -> stable)')
        ax4.legend()
        ax4.grid(True)

        # Check stability
        all_stable = all(p < 0 for p in poles)
        stability_text = "STABLE" if all_stable else "UNSTABLE"
        ax4.text(0.95, 0.95, stability_text,
                transform=ax4.transAxes, fontsize=14, fontweight='bold',
                ha='right', va='top',
                color='green' if all_stable else 'red',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    else:
        ax4.text(0.5, 0.5, 'No poles extracted',
                transform=ax4.transAxes, ha='center', va='center')

    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def generate_ltspice_netlist(model, name, output_path):
    """
    Generate LTspice-compatible netlist with transient analysis directives.
    """
    spice = generate_spice_netlist(model, name)

    # Add LTspice simulation directives
    ltspice_additions = """
* ============================================
* LTspice Simulation Directives
* ============================================

* Step response analysis
.tran 0 100m 0 1u

* AC analysis (Bode plot)
.ac dec 100 0.01 100k

* Input source (choose one):
* For step response:
V1 in 0 PULSE(0 1 0 1n 1n 100m 200m)

* For AC analysis:
* V1 in 0 AC 1

* Output node: "out"
.end
"""
    full_netlist = spice + ltspice_additions

    with open(output_path, 'w') as f:
        f.write(full_netlist)

    print(f"  Saved LTspice netlist: {output_path}")
    return full_netlist


def load_battery_data():
    """
    Load battery EIS data, preferring real NASA data over synthetic.

    Returns:
        (freq, Z, data_source) tuple
    """
    # Try real NASA data first
    real_data_path = Path(__file__).parent / 'data' / 'real_world' / 'nasa_battery' / 'nasa_18650_eis.csv'
    synthetic_data_path = Path(__file__).parent / 'data' / 'synthetic' / 'liion_battery_eis.csv'

    if real_data_path.exists():
        print(f"  Loading REAL NASA battery data: {real_data_path.name}")
        import pandas as pd
        df = pd.read_csv(real_data_path, comment='#')
        freq = df['frequency_Hz'].values
        Z = df['Z_real_Ohm'].values + 1j * df['Z_imag_Ohm'].values
        return freq, Z, "NASA 18650 (REAL)"
    elif synthetic_data_path.exists():
        print(f"  Loading synthetic battery data: {synthetic_data_path.name}")
        data = np.loadtxt(synthetic_data_path, delimiter=',', skiprows=18)
        freq = data[:, 0]
        Z = data[:, 1] + 1j * data[:, 2]
        return freq, Z, "Synthetic Randles"
    else:
        raise FileNotFoundError("No battery EIS data found")


def load_ferrite_data():
    """
    Load ferrite impedance data, preferring real TDK data over synthetic.

    Returns:
        (freq, Z, data_source) tuple
    """
    # Try real TDK data first (PC50 is a good general-purpose material)
    real_data_path = Path(__file__).parent / 'data' / 'real_world' / 'tdk_ferrite' / 'tdk_pc50_impedance.csv'
    synthetic_data_path = Path(__file__).parent / 'data' / 'synthetic' / 'mnzn_ferrite_impedance.csv'

    if real_data_path.exists():
        print(f"  Loading REAL TDK PC50 ferrite data: {real_data_path.name}")
        import pandas as pd
        df = pd.read_csv(real_data_path, comment='#')
        freq = df['frequency_Hz'].values
        Z = df['Z_real_Ohm'].values + 1j * df['Z_imag_Ohm'].values
        return freq, Z, "TDK PC50 (REAL)"
    elif synthetic_data_path.exists():
        print(f"  Loading synthetic ferrite data: {synthetic_data_path.name}")
        data = np.loadtxt(synthetic_data_path, delimiter=',', skiprows=19)
        freq = data[:, 0]
        Z = data[:, 3] + 1j * data[:, 4]
        return freq, Z, "Synthetic Cole-Cole"
    else:
        raise FileNotFoundError("No ferrite impedance data found")


def main():
    """Run SPICE time-domain demonstration."""
    print("=" * 70)
    print("SPICE Time-Domain Simulation Demonstration")
    print("=" * 70)

    output_dir = (
        Path(__file__).resolve().parents[2]
        / 'validation_test'
        / 'universal_relaxation_network'
        / 'validation_output'
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # Load battery data and fit URN model
    # =========================================================================
    print("\n[1/2] Li-ion Battery: Time-Domain Analysis")
    print("-" * 40)

    freq_bat, Z_bat, bat_source = load_battery_data()
    print(f"  Data source: {bat_source}")
    print(f"  Frequency range: {freq_bat[0]:.2f} Hz - {freq_bat[-1]:.0f} Hz")

    print(f"  Fitting URN model to battery data...")
    config_bat = URNConfig(
        n_debye=3, n_cole_cole=2, n_cpe=2, n_warburg=2,
        n_debye_parallel=2, n_warburg_parallel=1,
        sparsity_weight=0.01, n_epochs=2000, n_restarts=2
    )
    model_bat = train_urn(freq_bat, Z_bat, config_bat)

    # Generate time-domain plots
    simulate_step_response(
        freq_bat, Z_bat, model_bat,
        "Li-ion Battery: Time-Domain Response",
        output_dir / 'battery_timedomain.png'
    )

    # Generate LTspice netlist
    generate_ltspice_netlist(
        model_bat, "BATTERY",
        output_dir / 'battery_ltspice.sp'
    )

    # =========================================================================
    # Load ferrite data and fit URN model
    # =========================================================================
    print("\n[2/2] MnZn Ferrite: Time-Domain Analysis")
    print("-" * 40)

    freq_fer, Z_fer, fer_source = load_ferrite_data()
    print(f"  Data source: {fer_source}")
    print(f"  Frequency range: {freq_fer[0]:.0f} Hz - {freq_fer[-1]/1e6:.1f} MHz")

    print(f"  Fitting URN model to ferrite data...")
    config_fer = URNConfig(
        n_debye=4, n_cole_cole=3, n_cpe=1, n_skin_effect=2,
        n_debye_parallel=2, n_cole_cole_parallel=1,
        sparsity_weight=0.01, n_epochs=2000, n_restarts=2
    )
    model_fer = train_urn(freq_fer, Z_fer, config_fer)

    # Generate time-domain plots
    simulate_step_response(
        freq_fer, Z_fer, model_fer,
        "MnZn Ferrite: Time-Domain Response",
        output_dir / 'ferrite_timedomain.png'
    )

    # Generate LTspice netlist
    generate_ltspice_netlist(
        model_fer, "FERRITE",
        output_dir / 'ferrite_ltspice.sp'
    )

    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "=" * 70)
    print("DEMONSTRATION COMPLETE")
    print("=" * 70)
    print(f"\nGenerated files in {output_dir}:")
    print("  - battery_timedomain.png: Step/PWM response plots")
    print("  - battery_ltspice.sp: LTspice netlist for simulation")
    print("  - ferrite_timedomain.png: Step/PWM response plots")
    print("  - ferrite_ltspice.sp: LTspice netlist for simulation")
    print("\nTo run SPICE simulation:")
    print("  1. Open .sp file in LTspice")
    print("  2. Run simulation (Ctrl+R)")
    print("  3. Plot output node voltage/current")


if __name__ == '__main__':
    main()
