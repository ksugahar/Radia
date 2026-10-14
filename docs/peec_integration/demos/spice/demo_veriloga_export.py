"""
Verilog-A Export Demo for PEEC Models

This demo shows how to export PEEC analysis results to Verilog-A
for SPICE circuit simulation.

Features demonstrated:
1. Dowell skin effect model export
2. Debye/Cole-Cole dielectric model export
3. Complete PEEC segment export
4. Multi-pole Debye model export

Output files:
- demo_dowell_skin.va
- demo_debye_cap.va
- demo_cole_cole_cap.va
- demo_peec_segment.va
- demo_multi_debye.va

Author: Radia Development Team
Date: 2026-01-16
"""

import numpy as np

# Import Verilog-A generator
from radia.veriloga_generator import (
    VerilogAGenerator,
    DowellParams, DebyeParams, ColeColeParams, MultiDebyeParams,
    generate_dowell_veriloga,
    generate_debye_veriloga,
    generate_cole_cole_veriloga
)


def demo_dowell_skin_effect():
    """Demo: Dowell skin effect model for copper conductor."""
    print("=" * 60)
    print("Demo 1: Dowell Skin Effect Model")
    print("=" * 60)

    # Copper conductor parameters
    thickness = 0.1e-3   # 0.1 mm
    width = 10e-3        # 10 mm
    length = 0.1         # 100 mm
    sigma = 5.8e7        # Copper conductivity [S/m]
    mu_r = 1.0

    # Calculate DC parameters
    mu0 = 4 * np.pi * 1e-7
    R_dc = length / (sigma * width * thickness)
    L_int_dc = mu0 * mu_r * length * thickness / (3 * width)

    print(f"Conductor: {thickness*1e3:.2f} mm x {width*1e3:.1f} mm x {length*1e3:.1f} mm")
    print(f"Material: Copper (sigma = {sigma:.2e} S/m)")
    print(f"R_dc = {R_dc*1e3:.4f} mOhm")
    print(f"L_int_dc = {L_int_dc*1e9:.4f} nH")

    # Generate Verilog-A
    code = generate_dowell_veriloga(
        R_dc=R_dc,
        L_int_dc=L_int_dc,
        thickness=thickness,
        sigma=sigma,
        mu_r=mu_r,
        filename="demo_dowell_skin.va"
    )

    print(f"\nGenerated: demo_dowell_skin.va")
    print("\nFirst 30 lines of generated code:")
    print("-" * 40)
    for line in code.split('\n')[:30]:
        print(line)
    print("...")

    # Calculate skin depth at various frequencies
    print("\nSkin depth reference:")
    print(f"{'Frequency':>12} {'Skin depth':>15} {'xi = d/2/delta':>15}")
    print("-" * 45)
    for f in [1e3, 10e3, 100e3, 1e6, 10e6]:
        delta = np.sqrt(2 / (2*np.pi*f * mu0*mu_r * sigma))
        xi = (thickness/2) / delta
        print(f"{f/1e3:>10.0f} kHz {delta*1e6:>12.2f} um {xi:>15.3f}")


def demo_debye_capacitor():
    """Demo: Debye relaxation capacitor for FR4 substrate."""
    print("\n" + "=" * 60)
    print("Demo 2: Debye Relaxation Capacitor")
    print("=" * 60)

    # FR4 substrate parameters (typical values)
    eps_s = 4.5          # Static permittivity
    eps_inf = 4.2        # High-frequency permittivity
    tau = 10e-12         # Relaxation time (10 ps)
    sigma = 1e-12        # Very low DC conductivity

    # Parallel plate capacitor geometry
    area = 1e-6          # 1 mm^2
    thickness = 0.2e-3   # 0.2 mm substrate
    eps0 = 8.854e-12
    C_0 = eps0 * area / thickness  # Geometry factor

    print(f"Material: FR4 substrate")
    print(f"eps_s = {eps_s}, eps_inf = {eps_inf}")
    print(f"Relaxation time tau = {tau*1e12:.1f} ps")
    print(f"Geometry: {area*1e6:.1f} mm^2, {thickness*1e3:.2f} mm thick")
    print(f"C_0 = {C_0*1e15:.3f} fF")

    # Generate Verilog-A
    code = generate_debye_veriloga(
        C_0=C_0,
        eps_s=eps_s,
        eps_inf=eps_inf,
        tau=tau,
        sigma=sigma,
        filename="demo_debye_cap.va"
    )

    print(f"\nGenerated: demo_debye_cap.va")

    # Show frequency dependence
    print("\nComplex permittivity vs frequency:")
    print(f"{'Frequency':>12} {'eps_re':>12} {'eps_im':>12} {'tan(delta)':>12}")
    print("-" * 52)
    for f in [1e6, 10e6, 100e6, 1e9, 10e9]:
        omega = 2 * np.pi * f
        omega_tau = omega * tau
        denom = 1 + omega_tau**2
        eps_re = eps_inf + (eps_s - eps_inf) / denom
        eps_im = (eps_s - eps_inf) * omega_tau / denom
        tan_d = eps_im / eps_re
        freq_str = f"{f/1e6:.0f} MHz" if f < 1e9 else f"{f/1e9:.0f} GHz"
        print(f"{freq_str:>12} {eps_re:>12.4f} {eps_im:>12.6f} {tan_d:>12.6f}")


def demo_cole_cole_capacitor():
    """Demo: Cole-Cole model for ferrite material."""
    print("\n" + "=" * 60)
    print("Demo 3: Cole-Cole Relaxation Capacitor")
    print("=" * 60)

    # Ferrite material parameters (example: MnZn ferrite)
    eps_s = 100          # High static permittivity
    eps_inf = 15         # High-frequency permittivity
    tau = 1e-9           # Relaxation time (1 ns)
    alpha = 0.7          # Cole-Cole distribution parameter
    sigma = 0.01         # Low but non-zero conductivity

    # Geometry
    C_0 = 1e-12          # 1 pF geometry factor

    print(f"Material: MnZn Ferrite")
    print(f"eps_s = {eps_s}, eps_inf = {eps_inf}")
    print(f"tau = {tau*1e9:.1f} ns, alpha = {alpha}")
    print(f"C_0 = {C_0*1e12:.1f} pF")

    # Generate Verilog-A
    code = generate_cole_cole_veriloga(
        C_0=C_0,
        eps_s=eps_s,
        eps_inf=eps_inf,
        tau=tau,
        alpha=alpha,
        sigma=sigma,
        filename="demo_cole_cole_cap.va"
    )

    print(f"\nGenerated: demo_cole_cole_cap.va")

    # Compare Cole-Cole vs Debye
    print("\nCole-Cole vs Debye comparison (alpha = 0.7):")
    print(f"{'Frequency':>12} {'CC eps_re':>12} {'Debye eps_re':>12} {'Diff %':>10}")
    print("-" * 50)
    for f in [1e6, 10e6, 100e6, 1e9]:
        omega = 2 * np.pi * f

        # Debye
        omega_tau = omega * tau
        denom = 1 + omega_tau**2
        eps_debye = eps_inf + (eps_s - eps_inf) / denom

        # Cole-Cole
        x = (omega * tau)**alpha
        theta = alpha * np.pi / 2
        denom_re = 1 + x * np.cos(theta)
        denom_im = x * np.sin(theta)
        denom_mag2 = denom_re**2 + denom_im**2
        eps_cc = eps_inf + (eps_s - eps_inf) * denom_re / denom_mag2

        diff = (eps_cc - eps_debye) / eps_debye * 100

        freq_str = f"{f/1e6:.0f} MHz" if f < 1e9 else f"{f/1e9:.0f} GHz"
        print(f"{freq_str:>12} {eps_cc:>12.2f} {eps_debye:>12.2f} {diff:>9.1f}%")


def demo_multi_debye():
    """Demo: Multi-pole Debye model for water."""
    print("\n" + "=" * 60)
    print("Demo 4: Multi-Pole Debye Model")
    print("=" * 60)

    # Water (multiple relaxation processes)
    # Reference: Kaatze, J. Chem. Eng. Data, 1989
    eps_inf = 5.0
    poles = [
        (70.0, 8.3e-12),    # Main relaxation (~17 GHz)
        (4.0, 1.0e-13),     # Fast relaxation (~160 GHz)
    ]

    C_0 = 1e-12  # 1 pF

    print(f"Material: Water at 25C")
    print(f"eps_inf = {eps_inf}")
    print("Relaxation poles:")
    for i, (delta_eps, tau) in enumerate(poles):
        f_relax = 1 / (2 * np.pi * tau)
        print(f"  Pole {i+1}: Delta_eps = {delta_eps}, tau = {tau*1e12:.2f} ps, f = {f_relax/1e9:.1f} GHz")

    # Create generator
    gen = VerilogAGenerator()
    params = MultiDebyeParams(C_0=C_0, eps_inf=eps_inf, poles=poles)
    code = gen.generate_multi_debye_module(params, module_name="demo_multi_debye")
    gen.save_module(code, "demo_multi_debye.va")

    print(f"\nGenerated: demo_multi_debye.va")

    # Show frequency response
    print("\nPermittivity vs frequency:")
    print(f"{'Frequency':>12} {'eps_re':>12} {'eps_im':>12}")
    print("-" * 40)
    for f in [1e8, 1e9, 5e9, 10e9, 50e9, 100e9]:
        omega = 2 * np.pi * f
        eps_re = eps_inf
        eps_im = 0.0
        for delta_eps, tau in poles:
            omega_tau = omega * tau
            denom = 1 + omega_tau**2
            eps_re += delta_eps / denom
            eps_im += delta_eps * omega_tau / denom
        freq_str = f"{f/1e9:.1f} GHz"
        print(f"{freq_str:>12} {eps_re:>12.2f} {eps_im:>12.2f}")


def demo_peec_segment():
    """Demo: Complete PEEC segment with all effects."""
    print("\n" + "=" * 60)
    print("Demo 5: Complete PEEC Segment")
    print("=" * 60)

    # PCB trace parameters
    thickness = 35e-6    # 1 oz copper (35 um)
    width = 0.5e-3       # 0.5 mm trace width
    length = 10e-3       # 10 mm trace length
    sigma = 5.8e7        # Copper
    mu_r = 1.0

    # Substrate (FR4)
    h_sub = 0.2e-3       # 0.2 mm substrate
    eps_r = 4.3
    tan_d = 0.02         # Loss tangent

    # Calculate parameters
    mu0 = 4 * np.pi * 1e-7
    eps0 = 8.854e-12

    R_dc = length / (sigma * width * thickness)
    L_int = mu0 * mu_r * length * thickness / (3 * width)
    L_ext = mu0 * length / (2*np.pi) * (np.log(2*length/width) - 1)  # Approximate
    C_self = eps0 * eps_r * width * length / h_sub

    # Skin effect time constant
    tau_skin = (thickness/2)**2 * mu0 * mu_r * sigma

    # Debye parameters from loss tangent
    # tan_d = eps_im / eps_re at characteristic frequency
    # For FR4, relaxation around 1 GHz
    tau_debye = 1e-10  # 100 ps
    eps_s = eps_r * (1 + tan_d)
    eps_inf = eps_r * (1 - tan_d)

    print(f"PCB Trace: {width*1e3:.2f} mm x {length*1e3:.1f} mm x {thickness*1e6:.0f} um")
    print(f"Substrate: FR4, h = {h_sub*1e3:.2f} mm, eps_r = {eps_r}")
    print(f"\nExtracted parameters:")
    print(f"  R_dc = {R_dc*1e3:.4f} mOhm")
    print(f"  L_int = {L_int*1e12:.3f} pH")
    print(f"  L_ext = {L_ext*1e9:.3f} nH")
    print(f"  C_self = {C_self*1e15:.2f} fF")
    print(f"  tau_skin = {tau_skin*1e12:.3f} ps")

    # Create generator
    gen = VerilogAGenerator()

    debye = DebyeParams(
        C_0=C_self / eps_r,  # Geometry factor
        eps_s=eps_s,
        eps_inf=eps_inf,
        tau=tau_debye
    )

    code = gen.generate_peec_segment_module(
        R_dc=R_dc,
        L_self=L_int,
        C_self=C_self,
        tau_skin=tau_skin,
        debye_params=debye,
        L_ext=L_ext,
        module_name="demo_peec_segment"
    )
    gen.save_module(code, "demo_peec_segment.va")

    print(f"\nGenerated: demo_peec_segment.va")

    # Calculate impedance at various frequencies
    print("\nImpedance vs frequency (approximate):")
    print(f"{'Frequency':>12} {'R_ac':>12} {'X_L':>12} {'|Z|':>12}")
    print("-" * 52)
    for f in [1e3, 10e3, 100e3, 1e6, 10e6, 100e6]:
        omega = 2 * np.pi * f
        delta = np.sqrt(2 / (omega * mu0 * sigma))
        xi = (thickness/2) / delta

        if xi < 0.1:
            F_R = 1.0
            F_L = 1.0
        else:
            F_R = xi * (np.sinh(2*xi) + np.sin(2*xi)) / (np.cosh(2*xi) - np.cos(2*xi))
            F_L = 1.5 / xi * (np.sinh(2*xi) - np.sin(2*xi)) / (np.cosh(2*xi) - np.cos(2*xi))

        R_ac = R_dc * F_R
        X_L = omega * (L_int * F_L + L_ext)
        Z_mag = np.sqrt(R_ac**2 + X_L**2)

        freq_str = f"{f/1e3:.0f} kHz" if f < 1e6 else f"{f/1e6:.0f} MHz"
        print(f"{freq_str:>12} {R_ac*1e3:>10.4f} mOhm {X_L*1e3:>10.4f} mOhm {Z_mag*1e3:>10.4f} mOhm")


def main():
    """Run all demos."""
    print("=" * 60)
    print("Verilog-A Export Demo for PEEC Models")
    print("=" * 60)
    print("\nThis demo generates Verilog-A modules for SPICE simulation.")
    print("Generated files can be used with Spectre, ADS, HSPICE, etc.\n")

    demo_dowell_skin_effect()
    demo_debye_capacitor()
    demo_cole_cole_capacitor()
    demo_multi_debye()
    demo_peec_segment()

    print("\n" + "=" * 60)
    print("Summary: Generated Verilog-A files")
    print("=" * 60)
    print("  1. demo_dowell_skin.va    - Dowell skin effect impedance")
    print("  2. demo_debye_cap.va      - Debye relaxation capacitor")
    print("  3. demo_cole_cole_cap.va  - Cole-Cole relaxation capacitor")
    print("  4. demo_multi_debye.va    - Multi-pole Debye capacitor")
    print("  5. demo_peec_segment.va   - Complete PEEC segment")
    print("\nThese files can be included in SPICE simulations using:")
    print("  .HDLINCLUDE 'filename.va'")
    print("=" * 60)


if __name__ == '__main__':
    main()
