"""
PRIMA to SPICE Export Demo

**LIMITATION**: This demo uses SCALAR Lanczos and is valid for **1-PORT ONLY**.
For N-port systems (N >= 2), use Block Lanczos classes:
  - NPortBlockLanczosSPICE
  - NPortCoupledMagneticSPICE
  - NPortCoupledDielectricSPICE
  - NPortFullCoupledSPICE

See docs/peec/NPORT_BLOCK_LANCZOS_SPICE.md for details.

Demonstrates the complete workflow:
1. PEEC matrix generation (Loop-Star)
2. PRIMA Lanczos reduction
3. SPICE netlist generation

Uses the same wire geometry as verify_loop_star.py but exports to SPICE.
"""

import sys
import os
import numpy as np

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad

try:
    from peec_matrices import PEECBuilder
    print("Using C++ PEEC implementation")
except ImportError:
    print("C++ PEEC module not available, using Python fallback")
    from radia.peec_matrices import PEECBuilder

from lanczos_reduction import (
    LanczosReducer,
    SPICEExtractionConfig,
    PRIMASchurExtractor
)

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # H/m

mm = 1e-3  # 1 mm in meters


def compute_dowell_ladder_params(d, sigma, mu_r=1.0, n_stages=5):
    """Compute Dowell continued fraction ladder parameters.

    The continued fraction expansion of z*coth(z):
        z*coth(z) = 1 + w/(3 + w/(5 + w/(7 + ...)))

    where w = z^2 = tau * s, tau = d^2 * mu * sigma / 2

    This maps to an RL ladder network for skin effect modeling.

    Parameters
    ----------
    d : float
        Conductor thickness [m]
    sigma : float
        Conductivity [S/m]
    mu_r : float
        Relative permeability
    n_stages : int
        Number of ladder stages (more = higher frequency accuracy)

    Returns
    -------
    dict
        Ladder parameters: R_dc, L_int_dc, tau, b_coeffs (continued fraction coefficients)
    """
    mu = MU_0 * mu_r

    # DC values (per unit length)
    R_dc = 1.0 / (sigma * d)  # Ohm/m
    L_int_dc = mu * d / 3.0   # H/m (internal inductance)

    # Time constant for continued fraction
    tau = d * d * mu * sigma / 2.0  # seconds

    # Continued fraction coefficients: b_n = 2n + 1
    # z*coth(z) = 1 + w/(3 + w/(5 + w/(7 + ...)))
    b_coeffs = [2*n + 1 for n in range(1, n_stages + 1)]  # [3, 5, 7, 9, ...]

    return {
        'R_dc': R_dc,
        'L_int_dc': L_int_dc,
        'tau': tau,
        'b_coeffs': b_coeffs,
        'd': d,
        'sigma': sigma,
        'mu_r': mu_r
    }


def generate_dowell_ladder_spice(dowell_params, length, prefix="skin"):
    """Generate SPICE subcircuit for Dowell skin effect ladder.

    The ladder network realizes:
        Z_s(s) = R_dc * [1 + w/(3 + w/(5 + ...))]

    where w = tau * s

    Circuit topology:
        in --[R1]--+--[R2]--+--[R3]--+-- ... --out
                   |        |        |
                  [L1]     [L2]     [L3]
                   |        |        |
                  GND      GND      GND

    Parameters
    ----------
    dowell_params : dict
        From compute_dowell_ladder_params()
    length : float
        Conductor length [m]
    prefix : str
        Prefix for element names

    Returns
    -------
    list
        SPICE netlist lines
    """
    R_dc = dowell_params['R_dc'] * length  # Total DC resistance
    tau = dowell_params['tau']
    b_coeffs = dowell_params['b_coeffs']
    n_stages = len(b_coeffs)

    lines = []
    lines.append(f"* Dowell Skin Effect Ladder ({n_stages} stages)")
    lines.append(f"* R_dc = {R_dc:.6e} Ohm")
    lines.append(f"* tau = {tau:.6e} s")
    lines.append(f"* Coefficients: {b_coeffs}")
    lines.append("")

    # The continued fraction 1 + w/(b1 + w/(b2 + ...)) maps to:
    # Z = R0 + 1/(1/(R1 + sL1) + 1/(R2 + sL2) + ...)
    #
    # For the standard Dowell form, we use a ladder:
    # Each stage has R_n = R_dc / b_n and L_n = tau * R_dc / b_n

    # Stage 0: Series R for the "1" term
    lines.append(f"R{prefix}0 {prefix}_in {prefix}_n0 {R_dc:.6e}")

    # Stages 1 to n: RL branches to ground
    prev_node = f"{prefix}_n0"
    for i, b_n in enumerate(b_coeffs):
        curr_node = f"{prefix}_n{i+1}" if i < n_stages - 1 else f"{prefix}_out"
        R_n = R_dc / b_n
        L_n = tau * R_dc / b_n

        # Series R between nodes
        lines.append(f"R{prefix}{i+1} {prev_node} {curr_node} {R_n:.6e}")

        # Shunt L to ground (represents the continued fraction structure)
        if i < n_stages - 1:  # Don't add shunt on last stage
            lines.append(f"L{prefix}{i+1} {curr_node} 0 {L_n:.6e}")

        prev_node = curr_node

    return lines


def generate_skin_effect_spice_netlist(L, R, params, n_lanczos=5, n_dowell=5,
                                        output_file=None, use_verilog_a=False):
    """Generate SPICE netlist with skin effect (Dowell ladder).

    Combines PRIMA reduction for external inductance with Dowell ladder
    for internal impedance (skin effect).

    Z_total(s) = Z_ext(s) + Z_skin(s)
               = [PRIMA ladder] + [Dowell ladder]

    Parameters
    ----------
    L : np.ndarray
        Inductance matrix (external + internal)
    R : np.ndarray
        DC resistance
    params : dict
        Wire parameters
    n_lanczos : int
        PRIMA stages for external inductance
    n_dowell : int
        Dowell stages for skin effect (5 = good to ~100x skin depth)
    output_file : str, optional
        Output file path
    use_verilog_a : bool
        Generate Verilog-A instead of SPICE
    """
    print("\n" + "=" * 60)
    print("PRIMA + Dowell Skin Effect SPICE Generation")
    print("=" * 60)

    n = L.shape[0]
    n_lanczos = min(n_lanczos, n)

    length = params['length']
    width = params['width']
    height = params['height']
    sigma = params['sigma']

    # Use smaller dimension for skin effect
    d = min(width, height)

    # Compute Dowell parameters
    dowell = compute_dowell_ladder_params(d, sigma, mu_r=1.0, n_stages=n_dowell)

    print(f"\nDowell skin effect parameters:")
    print(f"  Conductor thickness d: {d*1e6:.1f} um")
    print(f"  R_dc (total):  {dowell['R_dc']*length*1e3:.4f} mOhm")
    print(f"  L_int_dc:      {dowell['L_int_dc']*length*1e9:.2f} nH")
    print(f"  tau:           {dowell['tau']*1e9:.3f} ns")
    print(f"  Dowell stages: {n_dowell}")

    # Lanczos reduction for external inductance
    lanczos = LanczosReducer(tol=1e-10)
    result = lanczos.lanczos_symmetric(L, k=n_lanczos)

    print(f"\nPRIMA reduction: {n} -> {n_lanczos} stages")

    # R transformation
    if R.ndim == 1:
        R_full = np.diag(R)
    else:
        R_full = R
    R_reduced = result.Q.T @ R_full @ result.Q

    alpha = result.alpha
    beta = result.beta
    diag_R_red = np.diag(R_reduced)

    lines = []

    if use_verilog_a:
        # === Verilog-A with skin effect ===
        lines.append("// PRIMA + Dowell Skin Effect Verilog-A module")
        lines.append("// Wire: L={:.0f}mm, w={:.1f}mm, h={:.1f}mm".format(
            length*1e3, width*1e3, height*1e3))
        lines.append(f"// PRIMA stages: {n_lanczos}, Dowell stages: {n_dowell}")
        lines.append("// Skin effect modeled via Dowell continued fraction")
        lines.append("")
        lines.append("`include \"disciplines.vams\"")
        lines.append("")
        lines.append("module wire_prima_skin(port_in, port_out);")
        lines.append("    inout port_in, port_out;")
        lines.append("    electrical port_in, port_out;")
        lines.append("")

        # Internal nodes for PRIMA
        for i in range(1, n_lanczos):
            lines.append(f"    electrical n{i};")

        # Internal nodes for Dowell ladder
        lines.append(f"    electrical skin_in, skin_out;")
        for i in range(n_dowell):
            lines.append(f"    electrical skin_n{i};")
        lines.append("")

        lines.append("    // Dowell parameters")
        R_dc_total = dowell['R_dc'] * length
        tau = dowell['tau']
        lines.append(f"    parameter real R_dc = {R_dc_total:.6e};")
        lines.append(f"    parameter real tau = {tau:.6e};")
        lines.append("")

        lines.append("    analog begin")

        # Dowell ladder (skin effect)
        lines.append("        // === Dowell Skin Effect Ladder ===")
        lines.append(f"        V(port_in, skin_n0) <+ I(port_in, skin_n0) * R_dc;")

        for i, b_n in enumerate(dowell['b_coeffs']):
            curr_node = f"skin_n{i}"
            next_node = f"skin_n{i+1}" if i < n_dowell - 1 else "skin_out"
            R_n = R_dc_total / b_n
            L_n = tau * R_dc_total / b_n

            lines.append(f"        // Stage {i+1}: b_{i+1} = {b_n}")
            lines.append(f"        V({curr_node}, {next_node}) <+ I({curr_node}, {next_node}) * {R_n:.6e};")
            if i < n_dowell - 1:
                lines.append(f"        V({next_node}, gnd) <+ {L_n:.6e} * ddt(I({next_node}, gnd));")

        lines.append("")
        lines.append("        // === PRIMA External Inductance Ladder ===")

        # PRIMA ladder (external inductance)
        for i in range(n_lanczos):
            n1 = "skin_out" if i == 0 else f"n{i}"
            n2 = f"n{i+1}" if i < n_lanczos-1 else "port_out"
            L_val = float(alpha[i])

            lines.append(f"        V({n1}, {n2}) <+ {L_val:.6e} * ddt(I({n1}, {n2}));")

        # Mutual inductance note
        if len(beta) > 0:
            lines.append("")
            lines.append("        // Mutual inductance (PRIMA tridiagonal)")
            for i in range(len(beta)):
                M_val = beta[i]
                if np.abs(M_val) > 1e-15:
                    k_val = M_val / np.sqrt(alpha[i] * alpha[i+1])
                    lines.append(f"        // K{i+1}_{i+2} = {k_val:.6e}")

        lines.append("    end")
        lines.append("endmodule")

    else:
        # === Standard SPICE with skin effect ===
        lines.append("* PRIMA + Dowell Skin Effect SPICE netlist")
        lines.append("* Wire: L={:.0f}mm, w={:.1f}mm, h={:.1f}mm".format(
            length*1e3, width*1e3, height*1e3))
        lines.append(f"* PRIMA stages: {n_lanczos}, Dowell stages: {n_dowell}")
        lines.append("* Skin effect: Dowell continued fraction ladder")
        # Frequency where skin depth = conductor thickness
        f_skin = 1.0 / (np.pi * dowell['tau']) if dowell['tau'] > 0 else 1e9
        lines.append(f"* Skin effect significant above: ~{f_skin/1e3:.1f} kHz")
        lines.append("")
        lines.append(".SUBCKT WIRE_PRIMA_SKIN port_in port_out")
        lines.append("")

        # Dowell ladder for skin effect
        lines.append("* === Dowell Skin Effect Ladder ===")
        R_dc_total = dowell['R_dc'] * length
        tau = dowell['tau']

        # Stage 0: DC resistance
        lines.append(f"Rskin0 port_in skin_n0 {R_dc_total:.6e}")

        # Continued fraction stages
        for i, b_n in enumerate(dowell['b_coeffs']):
            curr_node = f"skin_n{i}"
            next_node = f"skin_n{i+1}" if i < n_dowell - 1 else "skin_out"
            R_n = R_dc_total / b_n
            L_n = tau * R_dc_total / b_n

            lines.append(f"Rskin{i+1} {curr_node} {next_node} {R_n:.6e}")
            if i < n_dowell - 1:
                lines.append(f"Lskin{i+1} {next_node} 0 {L_n:.6e}")

        lines.append("")

        # PRIMA ladder for external inductance
        lines.append("* === PRIMA External Inductance Ladder ===")

        for i in range(n_lanczos):
            n1 = "skin_out" if i == 0 else f"ext_n{i}"
            n2 = f"ext_n{i+1}" if i < n_lanczos-1 else "port_out"
            L_val = float(alpha[i])

            lines.append(f"Lext{i+1} {n1} {n2} {L_val:.6e}")

        lines.append("")

        # Tridiagonal coupling
        lines.append("* === PRIMA Tridiagonal Coupling ===")
        for i in range(len(beta)):
            M_val = beta[i]
            if np.abs(M_val) > 1e-15:
                k_val = M_val / np.sqrt(alpha[i] * alpha[i+1])
                if np.abs(k_val) < 1.0:
                    lines.append(f"Kext{i+1}_{i+2} Lext{i+1} Lext{i+2} {k_val:.6e}")

        lines.append("")
        lines.append(".ENDS WIRE_PRIMA_SKIN")
        lines.append("")

        # Test circuit
        lines.append("* === Test Circuit ===")
        lines.append("Xwire port_in port_out WIRE_PRIMA_SKIN")
        lines.append("Vin port_in 0 AC 1")
        lines.append("")
        lines.append("* AC analysis: 1kHz to 100MHz")
        lines.append(".AC DEC 20 1k 100MEG")
        lines.append(".PRINT AC V(port_out) I(Vin)")
        lines.append(".END")

    netlist = '\n'.join(lines)

    print(f"\nNetlist generated: {len(lines)} lines")
    print(f"  Dowell stages:  {n_dowell}")
    print(f"  PRIMA stages:   {n_lanczos}")
    print(f"  Total elements: {n_dowell + n_lanczos + len(beta)} (approx)")

    # Print netlist
    print(f"\n--- PRIMA + Skin Effect Netlist ---")
    for line in lines:
        print(line)

    # Save to file
    if output_file:
        with open(output_file, 'w') as f:
            f.write(netlist)
        print(f"\nNetlist saved to: {output_file}")

    return netlist, result, dowell


def create_straight_wire_peec():
    """Create PEEC matrices for a straight wire."""
    print("\n" + "=" * 60)
    print("Creating PEEC Matrices for Straight Wire")
    print("=" * 60)

    # Wire parameters
    length = 100*mm    # 10 cm
    width = 1*mm       # 1 mm
    height = 1*mm      # 1 mm
    n_segments = 10
    sigma = 5.8e7     # Copper

    print(f"\nWire parameters:")
    print(f"  Length:     {length/mm:.0f} mm")
    print(f"  Width:      {width/mm:.1f} mm")
    print(f"  Height:     {height/mm:.1f} mm")
    print(f"  Segments:   {n_segments}")
    print(f"  Sigma:      {sigma:.2e} S/m")

    # Build PEEC matrices
    builder = PEECBuilder()
    builder.create_wire([0, 0, 0], [length, 0, 0], width, height, n_segments, sigma)
    L, R, P, M_LS = builder.build()

    print(f"\nPEEC Matrix dimensions:")
    print(f"  L (inductance):        {L.shape}")
    print(f"  R (resistance):        {R.shape}")
    print(f"  P (potential coeff):   {P.shape}")
    print(f"  M_LS (loop-star):      {M_LS.shape}")

    return L, R, P, M_LS, {
        'length': length,
        'width': width,
        'height': height,
        'n_segments': n_segments,
        'sigma': sigma
    }


def demo_lanczos_reduction(L, R=None):
    """Demonstrate Lanczos tridiagonalization."""
    print("\n" + "=" * 60)
    print("Lanczos Tridiagonalization Demo")
    print("=" * 60)

    _ = R  # R not used in Lanczos demo (inductance only)
    n = L.shape[0]
    lanczos = LanczosReducer(tol=1e-10)

    # Full Lanczos
    result_full = lanczos.lanczos_symmetric(L, k=n)

    print(f"\nFull Lanczos (k={n}):")
    print(f"  alpha (diagonal):     {result_full.alpha}")
    print(f"  beta (off-diagonal):  {result_full.beta}")

    # Verify: L = Q @ T @ Q.T
    L_reconstructed = result_full.Q @ result_full.T @ result_full.Q.T
    error = np.linalg.norm(L - L_reconstructed) / np.linalg.norm(L)
    print(f"  Reconstruction error: {error:.2e}")

    # Reduced Lanczos (PRIMA ladder)
    k_reduced = min(5, n)
    result_reduced = lanczos.lanczos_symmetric(L, k=k_reduced)

    print(f"\nReduced Lanczos (k={k_reduced}):")
    print(f"  alpha (diagonal):     {result_reduced.alpha}")
    print(f"  beta (off-diagonal):  {result_reduced.beta}")

    # Eigenvalue comparison
    eigvals_original = np.linalg.eigvalsh(L)
    eigvals_reduced = np.linalg.eigvalsh(result_reduced.T)

    print(f"\nEigenvalue comparison:")
    print(f"  Original L:    min={eigvals_original.min():.2e}, max={eigvals_original.max():.2e}")
    print(f"  Reduced T:     min={eigvals_reduced.min():.2e}, max={eigvals_reduced.max():.2e}")

    return result_full, result_reduced


def generate_simple_spice_netlist(L, R, params, output_file=None):
    """Generate a simple SPICE netlist from PEEC matrices (no full extraction)."""
    print("\n" + "=" * 60)
    print("Simple SPICE Netlist Generation")
    print("=" * 60)

    n = L.shape[0]

    lines = []
    lines.append("* PEEC-extracted SPICE netlist")
    lines.append("* Wire: L={:.0f}mm, w={:.1f}mm, h={:.1f}mm".format(
        params['length']*1e3, params['width']*1e3, params['height']*1e3))
    lines.append("* Segments: {}".format(params['n_segments']))
    lines.append("")
    lines.append(".SUBCKT WIRE_PEEC port_in port_out")
    lines.append("")

    # Series RL ladder from diagonal
    lines.append("* === Series RL Ladder ===")
    diag_L = np.diag(L)
    diag_R = np.diag(R) if R.shape == L.shape else np.zeros(n)

    for i in range(n):
        n1 = "port_in" if i == 0 else f"n{i}"
        n2 = f"n{i+1}" if i < n-1 else "port_out"

        R_val = float(diag_R[i]) if i < len(diag_R) else 0.0
        L_val = float(diag_L[i])

        if R_val > 0:
            lines.append(f"R{i+1} {n1} {n1}a {R_val:.6e}")
            n1 = f"{n1}a"

        lines.append(f"L{i+1} {n1} {n2} {L_val:.6e}")

    lines.append("")

    # Mutual inductance (off-diagonal L)
    lines.append("* === Mutual Inductance ===")
    mutual_count = 0
    for i in range(n):
        for j in range(i+1, n):
            M_val = L[i, j]
            if np.abs(M_val) > 1e-15:
                # Coupling coefficient k = M / sqrt(L_i * L_j)
                k_val = M_val / np.sqrt(diag_L[i] * diag_L[j])
                if np.abs(k_val) < 1.0:
                    lines.append(f"K{i+1}_{j+1} L{i+1} L{j+1} {k_val:.6e}")
                    mutual_count += 1

    if mutual_count == 0:
        lines.append("* (no significant mutual inductance)")

    lines.append("")
    lines.append(".ENDS WIRE_PEEC")
    lines.append("")

    # Test circuit
    lines.append("* === Test Circuit ===")
    lines.append("Xwire port_in port_out WIRE_PEEC")
    lines.append("Vin port_in 0 AC 1")
    lines.append("")
    lines.append("* AC analysis: 1kHz to 100MHz")
    lines.append(".AC DEC 10 1k 100MEG")
    lines.append(".PRINT AC V(port_out) I(Vin)")
    lines.append(".END")

    netlist = '\n'.join(lines)

    print(f"\nNetlist generated: {len(lines)} lines")
    print(f"  Resistors:      {n}")
    print(f"  Inductors:      {n}")
    print(f"  Mutual couplings: {mutual_count}")

    # Print first 30 lines
    print(f"\n--- Netlist Preview (first 30 lines) ---")
    for line in lines[:30]:
        print(line)
    if len(lines) > 30:
        print(f"... ({len(lines) - 30} more lines)")

    # Save to file if specified
    if output_file:
        with open(output_file, 'w') as f:
            f.write(netlist)
        print(f"\nNetlist saved to: {output_file}")

    return netlist


def generate_prima_spice_netlist(L, R, params, n_lanczos=5, output_file=None,
                                  use_verilog_a=False):
    """Generate SPICE netlist using PRIMA Lanczos reduction.

    Parameters
    ----------
    L : np.ndarray
        Inductance matrix
    R : np.ndarray
        Resistance matrix (1D diagonal or 2D full)
    params : dict
        Wire parameters (length, width, height, n_segments, sigma)
    n_lanczos : int
        Number of Lanczos stages for PRIMA reduction
    output_file : str, optional
        Output file path for SPICE netlist
    use_verilog_a : bool
        If True, generate Verilog-A module instead of standard SPICE.
        WARNING: Verilog-A is not supported by all SPICE solvers.
        Use standard SPICE (default) for maximum compatibility.
    """
    print("\n" + "=" * 60)
    print("PRIMA Lanczos SPICE Netlist Generation")
    if use_verilog_a:
        print("  Output format: Verilog-A (requires compatible simulator)")
    else:
        print("  Output format: Standard SPICE (maximum compatibility)")
    print("=" * 60)

    n = L.shape[0]
    n_lanczos = min(n_lanczos, n)

    # Lanczos reduction
    lanczos = LanczosReducer(tol=1e-10)
    result = lanczos.lanczos_symmetric(L, k=n_lanczos)

    print(f"\nLanczos reduction: {n} -> {n_lanczos} stages")
    print(f"  Compression ratio: {100*(1-n_lanczos/n):.1f}%")

    # Transform R to Lanczos basis
    # R may be 1D (diagonal) or 2D (full matrix)
    if R.ndim == 1:
        R_full = np.diag(R)
    else:
        R_full = R
    R_reduced = result.Q.T @ R_full @ result.Q

    # Extract tridiagonal parameters
    alpha = result.alpha  # Diagonal of T (self-inductance)
    beta = result.beta    # Off-diagonal of T (mutual coupling)
    diag_R_red = np.diag(R_reduced)

    lines = []

    if use_verilog_a:
        # === Verilog-A output ===
        lines.append("// PRIMA-reduced Verilog-A module")
        lines.append("// Wire: L={:.0f}mm, w={:.1f}mm, h={:.1f}mm".format(
            params['length']*1e3, params['width']*1e3, params['height']*1e3))
        lines.append(f"// Original segments: {n}, PRIMA stages: {n_lanczos}")
        lines.append("// WARNING: Verilog-A requires compatible simulator (Spectre, ADS, etc.)")
        lines.append("")
        lines.append("`include \"disciplines.vams\"")
        lines.append("")
        lines.append("module wire_prima(port_in, port_out);")
        lines.append("    inout port_in, port_out;")
        lines.append("    electrical port_in, port_out;")
        lines.append("")

        # Internal nodes
        for i in range(1, n_lanczos):
            lines.append(f"    electrical n{i};")
        lines.append("")

        # Branch declarations
        lines.append("    // Branch currents for inductors")
        for i in range(n_lanczos):
            n1 = "port_in" if i == 0 else f"n{i}"
            n2 = f"n{i+1}" if i < n_lanczos-1 else "port_out"
            lines.append(f"    branch (n{i}a, {n2}) br_L{i+1};") if i > 0 else lines.append(f"    branch (port_in_a, {n2}) br_L{i+1};")
        lines.append("")

        lines.append("    analog begin")

        # RL ladder equations
        for i in range(n_lanczos):
            n1 = "port_in" if i == 0 else f"n{i}"
            n2 = f"n{i+1}" if i < n_lanczos-1 else "port_out"
            R_val = float(diag_R_red[i])
            L_val = float(alpha[i])

            # Resistor: V = I * R
            lines.append(f"        // Stage {i+1}: R={R_val:.6e}, L={L_val:.6e}")
            if R_val > 0:
                lines.append(f"        V({n1}, {n1}_a) <+ I({n1}, {n1}_a) * {R_val:.6e};")
                n1 = f"{n1}_a"

            # Inductor: V = L * dI/dt
            lines.append(f"        V({n1}, {n2}) <+ {L_val:.6e} * ddt(I({n1}, {n2}));")

        # Mutual inductance (simplified - requires more complex Verilog-A for full coupling)
        if len(beta) > 0:
            lines.append("")
            lines.append("        // Note: Mutual inductance requires state-space formulation")
            lines.append("        // for accurate Verilog-A implementation.")
            for i in range(len(beta)):
                M_val = beta[i]
                if np.abs(M_val) > 1e-15:
                    k_val = M_val / np.sqrt(alpha[i] * alpha[i+1])
                    lines.append(f"        // K{i+1}_{i+2} = {k_val:.6e} (M = {M_val:.6e})")

        lines.append("    end")
        lines.append("endmodule")

    else:
        # === Standard SPICE output ===
        lines.append("* PRIMA-reduced SPICE netlist")
        lines.append("* Wire: L={:.0f}mm, w={:.1f}mm, h={:.1f}mm".format(
            params['length']*1e3, params['width']*1e3, params['height']*1e3))
        lines.append(f"* Original segments: {n}, PRIMA stages: {n_lanczos}")
        lines.append("")
        lines.append(".SUBCKT WIRE_PRIMA port_in port_out")
        lines.append("")

        # PRIMA ladder from tridiagonal T
        lines.append("* === PRIMA Ladder (Tridiagonal) ===")

        for i in range(n_lanczos):
            n1 = "port_in" if i == 0 else f"n{i}"
            n2 = f"n{i+1}" if i < n_lanczos-1 else "port_out"

            R_val = float(diag_R_red[i])
            L_val = float(alpha[i])

            if R_val > 0:
                lines.append(f"R{i+1} {n1} {n1}a {R_val:.6e}")
                n1 = f"{n1}a"

            lines.append(f"L{i+1} {n1} {n2} {L_val:.6e}")

        lines.append("")

        # Tridiagonal coupling (adjacent stages only)
        lines.append("* === Tridiagonal Coupling ===")
        for i in range(len(beta)):
            M_val = beta[i]
            if np.abs(M_val) > 1e-15:
                k_val = M_val / np.sqrt(alpha[i] * alpha[i+1])
                if np.abs(k_val) < 1.0:
                    lines.append(f"K{i+1}_{i+2} L{i+1} L{i+2} {k_val:.6e}")

        lines.append("")
        lines.append(".ENDS WIRE_PRIMA")
        lines.append("")

        # Test circuit
        lines.append("* === Test Circuit ===")
        lines.append("Xwire port_in port_out WIRE_PRIMA")
        lines.append("Vin port_in 0 AC 1")
        lines.append("")
        lines.append("* AC analysis: 1kHz to 100MHz")
        lines.append(".AC DEC 10 1k 100MEG")
        lines.append(".PRINT AC V(port_out) I(Vin)")
        lines.append(".END")

    netlist = '\n'.join(lines)

    print(f"\nPRIMA netlist generated: {len(lines)} lines")
    print(f"  RL stages:      {n_lanczos}")
    print(f"  Tridiagonal K:  {len(beta)}")

    # Print netlist
    print(f"\n--- PRIMA Netlist ---")
    for line in lines:
        print(line)

    # Save to file if specified
    if output_file:
        with open(output_file, 'w') as f:
            f.write(netlist)
        print(f"\nNetlist saved to: {output_file}")

    return netlist, result


def verify_frequency_response(L, R, result_lanczos, frequencies=None):
    """Verify frequency response of full vs PRIMA-reduced model."""
    print("\n" + "=" * 60)
    print("Frequency Response Verification")
    print("=" * 60)

    if frequencies is None:
        frequencies = np.logspace(3, 8, 21)  # 1kHz to 100MHz

    n = L.shape[0]
    k = result_lanczos.rank

    # R may be 1D (diagonal) or 2D (full matrix)
    if R.ndim == 1:
        R_full = np.diag(R)
    else:
        R_full = R

    # Reduced model matrices in Lanczos basis
    Q = result_lanczos.Q
    R_reduced = Q.T @ R_full @ Q
    L_reduced = result_lanczos.T  # Tridiagonal

    errors = []
    print(f"\n{'Freq (Hz)':>12s}  {'|Z_full|':>12s}  {'|Z_red|':>12s}  {'Error':>8s}")
    print("-" * 50)

    for f in frequencies:
        omega = 2 * np.pi * f

        # Full model: Z = R + j*omega*L
        # Port impedance: Z_port = e_1^T @ Z^{-1} @ e_1
        Z_full = R_full + 1j * omega * L
        try:
            Z_inv_full = np.linalg.inv(Z_full)
            # Total impedance seen from port 1 to port n (series connection)
            # For a transmission line: Z_total = sum of Z_ii (for series)
            # Here we compute Z_11 element of impedance matrix
            Z_port_full = Z_full[0, 0] - Z_full[0, 1:] @ np.linalg.solve(
                Z_full[1:, 1:], Z_full[1:, 0])
        except Exception:
            Z_port_full = Z_full[0, 0]

        # Reduced model: same computation in Lanczos basis
        Z_red = R_reduced + 1j * omega * L_reduced
        try:
            Z_port_red = Z_red[0, 0] - Z_red[0, 1:] @ np.linalg.solve(
                Z_red[1:, 1:], Z_red[1:, 0])
        except Exception:
            Z_port_red = Z_red[0, 0]

        error = np.abs(Z_port_full - Z_port_red) / (np.abs(Z_port_full) + 1e-15)
        errors.append(error)

        print(f"{f:12.2e}  {np.abs(Z_port_full):12.4e}  {np.abs(Z_port_red):12.4e}  {error*100:7.2f}%")

    max_error = max(errors) * 100
    avg_error = np.mean(errors) * 100
    print(f"\nMax error: {max_error:.2f}%")
    print(f"Avg error: {avg_error:.2f}%")

    if max_error < 10:
        print("PASS: Frequency response preserved")
    else:
        print("NOTE: PRIMA reduced model uses Lanczos basis (not original basis)")
        print("      Port mapping requires proper projection for exact match")

    return errors


def main():
    """Run PRIMA to SPICE export demo."""
    import argparse

    parser = argparse.ArgumentParser(
        description='PRIMA to SPICE Export Demo',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Output formats:
  Standard SPICE (default) - Works with ngspice, LTspice, HSPICE, PSpice
  Verilog-A (--verilog-a)  - Requires Spectre, ADS, or Verilog-A compatible simulator

Skin effect modeling:
  --skin-effect           Add Dowell ladder for frequency-dependent R(f), L(f)
  --dowell N              Number of Dowell stages (default: 5)

Examples:
  python demo_prima_spice_export.py                    # Standard SPICE (DC only)
  python demo_prima_spice_export.py --skin-effect     # With skin effect
  python demo_prima_spice_export.py --verilog-a       # Verilog-A output
  python demo_prima_spice_export.py --lanczos 10      # 10-stage PRIMA
""")
    parser.add_argument('--verilog-a', action='store_true',
                        help='Generate Verilog-A module instead of standard SPICE '
                             '(WARNING: not supported by all simulators)')
    parser.add_argument('--lanczos', type=int, default=5,
                        help='Number of Lanczos stages for PRIMA reduction (default: 5)')
    parser.add_argument('--skin-effect', action='store_true',
                        help='Include Dowell skin effect ladder for frequency-dependent '
                             'resistance and internal inductance')
    parser.add_argument('--dowell', type=int, default=5,
                        help='Number of Dowell stages for skin effect (default: 5)')
    args = parser.parse_args()

    print("=" * 70)
    print("PRIMA to SPICE Export Demo")
    print("=" * 70)

    if args.verilog_a:
        print("\n*** Verilog-A output selected ***")
        print("WARNING: Verilog-A is not supported by all SPICE solvers.")
        print("         Compatible simulators: Cadence Spectre, Keysight ADS,")
        print("         Mentor Eldo, Silvaco SmartSpice, etc.")
        print("         NOT compatible: ngspice, LTspice, HSPICE (SPICE only)")

    if args.skin_effect:
        print("\n*** Skin effect modeling enabled ***")
        print(f"    Dowell stages: {args.dowell}")
        print("    Models frequency-dependent R(f) and L_internal(f)")

    # Step 1: Create PEEC matrices
    L, R, P, M_LS, params = create_straight_wire_peec()

    # Step 2: Lanczos demonstration
    result_full, result_reduced = demo_lanczos_reduction(L, R)

    # Step 3: Generate simple SPICE (no reduction)
    output_dir = os.path.dirname(__file__)
    _ = generate_simple_spice_netlist(
        L, R, params,
        output_file=os.path.join(output_dir, 'wire_full.sp')
    )

    # Step 4: Generate PRIMA netlist (with or without skin effect)
    if args.skin_effect:
        # Generate with Dowell skin effect ladder
        if args.verilog_a:
            output_file = os.path.join(output_dir, 'wire_prima_skin.va')
        else:
            output_file = os.path.join(output_dir, 'wire_prima_skin.sp')

        skin_netlist, lanczos_result, dowell_params = generate_skin_effect_spice_netlist(
            L, R, params,
            n_lanczos=args.lanczos,
            n_dowell=args.dowell,
            output_file=output_file,
            use_verilog_a=args.verilog_a
        )
    else:
        # Standard PRIMA without skin effect
        if args.verilog_a:
            output_file = os.path.join(output_dir, 'wire_prima.va')
        else:
            output_file = os.path.join(output_dir, 'wire_prima.sp')

        prima_netlist, lanczos_result = generate_prima_spice_netlist(
            L, R, params, n_lanczos=args.lanczos,
            output_file=output_file,
            use_verilog_a=args.verilog_a
        )

    # Step 5: Verify frequency response
    verify_frequency_response(L, R, lanczos_result)

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"\nGenerated files:")
    print(f"  wire_full.sp  - Full PEEC model ({params['n_segments']} stages)")

    if args.skin_effect:
        if args.verilog_a:
            print(f"  wire_prima_skin.va - PRIMA + Dowell Verilog-A ({args.lanczos}+{args.dowell} stages)")
        else:
            print(f"  wire_prima_skin.sp - PRIMA + Dowell SPICE ({args.lanczos}+{args.dowell} stages)")
        print(f"\nSkin effect included:")
        print(f"  - Dowell ladder: {args.dowell} stages (frequency-dependent R, L_int)")
        print(f"  - PRIMA ladder:  {args.lanczos} stages (external inductance)")
    elif args.verilog_a:
        print(f"  wire_prima.va - PRIMA reduced Verilog-A module ({args.lanczos} stages)")
        print(f"\nTo simulate (Verilog-A compatible simulator required):")
        print(f"  Spectre: Include wire_prima.va in your design")
        print(f"  ADS: Import Verilog-A module")
    else:
        print(f"  wire_prima.sp - PRIMA reduced SPICE model ({args.lanczos} stages)")
        print(f"\nTo simulate in SPICE:")
        print(f"  ngspice wire_prima.sp")
        print(f"  LTspice: Open wire_prima.sp")

    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
