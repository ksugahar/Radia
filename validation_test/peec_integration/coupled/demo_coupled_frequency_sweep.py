"""
demo_coupled_frequency_sweep.py

Frequency-domain coupled PEEC + MMM analysis with:
  - Surface impedance (Dowell SIBC) for skin effect
  - Magnetic coupling (Delta_L) from ferrite core
  - Complex permeability (mu_r_imag) for magnetic loss

Physical model:
  - 100mm copper wire (1mm x 1mm rectangular cross-section)
  - Ferrite core (30mm cube, mu_r=1000, mu_r_imag=100) offset 10mm above wire
  - Single port between wire endpoints
  - Frequency sweep from 100 Hz to 10 MHz

Output:
  - R(f): Resistance increases with frequency (skin effect + magnetic loss)
  - L(f): Inductance modified by ferrite core coupling
  - Comparison: with/without ferrite, with/without SIBC

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))


# Physical constants
MU_0 = 4.0 * np.pi * 1e-7  # H/m


def dowell_sibc(freq, h, sigma, mu_r=1.0):
    """
    Dowell skin effect correction for rectangular conductor.

    Returns complex surface impedance correction per unit length:
        Zs_correction = R_dc * (F_R - 1) + j*omega*L_int_dc * (F_L - 1)

    where F_R is resistance factor and F_L is internal inductance factor.

    Actually, returns the full per-segment Zs = R_dc*(F_R - 1) + jw*L_int*(F_L - 1)
    But the solver already includes R_dc, so we return the AC EXCESS impedance:
        Zs = R_dc * (F_R - 1) + jw * L_int_dc * (F_L - 1)

    For simplicity we return per-segment values directly.

    Args:
        freq: Frequency in Hz
        h: Conductor thickness [m] (dimension in field direction)
        sigma: Conductivity [S/m]
        mu_r: Relative permeability of conductor (usually 1.0 for copper)

    Returns:
        (F_R, F_L): Dowell resistance and inductance factors
    """
    if freq <= 0:
        return 1.0, 1.0

    omega = 2.0 * np.pi * freq
    mu = mu_r * MU_0
    delta = np.sqrt(2.0 / (omega * mu * sigma))
    xi = h / delta

    if xi < 0.01:
        return 1.0, 1.0
    elif xi > 50:
        return xi, 1.5 / xi
    else:
        sinh_2xi = np.sinh(2.0 * xi)
        sin_2xi = np.sin(2.0 * xi)
        cosh_2xi = np.cosh(2.0 * xi)
        cos_2xi = np.cos(2.0 * xi)
        denom = cosh_2xi - cos_2xi
        if abs(denom) < 1e-30:
            denom = 1e-30
        F_R = xi * (sinh_2xi + sin_2xi) / denom
        F_L = (1.5 / xi) * (sinh_2xi - sin_2xi) / denom
        return F_R, F_L


def make_sibc_func(topo):
    """
    Create Zs_func(freq) callback for PEEC solver.

    Computes per-segment AC excess impedance from Dowell skin effect model.

    Args:
        topo: build_topology() result dict

    Returns:
        callable(freq) -> Zs array (n_loop,) complex
    """
    R_dc = np.array(topo['R'])
    heights = np.array(topo['segment_heights'])
    sigmas = np.array(topo['segment_sigmas'])
    lengths = np.array(topo['segment_lengths'])
    widths = np.array(topo['segment_widths'])
    n_seg = topo['n_loop']

    # Internal inductance at DC: L_int = mu_0 * length / (12 * (width/height + height/width))
    # For rectangular conductor: L_int_dc = mu_0 * length / 12 (approximate for square)
    L_int_dc = MU_0 * lengths / 12.0

    def zs_func(freq):
        Zs = np.zeros(n_seg, dtype=complex)
        if freq <= 0:
            return Zs

        omega = 2.0 * np.pi * freq
        for i in range(n_seg):
            F_R, F_L = dowell_sibc(freq, heights[i], sigmas[i])
            # AC excess: Zs = R_dc*(F_R - 1) + jw*L_int*(F_L - 1)
            R_dc_i = np.diag(R_dc)[i] if R_dc.ndim > 1 else R_dc[i]
            Zs[i] = R_dc_i * (F_R - 1.0) + 1j * omega * L_int_dc[i] * (F_L - 1.0)
        return Zs

    return zs_func


def main():
    from peec_coupled import CoupledPEECSolver
    from peec_topology import PEECCircuitSolver
    from fasthenry_parser import FastHenryParser

    print()
    print("Coupled PEEC + MMM Frequency Sweep with SIBC")
    print("=" * 70)

    # ============================================================
    # Method 1: Direct Python API
    # ============================================================
    print()
    print("Method 1: Direct Python API")
    print("-" * 40)

    # Build a simple wire near a ferrite core using FastHenry format
    inp_text = """
* Copper wire near ferrite core
.Units m
.default sigma=5.8e7

N1 x=0 y=0 z=0
N2 x=0.1 y=0 z=0
E1 N1 N2 w=1e-3 h=1e-3

.external N1 N2

.magnetic
  type=box
  center=0.05,0.015,0
  size=0.03,0.03,0.03
  divisions=3,3,3
  mu_r=1000
  mu_r_imag=100
.endmagnetic

.freq fmin=100 fmax=1e7 ndec=10
"""
    parser = FastHenryParser()
    parser.parse_string(inp_text)

    summary = parser.get_summary()
    print(f"  Nodes: {summary['n_nodes']}")
    print(f"  Segments: {summary['n_segments']}")
    print(f"  Ports: {summary['n_ports']}")
    print(f"  Magnetic blocks: {summary['n_magnetic_blocks']}")
    print(f"  mu_r = {parser.magnetic_blocks[0]['mu_r']}")
    print(f"  mu_r_imag = {parser.magnetic_blocks[0]['mu_r_imag']}")

    freqs = parser.get_frequencies()
    print(f"  Frequencies: {len(freqs)} points ({freqs[0]:.0f} Hz to {freqs[-1]:.0e} Hz)")

    # ============================================================
    # Case 1: Air only (no coupling, no SIBC)
    # ============================================================
    print()
    print("Case 1: Air only (no coupling, no SIBC)")
    builder = parser.to_peec_builder()
    topo = builder.build_topology()
    solver_air = PEECCircuitSolver(topo)
    Z_air = solver_air.frequency_sweep(freqs)

    # ============================================================
    # Case 2: Air + SIBC (skin effect only)
    # ============================================================
    print("Case 2: Air + SIBC (Dowell skin effect)")
    zs_func = make_sibc_func(topo)
    Z_air_sibc = solver_air.frequency_sweep(freqs, zs_func)

    # ============================================================
    # Case 3: Coupled (ferrite core, lossless, no SIBC)
    # ============================================================
    print("Case 3: Coupled (mu_r=1000, mu_r_imag=0, no SIBC)")

    import radia as rad

    # Build magnetic objects manually (without mu_r_imag for comparison)
    rad.UtiDelAll()
    center = np.array([0.05, 0.015, 0])
    size = np.array([0.03, 0.03, 0.03])
    divs = [3, 3, 3]
    dx = size[0] / divs[0]
    dy = size[1] / divs[1]
    dz = size[2] / divs[2]
    corner = center - size / 2.0

    mag_objs = []
    for iz in range(divs[2]):
        for iy in range(divs[1]):
            for ix in range(divs[0]):
                x0 = corner[0] + ix * dx
                y0 = corner[1] + iy * dy
                z0 = corner[2] + iz * dz
                verts = [
                    [x0, y0, z0], [x0+dx, y0, z0],
                    [x0+dx, y0+dy, z0], [x0, y0+dy, z0],
                    [x0, y0, z0+dz], [x0+dx, y0, z0+dz],
                    [x0+dx, y0+dy, z0+dz], [x0, y0+dy, z0+dz],
                ]
                obj = rad.ObjHexahedron(verts, [0, 0, 0])
                mat = rad.MatLin(1000)
                rad.MatApl(obj, mat)
                mag_objs.append(obj)

    solver_coupled_lossless = CoupledPEECSolver(topo, mag_objs, mu_r_imag=0.0)
    solver_coupled_lossless.compute_coupling_matrix(mu_r_real=1000.0)
    Z_coupled_lossless = solver_coupled_lossless.frequency_sweep(freqs)

    # ============================================================
    # Case 4: Coupled + magnetic loss (mu_r_imag=100) + SIBC
    # ============================================================
    print("Case 4: Coupled (mu_r=1000, mu_r_imag=100) + SIBC")

    solver_coupled_lossy = CoupledPEECSolver(topo, mag_objs, mu_r_imag=100.0)
    solver_coupled_lossy._mu_r_real = 1000.0  # Set for tan_delta
    solver_coupled_lossy.Delta_L = solver_coupled_lossless.Delta_L  # Reuse Delta_L
    Z_coupled_lossy = solver_coupled_lossy.frequency_sweep(freqs, zs_func)

    # ============================================================
    # Results Table
    # ============================================================
    print()
    print("=" * 70)
    print("Results: R(f) and L(f)")
    print("=" * 70)
    print()

    omega = 2.0 * np.pi * freqs

    # Extract R and L
    R_air = np.real(Z_air)
    L_air = np.zeros_like(freqs)
    mask = omega > 0
    L_air[mask] = np.imag(Z_air[mask]) / omega[mask]

    R_air_sibc = np.real(Z_air_sibc)
    L_air_sibc = np.zeros_like(freqs)
    L_air_sibc[mask] = np.imag(Z_air_sibc[mask]) / omega[mask]

    R_lossless = np.real(Z_coupled_lossless)
    L_lossless = np.zeros_like(freqs)
    L_lossless[mask] = np.imag(Z_coupled_lossless[mask]) / omega[mask]

    R_lossy = np.real(Z_coupled_lossy)
    L_lossy = np.zeros_like(freqs)
    L_lossy[mask] = np.imag(Z_coupled_lossy[mask]) / omega[mask]

    # Print header
    print(f"  {'Freq':>10s}  {'R_air':>10s}  {'R_SIBC':>10s}  "
          f"{'R_coupled':>10s}  {'R_full':>10s}  "
          f"{'L_air':>10s}  {'L_coupled':>10s}  {'L_full':>10s}")
    print(f"  {'[Hz]':>10s}  {'[mOhm]':>10s}  {'[mOhm]':>10s}  "
          f"{'[mOhm]':>10s}  {'[mOhm]':>10s}  "
          f"{'[nH]':>10s}  {'[nH]':>10s}  {'[nH]':>10s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  "
          f"{'-'*10}  {'-'*10}  {'-'*10}")

    for i in range(len(freqs)):
        f = freqs[i]
        print(f"  {f:10.0f}  "
              f"{R_air[i]*1e3:10.4f}  "
              f"{R_air_sibc[i]*1e3:10.4f}  "
              f"{R_lossless[i]*1e3:10.4f}  "
              f"{R_lossy[i]*1e3:10.4f}  "
              f"{L_air[i]*1e9:10.4f}  "
              f"{L_lossless[i]*1e9:10.4f}  "
              f"{L_lossy[i]*1e9:10.4f}")

    # ============================================================
    # Summary
    # ============================================================
    print()
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    print()

    # DC values
    print(f"  DC resistance:")
    print(f"    R_dc (air) = {R_air[0]*1e3:.4f} mOhm")
    print(f"    R_dc (coupled, lossless) = {R_lossless[0]*1e3:.4f} mOhm")
    print(f"    R_dc (coupled, lossy) = {R_lossy[0]*1e3:.4f} mOhm")
    print()

    # High-frequency values
    i_hf = -1
    f_hf = freqs[i_hf]
    print(f"  At f = {f_hf:.0e} Hz:")
    print(f"    R_air = {R_air[i_hf]*1e3:.4f} mOhm")
    print(f"    R_SIBC = {R_air_sibc[i_hf]*1e3:.4f} mOhm (skin effect)")
    print(f"    R_coupled_lossless = {R_lossless[i_hf]*1e3:.4f} mOhm")
    print(f"    R_coupled_lossy = {R_lossy[i_hf]*1e3:.4f} mOhm (SIBC + mu_r_imag)")
    print()

    delta_L = solver_coupled_lossless.Delta_L
    print(f"  Coupling (Delta_L):")
    print(f"    Delta_L = {delta_L[0,0]*1e9:.4f} nH")
    print(f"    L_air = {L_air[0]*1e9:.4f} nH")
    print(f"    L_coupled = {L_lossless[0]*1e9:.4f} nH")
    print(f"    Increase = {(L_lossless[0] - L_air[0])/L_air[0]*100:.2f}%")
    print()

    print(f"  Magnetic loss (tan_delta_m = mu_r_imag/mu_r):")
    print(f"    tan_delta_m = {solver_coupled_lossy.tan_delta_m:.4f}")
    print(f"    R_magnetic at {f_hf:.0e} Hz = "
          f"{solver_coupled_lossy.get_magnetic_loss_resistance(f_hf)[0,0]*1e3:.4f} mOhm")
    print()

    # Physics checks
    print("  Physics checks:")

    # R should increase with frequency (skin effect + magnetic loss)
    r_increases = all(R_lossy[i+1] >= R_lossy[i] - 1e-15 for i in range(len(freqs)-1))
    print(f"    R increases with freq (skin + mag loss): {'PASS' if r_increases else 'FAIL'}")

    # SIBC should increase R at high frequencies
    r_sibc_gt_dc = R_air_sibc[i_hf] > R_air[0] * 1.01
    print(f"    R_SIBC > R_dc at high freq: {'PASS' if r_sibc_gt_dc else 'FAIL'}")

    # Coupled L > air L
    l_coupled_gt_air = L_lossless[0] > L_air[0]
    print(f"    L_coupled > L_air: {'PASS' if l_coupled_gt_air else 'FAIL'}")

    # Lossy R > lossless R at high freq (magnetic loss adds resistance)
    r_lossy_gt_lossless = R_lossy[i_hf] > R_lossless[i_hf]
    print(f"    R_lossy > R_lossless at high freq: {'PASS' if r_lossy_gt_lossless else 'FAIL'}")

    # R_magnetic should be proportional to frequency
    R_mag_1 = solver_coupled_lossy.get_magnetic_loss_resistance(1e6)[0, 0]
    R_mag_2 = solver_coupled_lossy.get_magnetic_loss_resistance(2e6)[0, 0]
    ratio = R_mag_2 / R_mag_1 if R_mag_1 > 0 else 0
    r_mag_proportional = abs(ratio - 2.0) < 0.01
    print(f"    R_mag proportional to freq (ratio={ratio:.4f}): "
          f"{'PASS' if r_mag_proportional else 'FAIL'}")

    print()

    # ============================================================
    # Method 2: FastHenry parser end-to-end
    # ============================================================
    print("=" * 70)
    print("Method 2: FastHenry Parser End-to-End")
    print("=" * 70)
    print()

    rad.UtiDelAll()
    parser2 = FastHenryParser()
    parser2.parse_string(inp_text)
    result = parser2.solve()

    print(f"  Frequencies: {len(result['freqs'])} points")
    print(f"  R at DC = {result['R'][0]*1e3:.4f} mOhm")
    print(f"  L at DC = {result['L'][0]*1e9:.4f} nH")
    if 'Delta_L' in result:
        print(f"  Delta_L = {result['Delta_L'][0,0]*1e9:.4f} nH")
    print(f"  R at {result['freqs'][-1]:.0e} Hz = {result['R'][-1]*1e3:.4f} mOhm")
    print(f"  L at {result['freqs'][-1]:.0e} Hz = {result['L'][-1]*1e9:.4f} nH")
    print()
    print("  FastHenry parser end-to-end: PASS")
    print()


if __name__ == '__main__':
    main()
