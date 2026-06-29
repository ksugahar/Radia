"""
demo_circular_coil_4cases.py

4-case comparison of circular coil (R=20mm) with FastHenry-based PEEC:

  Case 1: Air core (baseline)
  Case 2: + iron core (mu_r=1000, lossless)   -- .magnetic
  Case 3: + Al shield                          -- .shield
  Case 4: + iron core + Al shield              -- .magnetic + .shield

Includes Bessel SIBC for wire skin effect (circular cross-section).

Demonstrates:
  - PEEC inductance and resistance for circular coil
  - Bessel SIBC: AC resistance increase from skin effect
  - Magnetic coupling (Delta_L > 0, L increases)
  - Shield coupling (Delta_Z, L decreases by Lenz's law)
  - Combined magnetic + shield coupling

Part of Radia project
"""

import sys
import os
import time
import numpy as np
from scipy.special import iv

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad
from fasthenry_parser import FastHenryParser

MU_0 = 4e-7 * np.pi


def bessel_impedance_per_length(freq, radius, sigma, mu_r=1.0):
    """Internal impedance per unit length of solid circular wire (exact).

    Z_i/l = k/(2*pi*a*sigma) * I0(ka)/I1(ka)

    where k = sqrt(j*omega*mu*sigma), I0/I1 are modified Bessel functions.

    DC limit: Z_i/l -> 1/(sigma*pi*a^2) + j*omega*mu/(8*pi)
    """
    if freq <= 0:
        return 1.0 / (sigma * np.pi * radius**2) + 0j

    omega = 2.0 * np.pi * freq
    mu = mu_r * MU_0
    k = np.sqrt(1j * omega * mu * sigma)
    ka = k * radius

    I0_ka = iv(0, ka)
    I1_ka = iv(1, ka)

    return k / (2.0 * np.pi * radius * sigma) * (I0_ka / I1_ka)


def make_bessel_Zs_func(n_filaments, seg_lengths, wire_w, wire_h, sigma):
    """Create Zs_func for Bessel SIBC on rectangular wire with equivalent radius.

    Zs_func(freq) returns (n_filaments,) complex array.
    Zs = Z_internal_per_seg - R_dc_per_seg  (because R_dc is already in PEEC).

    Args:
        n_filaments: Total number of filaments (after nwinc*nhinc expansion)
        seg_lengths: Array of segment lengths [m]
        wire_w, wire_h: Wire cross-section dimensions [m]
        sigma: Conductivity [S/m]

    Returns:
        Callable(freq) -> Zs array (n_filaments,) complex
    """
    # Equivalent circular radius for rectangular cross-section
    r_equiv = np.sqrt(wire_w * wire_h / np.pi)
    cross_area = wire_w * wire_h
    R_dc_per_length = 1.0 / (sigma * cross_area)

    def Zs_func(freq):
        if freq <= 0:
            return np.zeros(n_filaments, dtype=complex)

        Z_per_l = bessel_impedance_per_length(freq, r_equiv, sigma)
        # Zs = (Z_internal/l - R_dc/l) * length per segment
        Zs_per_l = Z_per_l - R_dc_per_length
        Zs = np.array([Zs_per_l * seg_lengths[i] for i in range(n_filaments)])
        return Zs

    return Zs_func


def generate_circular_coil_inp(R_mm=20.0, n_seg=16, w=1.0, h=1.0,
                                sigma=5.8e7, nwinc=1, nhinc=1):
    """Generate FastHenry .inp text for a circular coil.

    Coil lies in the x-y plane at z=0.

    Args:
        R_mm: Coil radius in mm
        n_seg: Number of straight segments
        w, h: Wire cross-section width and height in mm
        sigma: Conductivity in S/m
        nwinc, nhinc: Multi-filament subdivision
    Returns:
        str: FastHenry .inp node+segment+port definitions (no .freq, .end)
    """
    lines = []
    lines.append(f"* Circular coil: R={R_mm} mm, {n_seg} segments, "
                 f"{w}x{h} mm Cu wire")

    # Generate node positions around the circle
    angles = np.linspace(0, 2 * np.pi, n_seg + 1)[:-1]
    for i, theta in enumerate(angles):
        x = R_mm * np.cos(theta)
        y = R_mm * np.sin(theta)
        z = 0.0
        lines.append(f"N{i+1} x={x:.6f} y={y:.6f} z={z:.6f}")

    # Gap node for port (tiny z-offset from N1)
    x0 = R_mm * np.cos(0)
    y0 = R_mm * np.sin(0)
    gap_id = n_seg + 1
    lines.append(f"N{gap_id} x={x0:.6f} y={y0:.6f} z=-0.001")

    lines.append("")

    # Generate segments connecting consecutive nodes
    nw_str = f" nwinc={nwinc}" if nwinc > 1 else ""
    nh_str = f" nhinc={nhinc}" if nhinc > 1 else ""
    for i in range(n_seg - 1):
        lines.append(f"E{i+1} N{i+1} N{i+2} w={w} h={h} "
                     f"sigma={sigma:.2e}{nw_str}{nh_str}")
    # Last segment: N_{n_seg} -> N_{gap_id}
    lines.append(f"E{n_seg} N{n_seg} N{gap_id} w={w} h={h} "
                 f"sigma={sigma:.2e}{nw_str}{nh_str}")

    lines.append("")
    lines.append(f".external N1 N{gap_id}")

    return "\n".join(lines)


def case_1_air(coil_inp, freqs_str, wire_w_m, wire_h_m, sigma):
    """Case 1: Circular coil in air (baseline)."""
    print("=" * 65)
    print("Case 1: Circular coil in air (baseline)")
    print("=" * 65)

    rad.UtiDelAll()

    inp = f"""\
.Units mm
.default sigma=5.8e7

{coil_inp}

{freqs_str}
.end
"""
    parser = FastHenryParser()
    parser.parse_string(inp)

    # Build topology first to get segment info for Bessel SIBC
    builder = parser.to_peec_builder()
    topo = builder.build_topology()
    seg_lengths = np.asarray(topo['segment_lengths'])
    n_fil = topo['n_loop']

    Zs_func = make_bessel_Zs_func(n_fil, seg_lengths, wire_w_m, wire_h_m, sigma)

    t0 = time.perf_counter()
    result = parser.solve(Zs_func=Zs_func)
    t_solve = time.perf_counter() - t0

    freqs = result['freqs']
    R = result['R']
    L = result['L']
    result['t_solve'] = t_solve

    print(f"  Frequency points: {len(freqs)}")
    print(f"  DC resistance:  {R[0]*1e3:.4f} mOhm")
    print(f"  R at {freqs[-1]/1e6:.0f} MHz:  {R[-1]*1e3:.4f} mOhm")
    print(f"  L at {freqs[0]/1e3:.0f} kHz:   {L[0]*1e9:.2f} nH")
    print(f"  L at {freqs[-1]/1e6:.0f} MHz:  {L[-1]*1e9:.2f} nH")
    print(f"  Solve time:     {t_solve:.3f} s")

    # Skin effect info
    delta = np.sqrt(2.0 / (2 * np.pi * freqs[-1] * MU_0 * sigma))
    r_equiv = np.sqrt(wire_w_m * wire_h_m / np.pi)
    print(f"  Skin depth at {freqs[-1]/1e6:.0f} MHz: {delta*1e3:.4f} mm "
          f"(wire r_eq = {r_equiv*1e3:.3f} mm)")
    print()

    return result


def case_2_magnetic(coil_inp, freqs_str, wire_w_m, wire_h_m, sigma):
    """Case 2: Circular coil + iron core (lossless, mu_r=1000)."""
    print("=" * 65)
    print("Case 2: Circular coil + iron core (mu_r=1000)")
    print("=" * 65)

    rad.UtiDelAll()

    inp = f"""\
.Units mm
.default sigma=5.8e7

{coil_inp}

* Iron core: 15x15x10 mm box at center of coil
* Divisions: 2x2x1 for MSC resolution
.magnetic
  type=box
  center=0,0,0
  size=15,15,10
  divisions=2,2,1
  mu_r=1000
.endmagnetic

{freqs_str}
.end
"""
    parser = FastHenryParser()
    parser.parse_string(inp)

    builder = parser.to_peec_builder()
    topo = builder.build_topology()
    seg_lengths = np.asarray(topo['segment_lengths'])
    n_fil = topo['n_loop']

    Zs_func = make_bessel_Zs_func(n_fil, seg_lengths, wire_w_m, wire_h_m, sigma)

    t0 = time.perf_counter()
    result = parser.solve(Zs_func=Zs_func)
    t_solve = time.perf_counter() - t0

    freqs = result['freqs']
    R = result['R']
    L = result['L']
    result['t_solve'] = t_solve

    print(f"  DC resistance:  {R[0]*1e3:.4f} mOhm")
    print(f"  R at {freqs[-1]/1e6:.0f} MHz:  {R[-1]*1e3:.4f} mOhm")
    print(f"  L at {freqs[0]/1e3:.0f} kHz:   {L[0]*1e9:.2f} nH")
    print(f"  L at {freqs[-1]/1e6:.0f} MHz:  {L[-1]*1e9:.2f} nH")
    print(f"  Solve time:     {t_solve:.3f} s")

    if 'Delta_L' in result:
        dL = np.diag(result['Delta_L']).sum()
        print(f"  Delta_L (from core): {dL*1e9:.2f} nH")
    print()

    return result


def case_3_shield(coil_inp, freqs_str, wire_w_m, wire_h_m, sigma):
    """Case 3: Circular coil + Al shield."""
    print("=" * 65)
    print("Case 3: Circular coil + Al shield")
    print("=" * 65)

    rad.UtiDelAll()

    inp = f"""\
.Units mm
.default sigma=5.8e7

{coil_inp}

* Aluminum shield box: 50x50x10 mm close to coil
* z-faces at +/-5 mm from coil plane -> strong eddy current coupling
* sigma_Al = 3.7e7 S/m
.shield
  type=box
  center=0,0,0
  size=50,50,10
  sigma=3.7e7
.endshield

{freqs_str}
.end
"""
    parser = FastHenryParser()
    parser.parse_string(inp)

    builder = parser.to_peec_builder()
    topo = builder.build_topology()
    seg_lengths = np.asarray(topo['segment_lengths'])
    n_fil = topo['n_loop']

    Zs_func = make_bessel_Zs_func(n_fil, seg_lengths, wire_w_m, wire_h_m, sigma)

    t0 = time.perf_counter()
    result = parser.solve(Zs_func=Zs_func)
    t_solve = time.perf_counter() - t0

    freqs = result['freqs']
    R = result['R']
    L = result['L']
    result['t_solve'] = t_solve

    print(f"  DC resistance:  {R[0]*1e3:.4f} mOhm")
    print(f"  R at {freqs[-1]/1e6:.0f} MHz:  {R[-1]*1e3:.4f} mOhm")
    print(f"  L at {freqs[0]/1e3:.0f} kHz:   {L[0]*1e9:.2f} nH")
    print(f"  L at {freqs[-1]/1e6:.0f} MHz:  {L[-1]*1e9:.2f} nH")
    print(f"  Solve time:     {t_solve:.3f} s")
    print()

    return result


def case_4_magnetic_shield(coil_inp, freqs_str, wire_w_m, wire_h_m, sigma):
    """Case 4: Circular coil + iron core + Al shield (combined)."""
    print("=" * 65)
    print("Case 4: Circular coil + iron core + Al shield (combined)")
    print("=" * 65)

    rad.UtiDelAll()

    inp = f"""\
.Units mm
.default sigma=5.8e7

{coil_inp}

* Iron core: 15x15x10 mm box at center of coil
.magnetic
  type=box
  center=0,0,0
  size=15,15,10
  divisions=2,2,1
  mu_r=1000
.endmagnetic

* Aluminum shield box: 50x50x10 mm close to coil
.shield
  type=box
  center=0,0,0
  size=50,50,10
  sigma=3.7e7
.endshield

{freqs_str}
.end
"""
    parser = FastHenryParser()
    parser.parse_string(inp)

    builder = parser.to_peec_builder()
    topo = builder.build_topology()
    seg_lengths = np.asarray(topo['segment_lengths'])
    n_fil = topo['n_loop']

    Zs_func = make_bessel_Zs_func(n_fil, seg_lengths, wire_w_m, wire_h_m, sigma)

    t0 = time.perf_counter()
    result = parser.solve(Zs_func=Zs_func)
    t_solve = time.perf_counter() - t0

    freqs = result['freqs']
    R = result['R']
    L = result['L']
    result['t_solve'] = t_solve

    print(f"  DC resistance:  {R[0]*1e3:.4f} mOhm")
    print(f"  R at {freqs[-1]/1e6:.0f} MHz:  {R[-1]*1e3:.4f} mOhm")
    print(f"  L at {freqs[0]/1e3:.0f} kHz:   {L[0]*1e9:.2f} nH")
    print(f"  L at {freqs[-1]/1e6:.0f} MHz:  {L[-1]*1e9:.2f} nH")
    print(f"  Solve time:     {t_solve:.3f} s")

    if 'Delta_L' in result:
        dL = np.diag(result['Delta_L']).sum()
        print(f"  Delta_L (from core): {dL*1e9:.2f} nH")
    print()

    return result


def comparison_table(r1, r2, r3, r4):
    """Print comparison table across all 4 cases."""
    print("=" * 75)
    print("Comparison Table: Circular Coil (R=20mm)")
    print("=" * 75)

    freqs = r1['freqs']

    # Select representative frequencies
    idx_low = 0
    idx_mid = len(freqs) // 2
    idx_high = -1

    freq_labels = [
        (idx_low, f"{freqs[idx_low]/1e3:.0f} kHz"),
        (idx_mid, f"{freqs[idx_mid]/1e6:.1f} MHz"),
        (idx_high, f"{freqs[idx_high]/1e6:.0f} MHz"),
    ]

    # Inductance comparison
    print("\n  Inductance L [nH]:")
    print(f"  {'Freq':>10s}  {'Air':>8s}  {'+Core':>8s}  {'+Shield':>8s}  "
          f"{'+Both':>8s}  {'Core%':>8s}  {'Shld%':>8s}  {'Both%':>8s}")
    print(f"  {'-'*10}  {'-'*8}  {'-'*8}  {'-'*8}  "
          f"{'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")

    for idx, label in freq_labels:
        La = r1['L'][idx] * 1e9
        Lm = r2['L'][idx] * 1e9
        Ls = r3['L'][idx] * 1e9
        Lb = r4['L'][idx] * 1e9
        dm = (Lm - La) / La * 100 if La > 0 else 0
        ds = (Ls - La) / La * 100 if La > 0 else 0
        db = (Lb - La) / La * 100 if La > 0 else 0
        print(f"  {label:>10s}  {La:8.2f}  {Lm:8.2f}  {Ls:8.2f}  "
              f"{Lb:8.2f}  {dm:+7.1f}%  {ds:+7.1f}%  {db:+7.1f}%")

    # Resistance comparison
    print(f"\n  Resistance R [mOhm]:")
    print(f"  {'Freq':>10s}  {'Air':>10s}  {'+Core':>10s}  {'+Shield':>10s}  "
          f"{'+Both':>10s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")

    for idx, label in freq_labels:
        Ra = r1['R'][idx] * 1e3
        Rm = r2['R'][idx] * 1e3
        Rs = r3['R'][idx] * 1e3
        Rb = r4['R'][idx] * 1e3
        print(f"  {label:>10s}  {Ra:10.4f}  {Rm:10.4f}  {Rs:10.4f}  "
              f"{Rb:10.4f}")

    # Physics summary
    print("\n  Physics Summary:")
    L_air = r1['L'][0] * 1e9
    L_core = r2['L'][0] * 1e9
    L_shield = r3['L'][0] * 1e9
    L_both = r4['L'][0] * 1e9

    print(f"    L_air    = {L_air:.2f} nH (baseline)")
    print(f"    L_core   = {L_core:.2f} nH  -> Delta_L_mag = "
          f"{(L_core - L_air):+.2f} nH ({(L_core-L_air)/L_air*100:+.1f}%)")
    print(f"    L_shield = {L_shield:.2f} nH  -> Delta_L_shd = "
          f"{(L_shield - L_air):+.2f} nH ({(L_shield-L_air)/L_air*100:+.1f}%)")
    print(f"    L_both   = {L_both:.2f} nH  -> Delta_L_tot = "
          f"{(L_both - L_air):+.2f} nH ({(L_both-L_air)/L_air*100:+.1f}%)")

    # Check additivity
    delta_core = L_core - L_air
    delta_shld = L_shield - L_air
    delta_both = L_both - L_air
    delta_sum = delta_core + delta_shld
    if abs(delta_sum) > 0:
        print(f"\n    Additivity check:")
        print(f"      Delta_core + Delta_shield = "
              f"{delta_core:+.2f} + {delta_shld:+.2f} = {delta_sum:+.2f} nH")
        print(f"      Delta_both (actual)       = {delta_both:+.2f} nH")
        print(f"      Difference                = "
              f"{(delta_both - delta_sum):+.2f} nH "
              f"(cross-coupling)")

    # Timing summary
    t_total = sum(r['t_solve'] for r in [r1, r2, r3, r4])
    print(f"\n  Timing (PEEC+ngbem):")
    print(f"    Case 1 (Air):     {r1['t_solve']:6.2f} s")
    print(f"    Case 2 (+Core):   {r2['t_solve']:6.2f} s")
    print(f"    Case 3 (+Shield): {r3['t_solve']:6.2f} s")
    print(f"    Case 4 (+Both):   {r4['t_solve']:6.2f} s")
    print(f"    Total:            {t_total:6.2f} s")
    print()
    print(f"  For comparison, NGSolve FEM (PARDISO) takes ~15 min for same problem.")
    print(f"  Speedup: ~{15*60/t_total:.0f}x")

    print()


if __name__ == '__main__':
    print()
    print("*" * 75)
    print("  Circular Coil (R=20mm) -- 4-Case PEEC Comparison")
    print("  Case 1: Air | Case 2: +Core | Case 3: +Shield | Case 4: +Both")
    print("  With Bessel SIBC for wire skin effect")
    print("*" * 75)
    print()

    # Wire parameters
    wire_w_mm = 1.0   # mm
    wire_h_mm = 1.0   # mm
    sigma_cu = 5.8e7  # S/m

    # Convert to meters for Bessel SIBC
    wire_w_m = wire_w_mm * 1e-3
    wire_h_m = wire_h_mm * 1e-3
    r_equiv = np.sqrt(wire_w_m * wire_h_m / np.pi)

    # Generate circular coil definition (shared across all cases)
    coil_inp = generate_circular_coil_inp(
        R_mm=20.0, n_seg=16, w=wire_w_mm, h=wire_h_mm, sigma=sigma_cu)
    freqs_str = ".freq fmin=1e3 fmax=10e6 ndec=10"

    print("Coil definition:")
    print(f"  Radius: 20 mm, Segments: 16, Wire: {wire_w_mm}x{wire_h_mm} mm Cu")
    print(f"  Equivalent wire radius: {r_equiv*1e3:.3f} mm")
    print(f"  Bessel SIBC: scipy.special.iv (exact circular wire impedance)")
    print()

    # Run all 4 cases
    r1 = case_1_air(coil_inp, freqs_str, wire_w_m, wire_h_m, sigma_cu)
    r2 = case_2_magnetic(coil_inp, freqs_str, wire_w_m, wire_h_m, sigma_cu)
    r3 = case_3_shield(coil_inp, freqs_str, wire_w_m, wire_h_m, sigma_cu)
    r4 = case_4_magnetic_shield(coil_inp, freqs_str, wire_w_m, wire_h_m, sigma_cu)

    # Comparison table
    comparison_table(r1, r2, r3, r4)

    rad.UtiDelAll()
