"""
demo_fasthenry_magnetic_core.py

Demonstration of FastHenry .inp + Linear Magnetic Material coupling.

Shows how to use FastHenry format to define conductors and .magnetic blocks
to define ferrite/iron cores, then solve with CoupledPEECSolver.

Three scenarios:
  1. Wire in air (baseline)
  2. Wire + iron core (mu_r=1000, lossless)
  3. Wire + ferrite core (mu_r=2000, lossy mu"_r=200)

Key results:
  - Inductance enhancement from core material (Delta_L)
  - Magnetic loss from complex permeability (R_mag ~ omega*tan_delta*Delta_L)
  - Frequency sweep showing L(f) and R(f) behavior

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad
from fasthenry_parser import FastHenryParser

MU_0 = 4e-7 * np.pi


def scenario_1_air():
    """Scenario 1: Single wire in air (baseline)."""
    print("=" * 65)
    print("Scenario 1: Wire in air (baseline)")
    print("=" * 65)

    rad.UtiDelAll()

    inp = """\
.Units mm
.default sigma=5.8e7

* Straight copper wire: 100 mm long, 1x1 mm cross-section
N1 x=0 y=0 z=0
N2 x=100 y=0 z=0

E1 N1 N2 w=1 h=1

.external N1 N2
.freq fmin=1e3 fmax=10e6 ndec=10
.end
"""
    parser = FastHenryParser()
    parser.parse_string(inp)
    result = parser.solve()

    freqs = result['freqs']
    R = result['R']
    omega = 2.0 * np.pi * freqs
    L = result['L']

    print(f"  Frequency points: {len(freqs)}")
    print(f"  DC resistance: {R[0]*1e3:.4f} mOhm")
    print(f"  L at 1 kHz:  {L[0]*1e9:.2f} nH")
    print(f"  L at 10 MHz: {L[-1]*1e9:.2f} nH")
    print()

    return result


def scenario_2_iron_core():
    """Scenario 2: Wire + iron core (lossless, mu_r=1000)."""
    print("=" * 65)
    print("Scenario 2: Wire + iron core (mu_r=1000, lossless)")
    print("=" * 65)

    rad.UtiDelAll()

    inp = """\
.Units mm
.default sigma=5.8e7

* Straight copper wire: 100 mm long, 1x1 mm cross-section
N1 x=0 y=0 z=0
N2 x=100 y=0 z=0

E1 N1 N2 w=1 h=1

.external N1 N2

* Iron core block: 60x10x10 mm, centered 10 mm above wire
* Divisions: 3x1x1 for sufficient MSC resolution
.magnetic
  type=box
  center=50,10,0
  size=60,10,10
  divisions=3,1,1
  mu_r=1000
.endmagnetic

.freq fmin=1e3 fmax=10e6 ndec=10
.end
"""
    parser = FastHenryParser()
    parser.parse_string(inp)
    result = parser.solve()

    freqs = result['freqs']
    R = result['R']
    L = result['L']

    print(f"  DC resistance: {R[0]*1e3:.4f} mOhm")
    print(f"  L at 1 kHz:  {L[0]*1e9:.2f} nH")
    print(f"  L at 10 MHz: {L[-1]*1e9:.2f} nH")

    if 'Delta_L' in result:
        dL = result['Delta_L'][0, 0]
        print(f"  Delta_L (from core): {dL*1e9:.2f} nH")
    print()

    return result


def scenario_3_ferrite_core():
    """Scenario 3: Wire + ferrite core (lossy, mu_r=2000, mu"_r=200)."""
    print("=" * 65)
    print("Scenario 3: Wire + ferrite core (mu_r=2000, mu\"_r=200)")
    print("=" * 65)

    rad.UtiDelAll()

    inp = """\
.Units mm
.default sigma=5.8e7

* Straight copper wire: 100 mm long, 1x1 mm cross-section
N1 x=0 y=0 z=0
N2 x=100 y=0 z=0

E1 N1 N2 w=1 h=1

.external N1 N2

* Ferrite core with complex permeability:
*   mu = mu'_r - j*mu"_r = 2000 - j*200
*   Loss tangent tan(delta_m) = 200/2000 = 0.1
.magnetic
  type=box
  center=50,10,0
  size=60,10,10
  divisions=3,1,1
  mu_r=2000
  mu_r_imag=200
.endmagnetic

.freq fmin=1e3 fmax=10e6 ndec=10
.end
"""
    parser = FastHenryParser()
    parser.parse_string(inp)
    result = parser.solve()

    freqs = result['freqs']
    R = result['R']
    L = result['L']

    print(f"  DC resistance: {R[0]*1e3:.4f} mOhm")
    print(f"  L at 1 kHz:  {L[0]*1e9:.2f} nH")
    print(f"  L at 10 MHz: {L[-1]*1e9:.2f} nH")
    print(f"  R at 1 kHz:  {R[0]*1e3:.4f} mOhm")
    print(f"  R at 10 MHz: {R[-1]*1e3:.4f} mOhm")

    if 'Delta_L' in result:
        dL = result['Delta_L'][0, 0]
        tan_d = 200.0 / 2000.0
        print(f"  Delta_L (from core): {dL*1e9:.2f} nH")
        print(f"  tan(delta_m) = {tan_d:.2f}")
        print(f"  R_mag at 1 MHz = omega*tan_delta*Delta_L = "
              f"{2*np.pi*1e6*tan_d*dL*1e3:.4f} mOhm")
    print()

    return result


def scenario_4_multi_segment_coil():
    """Scenario 4: 4-segment rectangular coil around ferrite core."""
    print("=" * 65)
    print("Scenario 4: Rectangular coil around ferrite core")
    print("=" * 65)

    rad.UtiDelAll()

    # Rectangular coil: 4 segments forming a near-loop around the core
    # Port at a gap between N5 (=N1 position) and N1:
    #
    #   N1 ----E1---- N2
    #   |               |
    #   E4   [core]    E2
    #   |               |
    #   N5 ----E3---- N3
    #
    # N5 is at (-20, 0, -10), N4 at (-20, 0, 10) = N1
    # This creates a coil with port at the bottom-left corner
    #
    # Core at center (0, 0, 0), coil in x-z plane

    inp = """\
.Units mm
.default sigma=5.8e7 w=1 h=1

* Rectangular coil (40x20 mm) in x-z plane
* Gap at bottom-left corner: port = N1 to N5
N1 x=-20 y=0 z=-10
N2 x=-20 y=0 z=10
N3 x=20  y=0 z=10
N4 x=20  y=0 z=-10
N5 x=-20 y=0 z=-10.001

* Left segment (z-direction, upward)
E1 N1 N2 w=1 h=1
* Top segment (x-direction, right)
E2 N2 N3 w=1 h=1
* Right segment (z-direction, downward)
E3 N3 N4 w=1 h=1
* Bottom segment (x-direction, left)
E4 N4 N5 w=1 h=1

.external N1 N5

* Ferrite core inside the coil
* 30x8x8 mm block at center
.magnetic
  type=box
  center=0,0,0
  size=30,8,8
  divisions=3,1,1
  mu_r=2000
  mu_r_imag=100
.endmagnetic

.freq fmin=1e3 fmax=10e6 ndec=10
.end
"""
    parser = FastHenryParser()
    parser.parse_string(inp)
    result = parser.solve()

    freqs = result['freqs']
    R = result['R']
    L = result['L']

    print(f"  Coil: 4 segments, 40x20 mm rectangular")
    print(f"  Core: 30x8x8 mm ferrite (mu_r=2000, mu\"_r=100)")
    print(f"  DC resistance: {R[0]*1e3:.4f} mOhm")
    print(f"  L at 1 kHz:  {L[0]*1e9:.2f} nH")
    print(f"  L at 1 MHz:  {L[len(freqs)//2]*1e9:.2f} nH")
    print(f"  L at 10 MHz: {L[-1]*1e9:.2f} nH")
    print(f"  R at 1 kHz:  {R[0]*1e3:.4f} mOhm")
    print(f"  R at 1 MHz:  {R[len(freqs)//2]*1e3:.4f} mOhm")
    print(f"  R at 10 MHz: {R[-1]*1e3:.4f} mOhm")
    print()

    return result


def comparison_table(r_air, r_iron, r_ferrite):
    """Print comparison table across scenarios."""
    print("=" * 65)
    print("Comparison: Air vs Iron vs Ferrite Core")
    print("=" * 65)

    freqs = r_air['freqs']

    # Select representative frequencies
    idx_low = 0   # lowest freq
    idx_mid = len(freqs) // 2
    idx_high = -1  # highest freq

    freq_labels = [
        (idx_low, f"{freqs[idx_low]/1e3:.0f} kHz"),
        (idx_mid, f"{freqs[idx_mid]/1e6:.1f} MHz"),
        (idx_high, f"{freqs[idx_high]/1e6:.0f} MHz"),
    ]

    # Inductance comparison
    print("\n  Inductance L [nH]:")
    print(f"  {'Freq':>10s}  {'Air':>10s}  {'Iron':>10s}  {'Ferrite':>10s}  "
          f"{'Iron/Air':>10s}  {'Ferr/Air':>10s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")

    for idx, label in freq_labels:
        L_a = r_air['L'][idx] * 1e9
        L_i = r_iron['L'][idx] * 1e9
        L_f = r_ferrite['L'][idx] * 1e9
        ratio_i = L_i / L_a if L_a > 0 else 0
        ratio_f = L_f / L_a if L_a > 0 else 0
        print(f"  {label:>10s}  {L_a:10.2f}  {L_i:10.2f}  {L_f:10.2f}  "
              f"{ratio_i:10.2f}x  {ratio_f:10.2f}x")

    # Resistance comparison
    print(f"\n  Resistance R [mOhm]:")
    print(f"  {'Freq':>10s}  {'Air':>10s}  {'Iron':>10s}  {'Ferrite':>10s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")

    for idx, label in freq_labels:
        R_a = r_air['R'][idx] * 1e3
        R_i = r_iron['R'][idx] * 1e3
        R_f = r_ferrite['R'][idx] * 1e3
        print(f"  {label:>10s}  {R_a:10.4f}  {R_i:10.4f}  {R_f:10.4f}")

    # Delta_L summary
    if 'Delta_L' in r_iron:
        dL_iron = r_iron['Delta_L'][0, 0]
        L_air = r_air['L'][0]
        print(f"\n  Coupling (Delta_L):")
        print(f"    Iron core:    Delta_L = {dL_iron*1e9:.2f} nH "
              f"({dL_iron/L_air*100:.1f}% of L_air)")

    if 'Delta_L' in r_ferrite:
        dL_ferr = r_ferrite['Delta_L'][0, 0]
        L_air = r_air['L'][0]
        print(f"    Ferrite core: Delta_L = {dL_ferr*1e9:.2f} nH "
              f"({dL_ferr/L_air*100:.1f}% of L_air)")
        tan_d = 200.0 / 2000.0
        print(f"    Ferrite loss: tan(delta_m) = {tan_d:.2f}")
        print(f"    R_mag(1MHz) = omega * tan_delta * Delta_L = "
              f"{2*np.pi*1e6*tan_d*dL_ferr*1e3:.4f} mOhm")

    print()


def print_key_takeaways():
    """Print key physics and method selection guide."""
    print("=" * 65)
    print("KEY TAKEAWAYS")
    print("=" * 65)
    print("""
  FastHenry .inp + Radia MMM Coupling:
  ------------------------------------
  1. PEEC defines conductor geometry (nodes, segments, ports)
  2. .magnetic block defines iron/ferrite core (box, hex, or GMSH mesh)
  3. CoupledPEECSolver computes Delta_L via:
     Unit current -> Biot-Savart H -> Radia Solve -> A-field -> Delta_L
  4. Delta_L is frequency-independent for linear materials (computed once)

  Complex Permeability (mu = mu'_r - j*mu"_r):
  ---------------------------------------------
  - mu'_r: Energy storage (real inductance Delta_L)
  - mu"_r: Magnetic loss (core heating)
  - Loss resistance: R_mag(f) = omega * (mu"_r / mu'_r) * Delta_L
  - Total Z(f) = R_dc + R_mag(f) + j*omega*(L_air + Delta_L)

  Method Selection:
  -----------------
  | Application              | Method                     |
  |--------------------------|----------------------------|
  | Inductor with core       | FastHenry + .magnetic box  |
  | Transformer coupling     | FastHenry + multiple coils |
  | Ferrite loss estimation  | .magnetic mu_r_imag > 0    |
  | Complex core geometry    | .magnetic type=mesh (.msh) |
  | Skin effect in conductor | nwinc/nhinc on segments    |
""")


if __name__ == '__main__':
    print()
    print("*" * 65)
    print("  FastHenry + Linear Magnetic Material Demo")
    print("  PEEC conductor + Radia MMM core coupling")
    print("*" * 65)
    print()

    # Run three comparable scenarios
    r_air = scenario_1_air()
    r_iron = scenario_2_iron_core()
    r_ferrite = scenario_3_ferrite_core()

    # Comparison table
    comparison_table(r_air, r_iron, r_ferrite)

    # Rectangular coil demo
    r_coil = scenario_4_multi_segment_coil()

    # Summary
    print_key_takeaways()
