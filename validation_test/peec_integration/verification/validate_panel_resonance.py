"""
Validation: PEEC Panel + Resonance (Phase A-D Integration Test)

Tests the full RLCM PEEC system with capacitive panels:
1. Loop-only vs RLCM consistency at low frequency
2. Panel effects at high frequency (capacitive regime)
3. Resonance detection (Im(Z) = 0 crossing)
4. FastHenry .panel block end-to-end
5. Surface mesh -> PEEC with panels

Author: Radia Development Team
Date: 2026-02-16
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

from peec_matrices import PEECBuilder
from peec_topology import PEECCircuitSolver
from fasthenry_parser import FastHenryParser
from peec_mesh_import import create_plate_mesh

MU_0 = 4 * np.pi * 1e-7

print("=" * 70)
print("PEEC Panel + Resonance Validation")
print("=" * 70)

n_pass = 0
n_fail = 0
n_skip = 0


def check(condition, name, detail=""):
    global n_pass, n_fail
    if condition:
        n_pass += 1
        print("  PASS: %s" % name)
    else:
        n_fail += 1
        print("  FAIL: %s %s" % (name, detail))


# =========================================================================
# Test 1: Loop-only backward compatibility
# =========================================================================
print("\n--- Test 1: Loop-only backward compatibility ---")

builder = PEECBuilder()
n1 = builder.add_node_at(0, 0, 0)
n2 = builder.add_node_at(0.1, 0, 0)
builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
builder.add_port(n1, n2)

topo_nostar = builder.build_topology(include_star=False)
solver_nostar = PEECCircuitSolver(topo_nostar)
Z_nostar = solver_nostar.compute_port_impedance(1e6)

check(not solver_nostar.has_panels, "No panels when include_star=False")
check(Z_nostar.imag > 0, "Inductive at 1 MHz",
      "X=%.4e" % Z_nostar.imag)

# =========================================================================
# Test 2: With panels, RLCM mode activated
# =========================================================================
print("\n--- Test 2: RLCM mode with panels ---")

builder2 = PEECBuilder()
n1 = builder2.add_node_at(0, 0, 0)
n2 = builder2.add_node_at(0.1, 0, 0)
builder2.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
builder2.add_port(n1, n2)

# Add panel plates above the wire
for i in range(4):
    x0 = i * 0.025
    x1 = (i + 1) * 0.025
    builder2.add_panel([[x0, -0.01, 0.003],
                         [x1, -0.01, 0.003],
                         [x1, 0.01, 0.003]])
    builder2.add_panel([[x0, -0.01, 0.003],
                         [x1, 0.01, 0.003],
                         [x0, 0.01, 0.003]])

topo2 = builder2.build_topology(include_star=True)
solver2 = PEECCircuitSolver(topo2)

check(solver2.has_panels, "has_panels=True with panels")
check(topo2['n_star'] == 8, "n_star=8 (8 triangles)",
      "got %d" % topo2['n_star'])
check(topo2['P'].shape == (8, 8), "P matrix 8x8",
      "got %s" % str(topo2['P'].shape))
check(topo2['M_LS'].shape[1] == 8, "M_LS n_cols=8",
      "got %d" % topo2['M_LS'].shape[1])

# =========================================================================
# Test 3: Low-frequency consistency
# =========================================================================
print("\n--- Test 3: Low-frequency consistency ---")

# At low frequency, panels should have minimal effect
# (capacitors are open circuit at low freq)
Z_rlcm_low = solver2.compute_port_impedance(100.0)  # 100 Hz
Z_loop_low = solver_nostar.compute_port_impedance(100.0)

# R should be similar (same wire)
R_ratio = abs(Z_rlcm_low.real - Z_loop_low.real) / abs(Z_loop_low.real)
check(R_ratio < 0.01, "R consistent at 100 Hz (%.2f%% diff)" % (R_ratio * 100))

# X should be inductive (positive) at low freq
check(Z_rlcm_low.imag > 0, "Inductive at 100 Hz",
      "X=%.4e" % Z_rlcm_low.imag)

# =========================================================================
# Test 4: DC fallback
# =========================================================================
print("\n--- Test 4: DC fallback ---")

Z_dc = solver2.compute_port_impedance(0.0)
check(abs(Z_dc.imag) < 1e-10, "DC has zero reactance",
      "X=%.4e" % Z_dc.imag)
check(Z_dc.real > 0, "DC has positive resistance",
      "R=%.4e" % Z_dc.real)

# =========================================================================
# Test 5: Frequency sweep physical behavior
# =========================================================================
print("\n--- Test 5: Frequency sweep physical behavior ---")

freqs = np.logspace(3, 10, 100)
Z_sweep = solver2.frequency_sweep(freqs)

# At low frequency: inductive (X > 0)
check(Z_sweep[0].imag > 0, "Inductive at f_min",
      "X=%.4e at f=%.0e" % (Z_sweep[0].imag, freqs[0]))

# Impedance magnitude increases with frequency (at least at low freq)
Z_mag_low = abs(Z_sweep[10])
Z_mag_mid = abs(Z_sweep[50])
check(Z_mag_mid > Z_mag_low, "Z increases with frequency",
      "|Z_mid|=%.4e > |Z_low|=%.4e" % (Z_mag_mid, Z_mag_low))

# =========================================================================
# Test 6: Resonance detection
# =========================================================================
print("\n--- Test 6: Resonance detection ---")

res = solver2.get_resonant_frequencies(1e3, 1e12, 10000)
if len(res) > 0:
    for i, fr in enumerate(res[:3]):
        Z_res = solver2.compute_port_impedance(fr)
        check(abs(Z_res.imag) < abs(Z_res.real) * 10,
              "Near resonance: Im(Z)/Re(Z) < 10 at f=%.3e" % fr,
              "Im=%.4e, Re=%.4e" % (Z_res.imag, Z_res.real))
    print("  INFO: Found %d resonance(s)" % len(res))
else:
    print("  INFO: No resonance found in scan range (acceptable)")
    n_skip += 1

# =========================================================================
# Test 7: FastHenry .panel block end-to-end
# =========================================================================
print("\n--- Test 7: FastHenry .panel block ---")

parser = FastHenryParser()
parser.parse_string("""
.Units mm
.default sigma=5.8e7

N1 x=0 y=0 z=0
N2 x=100 y=0 z=0

E1 N1 N2 w=1 h=1 nwinc=3 nhinc=3

.external N1 N2

.panel
  type=plate
  center=50,0,5
  size=80,20
  divisions=4,2
.endpanel

.freq fmin=1e3 fmax=1e9 ndec=5
.end
""")

summary = parser.get_summary()
check(summary['n_panel_blocks'] == 1, "1 panel block parsed")
check(summary['n_nodes'] == 2, "2 nodes parsed")
check(summary['n_segments'] == 1, "1 segment parsed")

result = parser.solve()
check(len(result['freqs']) > 0, "Frequency sweep computed")
check(result['Z_port'].shape == result['freqs'].shape, "Z_port shape matches freqs")

# Check that inductance is reasonable
L_low = result['L'][0]
check(L_low > 0, "Positive inductance",
      "L=%.4e H" % L_low)

# =========================================================================
# Test 8: Surface mesh -> PEEC with panels
# =========================================================================
print("\n--- Test 8: Surface mesh -> PEEC with panels ---")

builder8 = create_plate_mesh(
    center=[0.0, 0.0, 0.005],
    size=[0.08, 0.02],
    divisions=[4, 2],
    sigma=5.8e7,
    thickness=1e-3,
    ports=[(0, 14)],
    include_panels=True
)

topo8 = builder8.build_topology(include_star=True)
solver8 = PEECCircuitSolver(topo8)

check(solver8.has_panels, "Mesh-based panels detected")
check(topo8['n_star'] > 0, "Star DOFs from mesh panels",
      "n_star=%d" % topo8['n_star'])

Z8 = solver8.compute_port_impedance(1e6)
check(Z8.imag != 0, "Non-zero reactance at 1 MHz",
      "X=%.4e" % Z8.imag)

# =========================================================================
# Test 9: Multiple .panel blocks
# =========================================================================
print("\n--- Test 9: Multiple .panel blocks ---")

parser9 = FastHenryParser()
parser9.parse_string("""
.Units m
N1 x=0 y=0 z=0
N2 x=0.1 y=0 z=0
E1 N1 N2 w=1e-3 h=1e-3 sigma=5.8e7
.external N1 N2

.panel
  type=triangle
  vertices=0,0,0.005, 0.05,0,0.005, 0.025,0.02,0.005
.endpanel

.panel
  type=triangle
  vertices=0.05,0,0.005, 0.1,0,0.005, 0.075,0.02,0.005
.endpanel

.freq fmin=1e6 fmax=1e6 ndec=1
.end
""")

check(parser9.get_summary()['n_panel_blocks'] == 2, "2 panel blocks parsed")
result9 = parser9.solve()
check(len(result9['Z_port']) == 1, "Single frequency solved")

# =========================================================================
# Summary
# =========================================================================
print("\n" + "=" * 70)
total = n_pass + n_fail
print("Results: %d/%d PASSED, %d FAILED, %d SKIPPED" %
      (n_pass, total, n_fail, n_skip))
if n_fail == 0:
    print("ALL TESTS PASSED")
else:
    print("SOME TESTS FAILED")
print("=" * 70)

sys.exit(n_fail)
