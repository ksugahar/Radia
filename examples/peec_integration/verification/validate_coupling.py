"""
Validation: PEEC Multi-Port Z-Matrix and Coupling Coefficient

Tests:
1. Two parallel wires - coupling depends on separation
2. Z-matrix reciprocity (Z12 = Z21)
3. Single-port consistency (Z11 matches compute_port_impedance)
4. k range: 0 < k < 1
5. k increases as wire separation decreases
6. FastHenry two-coil .inp parsing
7. Two coaxial square loops (WPT-like geometry)
8. Multi-filament coupling

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

MU_0 = 4 * np.pi * 1e-7

print("=" * 70)
print("PEEC Multi-Port Z-Matrix and Coupling Coefficient Validation")
print("=" * 70)

n_pass = 0
n_fail = 0


def check(condition, name, detail=""):
    global n_pass, n_fail
    if condition:
        n_pass += 1
        print("  PASS: %s" % name)
    else:
        n_fail += 1
        print("  FAIL: %s %s" % (name, detail))


# =========================================================================
# Test 1: Two parallel wires - basic coupling
# =========================================================================
print("\n--- Test 1: Two parallel wires - basic coupling ---")

# Wire 1: (0,0,0) -> (0.1,0,0)
# Wire 2: (0,d,0) -> (0.1,d,0)  with d = 5mm
d_sep = 0.005  # 5mm separation
wire_len = 0.1  # 100mm
w = 1e-3  # 1mm cross-section
sigma = 5.8e7

builder = PEECBuilder()
# Wire 1 nodes
n1 = builder.add_node_at(0, 0, 0)
n2 = builder.add_node_at(wire_len, 0, 0)
# Wire 2 nodes
n3 = builder.add_node_at(0, d_sep, 0)
n4 = builder.add_node_at(wire_len, d_sep, 0)

builder.add_connected_segment(n1, n2, w, w, sigma)
builder.add_connected_segment(n3, n4, w, w, sigma)

# Two ports: one per wire
builder.add_port(n1, n2)  # Port 0: wire 1
builder.add_port(n3, n4)  # Port 1: wire 2

topo = builder.build_topology()
solver = PEECCircuitSolver(topo)

result = solver.compute_coupling_coefficient(1e6)

check(len(solver.ports) == 2, "Two ports defined")
check(result['Z_matrix'].shape == (2, 2), "Z-matrix is 2x2",
      "got %s" % str(result['Z_matrix'].shape))
check(result['k'] > 0, "k > 0 (positive coupling)",
      "k=%.4f" % result['k'])
check(result['k'] < 1, "k < 1 (physically valid)",
      "k=%.4f" % result['k'])
check(result['L1'] > 0, "L1 > 0",
      "L1=%.4e H" % result['L1'])
check(result['L2'] > 0, "L2 > 0",
      "L2=%.4e H" % result['L2'])
check(result['M'] > 0, "M > 0 (positive mutual inductance)",
      "M=%.4e H" % result['M'])

print("  INFO: k=%.4f, L1=%.2f nH, L2=%.2f nH, M=%.2f nH" %
      (result['k'], result['L1'] * 1e9, result['L2'] * 1e9, result['M'] * 1e9))

# =========================================================================
# Test 2: Z-matrix reciprocity (Z12 = Z21)
# =========================================================================
print("\n--- Test 2: Z-matrix reciprocity ---")

Z = result['Z_matrix']
Z12 = Z[0, 1]
Z21 = Z[1, 0]
recip_err = abs(Z12 - Z21) / abs(Z12) if abs(Z12) > 1e-30 else 0.0

check(recip_err < 1e-6, "Z12 = Z21 (reciprocity, %.2e%% error)" % (recip_err * 100),
      "Z12=%.4e, Z21=%.4e" % (abs(Z12), abs(Z21)))

# =========================================================================
# Test 3: Single-port consistency (Z11 matches compute_port_impedance)
# =========================================================================
print("\n--- Test 3: Single-port consistency ---")

# Z11 from Z-matrix should match single-port impedance
Z11_matrix = Z[0, 0]
Z11_single = solver.compute_port_impedance(1e6)

# These should be very close (both use the same MNA formulation)
z_diff = abs(Z11_matrix - Z11_single) / abs(Z11_single) if abs(Z11_single) > 1e-30 else 0.0
check(z_diff < 1e-6, "Z11(matrix) = Z_port(single) (%.2e%% diff)" % (z_diff * 100),
      "|Z11|=%.4e, |Z_single|=%.4e" % (abs(Z11_matrix), abs(Z11_single)))

# =========================================================================
# Test 4: L1 ~= L2 for identical wires
# =========================================================================
print("\n--- Test 4: Symmetric self-inductance ---")

L_ratio = abs(result['L1'] - result['L2']) / result['L1']
check(L_ratio < 0.01, "L1 ~= L2 for identical wires (%.2f%% diff)" % (L_ratio * 100))

# =========================================================================
# Test 5: k increases as separation decreases
# =========================================================================
print("\n--- Test 5: k vs separation distance ---")

k_values = []
separations = [0.002, 0.005, 0.01, 0.02, 0.05]

for d in separations:
    b = PEECBuilder()
    na = b.add_node_at(0, 0, 0)
    nb = b.add_node_at(wire_len, 0, 0)
    nc = b.add_node_at(0, d, 0)
    nd = b.add_node_at(wire_len, d, 0)
    b.add_connected_segment(na, nb, w, w, sigma)
    b.add_connected_segment(nc, nd, w, w, sigma)
    b.add_port(na, nb)
    b.add_port(nc, nd)
    t = b.build_topology()
    s = PEECCircuitSolver(t)
    r = s.compute_coupling_coefficient(1e6)
    k_values.append(r['k'])
    print("  INFO: d=%.1f mm -> k=%.4f" % (d * 1e3, r['k']))

# k should decrease monotonically with distance
monotonic = all(k_values[i] > k_values[i + 1] for i in range(len(k_values) - 1))
check(monotonic, "k decreases monotonically with separation")

# Closest wires should have highest coupling
check(k_values[0] > 0.5, "Close wires (2mm) have k > 0.5",
      "k=%.4f" % k_values[0])

# Farthest wires should have lowest coupling
check(k_values[-1] < k_values[0], "Far wires have lower k")

# =========================================================================
# Test 6: FastHenry two-coil .inp
# =========================================================================
print("\n--- Test 6: FastHenry two-coil model ---")

parser = FastHenryParser()
parser.parse_string("""
.Units mm
.default sigma=5.8e7

* Coil 1: U-shaped open loop at y=0 (3 segments)
N1 x=0 y=0 z=0
N2 x=50 y=0 z=0
N3 x=50 y=0 z=50
N4 x=0 y=0 z=50

E1 N1 N2 w=1 h=1
E2 N2 N3 w=1 h=1
E3 N3 N4 w=1 h=1

* Coil 2: U-shaped open loop at y=10mm (3 segments)
N5 x=0 y=10 z=0
N6 x=50 y=10 z=0
N7 x=50 y=10 z=50
N8 x=0 y=10 z=50

E5 N5 N6 w=1 h=1
E6 N6 N7 w=1 h=1
E7 N7 N8 w=1 h=1

.external N1 N4
.external N5 N8

.freq fmin=1e6 fmax=1e6 ndec=1
.end
""")

summary = parser.get_summary()
check(summary['n_nodes'] == 8, "8 nodes parsed",
      "got %d" % summary['n_nodes'])
check(summary['n_segments'] == 6, "6 segments parsed",
      "got %d" % summary['n_segments'])
check(len(parser.ports) == 2, "2 ports parsed",
      "got %d" % len(parser.ports))

fh_builder = parser.to_peec_builder()
fh_topo = fh_builder.build_topology()
fh_solver = PEECCircuitSolver(fh_topo)

fh_result = fh_solver.compute_coupling_coefficient(1e6)

check(fh_result['Z_matrix'].shape == (2, 2), "2x2 Z-matrix from FastHenry")
check(fh_result['k'] > 0, "Positive coupling",
      "k=%.4f" % fh_result['k'])
check(fh_result['k'] < 1, "k < 1",
      "k=%.4f" % fh_result['k'])

print("  INFO: k=%.4f, L1=%.2f nH, L2=%.2f nH, M=%.2f nH" %
      (fh_result['k'], fh_result['L1'] * 1e9, fh_result['L2'] * 1e9,
       fh_result['M'] * 1e9))

# =========================================================================
# Test 7: Two coaxial square loops (WPT geometry)
# =========================================================================
print("\n--- Test 7: Coaxial square loops (WPT-like) ---")

# Tx loop: U-shaped 50mm square at z=0 (3 segments, open at one side)
# Rx loop: U-shaped 50mm square at z=gap
gap = 0.02  # 20mm gap
side = 0.05  # 50mm side

builder7 = PEECBuilder()

# Tx loop (z=0): N0->N1->N2->N3 (open loop, 3 segments)
tx_n = []
coords_tx = [
    [0, 0, 0], [side, 0, 0], [side, side, 0], [0, side, 0]
]
for c in coords_tx:
    tx_n.append(builder7.add_node_at(c[0], c[1], c[2]))

for i in range(3):  # Only 3 segments (leave one side open)
    builder7.add_connected_segment(tx_n[i], tx_n[i + 1], w, w, sigma)

# Rx loop (z=gap): N4->N5->N6->N7 (open loop, 3 segments)
rx_n = []
coords_rx = [
    [0, 0, gap], [side, 0, gap], [side, side, gap], [0, side, gap]
]
for c in coords_rx:
    rx_n.append(builder7.add_node_at(c[0], c[1], c[2]))

for i in range(3):
    builder7.add_connected_segment(rx_n[i], rx_n[i + 1], w, w, sigma)

builder7.add_port(tx_n[0], tx_n[3])  # Tx port: open ends
builder7.add_port(rx_n[0], rx_n[3])  # Rx port: open ends

topo7 = builder7.build_topology()
solver7 = PEECCircuitSolver(topo7)
result7 = solver7.compute_coupling_coefficient(1e6)

check(result7['k'] > 0, "WPT: positive coupling",
      "k=%.4f" % result7['k'])
check(result7['k'] < 1, "WPT: k < 1",
      "k=%.4f" % result7['k'])

# L1 ~= L2 for identical loops
L_sym = abs(result7['L1'] - result7['L2']) / max(result7['L1'], 1e-30)
check(L_sym < 0.05, "WPT: L1 ~= L2 (%.1f%% diff)" % (L_sym * 100),
      "L1=%.2e, L2=%.2e" % (result7['L1'], result7['L2']))

# Z12 = Z21 reciprocity
Z7 = result7['Z_matrix']
r7 = abs(Z7[0, 1] - Z7[1, 0]) / abs(Z7[0, 1]) if abs(Z7[0, 1]) > 1e-30 else 0.0
check(r7 < 1e-6, "WPT: Z12 = Z21 reciprocity (%.2e%%)" % (r7 * 100))

print("  INFO: gap=%dmm, k=%.4f, L1=%.2f nH, L2=%.2f nH, M=%.2f nH" %
      (gap * 1e3, result7['k'],
       result7['L1'] * 1e9, result7['L2'] * 1e9, result7['M'] * 1e9))

# =========================================================================
# Test 8: Frequency sweep multiport
# =========================================================================
print("\n--- Test 8: Frequency sweep multiport ---")

freqs = np.logspace(3, 8, 20)
sweep = solver.frequency_sweep_multiport(freqs)

check(sweep['Z_matrix'].shape == (20, 2, 2), "Z_matrix shape (20, 2, 2)",
      "got %s" % str(sweep['Z_matrix'].shape))
check(sweep['k'].shape == (20,), "k array shape (20,)")
check(all(sweep['k'] > 0), "All k > 0 in sweep")
check(all(sweep['k'] < 1), "All k < 1 in sweep")

# k should be roughly constant (frequency-independent for DC inductance)
k_mean = np.mean(sweep['k'])
k_std = np.std(sweep['k'])
k_variation = k_std / k_mean if k_mean > 0 else 0
check(k_variation < 0.1, "k roughly frequency-independent (CV=%.2f%%)" % (k_variation * 100),
      "k_mean=%.4f, k_std=%.4e" % (k_mean, k_std))

print("  INFO: k_mean=%.4f, k_std=%.2e" % (k_mean, k_std))

# =========================================================================
# Test 9: Multi-filament coupling
# =========================================================================
print("\n--- Test 9: Multi-filament coupling ---")

builder9 = PEECBuilder()
na = builder9.add_node_at(0, 0, 0)
nb = builder9.add_node_at(0.1, 0, 0)
nc = builder9.add_node_at(0, 0.005, 0)
nd = builder9.add_node_at(0.1, 0.005, 0)

# Multi-filament segments (3x3)
builder9.add_connected_segment(na, nb, 3e-3, 3e-3, sigma, nwinc=3, nhinc=3)
builder9.add_connected_segment(nc, nd, 3e-3, 3e-3, sigma, nwinc=3, nhinc=3)

builder9.add_port(na, nb)
builder9.add_port(nc, nd)

topo9 = builder9.build_topology()
solver9 = PEECCircuitSolver(topo9)
result9 = solver9.compute_coupling_coefficient(1e6)

check(result9['k'] > 0, "Multi-filament: k > 0",
      "k=%.4f" % result9['k'])
check(result9['k'] < 1, "Multi-filament: k < 1",
      "k=%.4f" % result9['k'])

print("  INFO: 3x3 filaments, k=%.4f, L1=%.2f nH, M=%.2f nH" %
      (result9['k'], result9['L1'] * 1e9, result9['M'] * 1e9))

# =========================================================================
# Test 10: Decoupled wires (large separation -> k ~= 0)
# =========================================================================
print("\n--- Test 10: Decoupled wires (large separation) ---")

builder10 = PEECBuilder()
na = builder10.add_node_at(0, 0, 0)
nb = builder10.add_node_at(0.01, 0, 0)  # Short 10mm wire
nc = builder10.add_node_at(0, 10.0, 0)  # 10m away!
nd = builder10.add_node_at(0.01, 10.0, 0)

builder10.add_connected_segment(na, nb, w, w, sigma)
builder10.add_connected_segment(nc, nd, w, w, sigma)
builder10.add_port(na, nb)
builder10.add_port(nc, nd)

topo10 = builder10.build_topology()
solver10 = PEECCircuitSolver(topo10)
result10 = solver10.compute_coupling_coefficient(1e6)

check(result10['k'] < 0.01, "Far wires: k < 0.01 (nearly decoupled)",
      "k=%.6f" % result10['k'])

print("  INFO: 10m separation, k=%.6f" % result10['k'])

# =========================================================================
# Summary
# =========================================================================
print("\n" + "=" * 70)
total = n_pass + n_fail
print("Results: %d/%d PASSED, %d FAILED" % (n_pass, total, n_fail))
if n_fail == 0:
    print("ALL TESTS PASSED")
else:
    print("SOME TESTS FAILED")
print("=" * 70)

sys.exit(n_fail)
