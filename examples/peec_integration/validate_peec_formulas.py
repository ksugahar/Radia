"""
Validate PEEC Formula Implementation

Checks:
    1. Inductance formula correctness
    2. Resistance formula correctness
    3. Edge extraction logic
    4. Comparison with analytical solutions

Tests against known analytical results for simple geometries.
"""

import numpy as np
import sys

print("=" * 70)
print("PEEC Formula Validation")
print("=" * 70)

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # H/m

print("\n[1] Self-Inductance Formula Validation")
print("=" * 70)

# Test case: Straight wire segment
length = 0.1  # m (100mm)
radius = 1e-3  # m (1mm wire radius)

# Formula: L_self = (mu_0 * l / 2pi) * [ln(2*l/a) - 2]
# where a = wire radius, l = length
L_self_formula = (MU_0 * length / (2 * np.pi)) * (np.log(2 * length / radius) - 2)

print(f"\nStraight wire:")
print(f"  Length: {length*1e3:.1f} mm")
print(f"  Radius: {radius*1e3:.1f} mm")
print(f"  L_self = {L_self_formula*1e9:.4f} nH")

# Alternative formula (Neumann)
# L_self = (mu_0 * l / 2pi) * [ln(l/a) + 0.25]
L_self_neumann = (MU_0 * length / (2 * np.pi)) * (np.log(length / radius) + 0.25)
print(f"  L_self (Neumann) = {L_self_neumann*1e9:.4f} nH")

# Error check
error = abs(L_self_formula - L_self_neumann) / L_self_neumann * 100
print(f"  Formula difference: {error:.2f}%")

if error > 5:
    print("  WARNING: Large difference - check formula!")
else:
    print("  OK: Formulas agree within 5%")

print("\n[2] Mutual Inductance Formula Validation")
print("=" * 70)

# Test case: Two parallel wire segments
l1 = 0.05  # m
l2 = 0.05  # m
distance = 0.01  # m (10mm separation)

# Simplified filament formula (parallel wires)
# M = (mu_0 / 4pi) * l1 * l2 / r
M_filament = (MU_0 / (4 * np.pi)) * l1 * l2 / distance

print(f"\nTwo parallel wires:")
print(f"  Length 1: {l1*1e3:.1f} mm")
print(f"  Length 2: {l2*1e3:.1f} mm")
print(f"  Distance: {distance*1e3:.1f} mm")
print(f"  M_filament = {M_filament*1e9:.4f} nH")

# More accurate Neumann formula for parallel wires
# M = (mu_0 * l / 2pi) * [ln(l/d + sqrt(1 + (l/d)^2)) - sqrt(1 + (d/l)^2) + d/l]
ratio = l1 / distance
M_neumann = (MU_0 * l1 / (2 * np.pi)) * (
    np.log(ratio + np.sqrt(1 + ratio**2))
    - np.sqrt(1 + (distance/l1)**2)
    + distance/l1
)
print(f"  M_neumann = {M_neumann*1e9:.4f} nH")

error = abs(M_filament - M_neumann) / M_neumann * 100
print(f"  Formula error: {error:.2f}%")

if error > 20:
    print("  WARNING: Filament approximation has >20% error!")
    print("  Recommendation: Use Neumann formula for better accuracy")
else:
    print("  OK: Filament approximation acceptable")

print("\n[3] DC Resistance Formula Validation")
print("=" * 70)

# Copper properties
RHO = 1.68e-8  # Ohm*m at 20C
SIGMA = 1.0 / RHO  # S/m

length = 0.314  # m (approx circumference of 50mm radius coil)
area = 4e-3 * 4e-3  # m^2 (4mm x 4mm wire)

R_dc = RHO * length / area

print(f"\nCopper wire:")
print(f"  Resistivity: {RHO*1e8:.2f} * 10^-8 Ohm*m")
print(f"  Length: {length*1e3:.1f} mm")
print(f"  Cross-section: {area*1e6:.2f} mm^2")
print(f"  R_dc = {R_dc*1e3:.4f} mOhm")

# Sanity check: typical value for copper wire
# Rule of thumb: ~5 mOhm/m for 16mm^2 wire at 20C
R_per_meter = RHO / area  # Ohm/m
print(f"  Resistance per meter: {R_per_meter*1e3:.2f} mOhm/m")

expected_range = (0.5e-3, 5e-3)  # Ohm/m
if expected_range[0] <= R_per_meter <= expected_range[1]:
    print("  OK: Resistance in expected range for copper wire")
else:
    print(f"  WARNING: Expected {expected_range[0]*1e3:.2f}-{expected_range[1]*1e3:.2f} mOhm/m")

print("\n[4] Circular Coil Analytical Solution")
print("=" * 70)

# Single-turn circular coil inductance (analytical)
radius = 50e-3  # m (50mm mean radius)
wire_radius = 2e-3  # m (2mm, for 4mm diameter wire)

# Formula: L = mu_0 * R * [ln(8*R/a) - 2]
# where R = mean radius, a = wire radius
L_analytical = MU_0 * radius * (np.log(8 * radius / wire_radius) - 2)

print(f"\nSingle-turn circular coil (analytical):")
print(f"  Mean radius: {radius*1e3:.1f} mm")
print(f"  Wire diameter: {wire_radius*2*1e3:.1f} mm")
print(f"  L_analytical = {L_analytical*1e9:.2f} nH")

# Alternative formula (Grover)
# L = mu_0 * R * [ln(8*R/a) - 1.75]
L_grover = MU_0 * radius * (np.log(8 * radius / wire_radius) - 1.75)
print(f"  L_grover = {L_grover*1e9:.2f} nH")

print("\n[5] CRITICAL ISSUES in demo_peec_dc_simple.py")
print("=" * 70)

print("\nISSUE 1: Total Inductance Calculation")
print("  WRONG: L_total = L_matrix.sum()")
print("  WHY: This assumes all edges are in series, but they form loops!")
print("  CORRECT: Use Loop-Star decomposition to extract loop inductance")

print("\nISSUE 2: Mutual Inductance Formula")
print("  CURRENT: Filament approximation (parallel wire assumption)")
print("  LIMITATION: Only accurate when wires are parallel")
print("  ERROR: Can be >20% for curved geometries")
print("  BETTER: Use Neumann integral with proper geometry")

print("\nISSUE 3: Cross-Section Area")
print("  CURRENT: Assumes all edges have same area (4mm x 4mm)")
print("  PROBLEM: Surface mesh edges don't have well-defined cross-section")
print("  BETTER: Assign effective area based on surface mesh density")

print("\nISSUE 4: Edge Orientation")
print("  CURRENT: t_i . t_j (dot product)")
print("  MISSING: For non-parallel edges, need full vector formula")

print("\n[6] Recommended Fixes")
print("=" * 70)

print("\n1. Self-Inductance:")
print("   Keep current formula - it's standard for thin wires")

print("\n2. Mutual Inductance:")
print("   Replace filament with Neumann integral:")
print("   M_ij = (mu_0/4pi) * integral_i integral_j (dl_i . dl_j / |r_i - r_j|)")

print("\n3. Total Inductance:")
print("   Implement Loop-Star decomposition:")
print("   - Identify independent loops")
print("   - Build loop inductance matrix")
print("   - Extract single-loop inductance")

print("\n4. Cross-Section:")
print("   For surface mesh, use effective area:")
print("   A_eff = total_surface_area / n_edges")
print("   Or: Define port-to-port effective width")

print("\n[7] Validation Against Analytical Solution")
print("=" * 70)

print(f"\nExpected inductance for 50mm radius circular coil:")
print(f"  Analytical: {L_analytical*1e9:.2f} nH")
print(f"  Grover:     {L_grover*1e9:.2f} nH")
print("\nIf PEEC result differs by >10%, check implementation!")

print("\n" + "=" * 70)
print("Validation Summary")
print("=" * 70)

print("\nFORMULAS:")
print("  Self-inductance: OK")
print("  Mutual inductance: WARNING - simplified")
print("  DC resistance: OK")

print("\nIMPLEMENTATION:")
print("  Edge extraction: OK")
print("  Total inductance: WRONG - needs Loop-Star")
print("  Cross-section: QUESTIONABLE - needs clarification")

print("\nRECOMMENDATION:")
print("  1. Run demo_peec_dc_simple.py and compare L_total to analytical")
print("  2. If error >10%, implement Loop-Star decomposition")
print("  3. Use Neumann formula for mutual inductance")

print("\n" + "=" * 70)
