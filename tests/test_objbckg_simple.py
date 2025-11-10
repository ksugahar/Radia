#!/usr/bin/env python
"""
Test ObjBckg B->H conversion

Verifies that uniform background field (ObjBckg) correctly converts B (Tesla) to H (A/m)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../build/Release'))

import radia as rd

print("=" * 70)
print("ObjBckg B->H Conversion Test")
print("=" * 70)

# Create uniform background field
B_bg = [1.0, 0, 0]  # 1 Tesla in X direction
bg = rd.ObjBckg(B_bg)

print(f"\nBackground field: B = {B_bg} T")

# Test point
pt = [0, 0, 0]

# Get B and H fields
B = rd.Fld(bg, 'b', pt)
H = rd.Fld(bg, 'h', pt)

print(f"\nAt point {pt}:")
print(f"  B = {B} T")
print(f"  H = {H} A/m")

# Calculate expected H
mu_0 = 1.25663706212e-6  # T/(A/m)
H_expected = [B_bg[i] / mu_0 for i in range(3)]

print(f"\nExpected:")
print(f"  H = {H_expected} A/m")

# Check B/H ratio
if abs(H[0]) > 1e-10:
	ratio = B[0] / H[0]
	print(f"\nB_x / H_x = {ratio:.15e} T/(A/m)")
	print(f"mu_0      = {mu_0:.15e} T/(A/m)")
	print(f"Relative error: {abs(ratio - mu_0)/mu_0 * 100:.10f}%")

	if abs(ratio - mu_0) / mu_0 < 1e-6:
		print(f"\n[OK] B/H = mu_0 (conversion working correctly)")
	else:
		print(f"\n[ERROR] B/H != mu_0 (conversion NOT working)")
else:
	print("\n[ERROR] H_x is zero!")

rd.UtiDelAll()
