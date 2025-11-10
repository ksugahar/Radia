#!/usr/bin/env python
"""
Test ObjBckgCF alone (no magnetic material)

Tests the callback-based background field source in isolation
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../build/Release'))

import radia as rd

print("=" * 70)
print("ObjBckgCF Alone - B->H Conversion Test")
print("=" * 70)

# Simple uniform field callback (returns B in Tesla)
def uniform_field_callback():
	"""Returns uniform B = [1, 0, 0] T"""
	call_count = [0]
	def field(pos):
		call_count[0] += 1
		result = [1.0, 0.0, 0.0]  # 1 T in X direction
		if call_count[0] <= 2:
			print(f"  [Callback #{call_count[0]}] pos={pos} -> B={result} T")
		return result
	return field

# Create ObjBckgCF
callback = uniform_field_callback()
bckg_cf = rd.ObjBckgCF(callback)
print(f"\nObjBckgCF created")

# Test point
pt = [0, 0, 0]

# Get B and H fields FROM ObjBckgCF ALONE
print(f"\nCalling Fld on ObjBckgCF source alone...")
B = rd.Fld(bckg_cf, 'b', pt)
H = rd.Fld(bckg_cf, 'h', pt)

print(f"\nAt point {pt}:")
print(f"  B = {B} T")
print(f"  H = {H} A/m")

# Expected values
mu_0 = 1.25663706212e-6
H_expected = 1.0 / mu_0

print(f"\nExpected:")
print(f"  H_x = {H_expected} A/m")

# Check ratio
if abs(H[0]) > 1e-10:
	ratio = B[0] / H[0]
	print(f"\nB_x / H_x = {ratio:.15e} T/(A/m)")
	print(f"mu_0      = {mu_0:.15e} T/(A/m)")
	print(f"Relative error: {abs(ratio - mu_0)/mu_0 * 100:.10f}%")

	if abs(ratio - mu_0) / mu_0 < 1e-6:
		print(f"\n[OK] B/H = mu_0 (conversion working)")
	else:
		print(f"\n[ERROR] B/H != mu_0 (conversion NOT working)")

		if abs(B[0] - H[0]) < 1e-10:
			print(f"[ERROR] H equals B! No conversion applied!")
else:
	print("\n[ERROR] H_x is zero!")

rd.UtiDelAll()
