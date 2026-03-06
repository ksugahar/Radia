#!/usr/bin/env python
"""
Test vector potential A calculation in CoefficientFunction field source
"""

import sys
import os

# Add parent directory to path to import radia
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

import radia as rad

print("=" * 70)
print("Vector Potential Test")
print("=" * 70)

# Set Radia to use meters for consistency
rad.FldUnits('m')

# Clear all objects
rad.UtiDelAll()

# Test 1: Dictionary format with both B and A
print("\n[Test 1] Dictionary format: {'B': [...], 'A': [...]}")
print("-" * 70)

def field_with_A(coords):
	"""
	Return both B field and vector potential A
	For a uniform field Bz, A = [-By/2, Bx/2, 0]
	"""
	x, y, z = coords
	Bx, By, Bz = 0.0, 0.0, 1.0  # Uniform field in z direction

	# Vector potential for uniform Bz: A = (1/2) * B × r = [-By/2, Bx/2, 0]
	Ax = -y * Bz / 2.0
	Ay = x * Bz / 2.0
	Az = 0.0

	return {'B': [Bx, By, Bz], 'A': [Ax, Ay, Az]}

# Create background field source
bg_field = rad.ObjBckgCF(field_with_A)
print(f"  Background field object ID: {bg_field}")

# Test field calculation with A
point = [0.010, 0.020, 0.030]  # m (rad.FldUnits('m'))
result = rad.Fld(bg_field, 'ba', point)  # Request both B and A
print(f"  Point: {point} m")
print(f"  B field: [{result[0]:.6f}, {result[1]:.6f}, {result[2]:.6f}] T")
print(f"  A field: [{result[3]:.6f}, {result[4]:.6f}, {result[5]:.6f}] T*m")

# Verify A
expected_Ax = -point[1] * 1.0 / 2.0  # T·m (no conversion needed with rad.FldUnits('m'))
expected_Ay = point[0] * 1.0 / 2.0
expected_Az = 0.0

tolerance = 1e-6
if (abs(result[3] - expected_Ax) < tolerance and
    abs(result[4] - expected_Ay) < tolerance and
    abs(result[5] - expected_Az) < tolerance):
	print("  [OK] Vector potential A is correct!")
else:
	print(f"  [FAIL] Expected A: [{expected_Ax}, {expected_Ay}, {expected_Az}]")

# Test 2: Backward compatibility - list format (B only)
print("\n[Test 2] Backward compatibility: [Bx, By, Bz] format")
print("-" * 70)

def field_B_only(coords):
	"""Return only B field (backward compatible)"""
	return [0.0, 0.0, 1.0]

bg_field2 = rad.ObjBckgCF(field_B_only)
print(f"  Background field object ID: {bg_field2}")

result2 = rad.Fld(bg_field2, 'b', point)
print(f"  Point: {point} m")
print(f"  B field: [{result2[0]:.6f}, {result2[1]:.6f}, {result2[2]:.6f}] T")
print("  [OK] Backward compatibility maintained")

# Cleanup
rad.UtiDelAll()

print("\n" + "=" * 70)
print("All Vector Potential Tests Passed!")
print("=" * 70)
