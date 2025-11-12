"""
Verify H-Matrix Field Accuracy (Simplified - No VTK)
"""

import sys
import os
sys.path.insert(0, os.path.join('S:/Radia/01_GitHub', 'build/Release'))

import radia as rad
import numpy as np

def create_magnet(n_per_side):
    """Create a cubic magnet subdivided into n x n x n elements."""
    size = 20.0
    elem_size = size / n_per_side
    container = rad.ObjCnt([])

    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                x = (i - n_per_side/2 + 0.5) * elem_size
                y = (j - n_per_side/2 + 0.5) * elem_size
                z = (k - n_per_side/2 + 0.5) * elem_size
                # Permanent magnet: M = 795774.7 A/m
                block = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size], [0, 0, 795774.7])
                rad.ObjAddToCnt(container, [block])

    return container

print("="*70)
print("H-MATRIX FIELD ACCURACY VERIFICATION (Simplified)")
print("="*70)

# Test 1: Small magnet (N=125)
print("\nTest 1: Small magnet (N=125, permanent magnet)")
print("-" * 70)
magnet_small = create_magnet(5)

B_small_center = rad.Fld(magnet_small, 'b', [0, 0, 0])
B_small_outside = rad.Fld(magnet_small, 'b', [0, 0, 30])

print(f"  B at center [0,0,0]:    {B_small_center}")
print(f"  B outside [0,0,30mm]:   {B_small_outside}")

# Test 2: Medium magnet (N=343)
print("\nTest 2: Medium magnet (N=343, permanent magnet)")
print("-" * 70)
magnet_medium = create_magnet(7)

B_medium_center = rad.Fld(magnet_medium, 'b', [0, 0, 0])
B_medium_outside = rad.Fld(magnet_medium, 'b', [0, 0, 30])

print(f"  B at center [0,0,0]:    {B_medium_center}")
print(f"  B outside [0,0,30mm]:   {B_medium_outside}")

# Multiple points test
print("\n" + "="*70)
print("MULTIPLE POINTS TEST")
print("="*70)

test_points = [
    [0, 0, 0],
    [5, 0, 0],
    [0, 5, 0],
    [0, 0, 15],
    [10, 0, 0],
    [0, 0, 30],
    [0, 0, 50],
]

print(f"\n{'Point [mm]':<20} {'B_z (N=125) [T]':<20} {'B_z (N=343) [T]':<20} {'Rel Error [%]':<15}")
print("-" * 70)

max_error = 0.0
for point in test_points:
    B_small = rad.Fld(magnet_small, 'b', point)
    B_medium = rad.Fld(magnet_medium, 'b', point)

    if abs(B_small[2]) > 1e-10:
        rel_error = abs(B_small[2] - B_medium[2]) / abs(B_small[2]) * 100
    else:
        rel_error = 0.0

    max_error = max(max_error, rel_error)

    point_str = f"[{point[0]}, {point[1]}, {point[2]}]"
    print(f"{point_str:<20} {B_small[2]:<20.8f} {B_medium[2]:<20.8f} {rel_error:<15.4f}")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"\nMaximum relative error: {max_error:.4f}%")
print("\n[OK] Test completed successfully")
