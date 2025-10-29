#!/usr/bin/env python
"""
Test 2024_01_30_Radiaの練習_02.py
Modified version without OpenGL and GUI dependencies
"""

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))

import radia as rad
from numpy import *

print("="*60)
print("Practice 02: Particle trajectory in arc current field")
print("="*60)

# Create arc current with larger height
p = []
p.append(rad.ObjArcCur([0,0,0], [95, 105], [0,2*pi], 2000, 36, 1e6/2000/10, "man"))
g1 = rad.ObjCnt(p)
print(f"Arc current created: ID {g1}")

g = rad.ObjCnt([g1])
print(f"Container object: ID {g}")

# Calculate field along x-axis
x1 =  -50
y1 =    0
z1 =    0
x2 =   50
y2 =    0
z2 =    0
np_points =  100

print(f"\nCalculating field from x={x1} to x={x2} mm ({np_points} points)")
Bz = rad.FldLst(g,"Bz", [x1,y1,z1], [x2,y2,z2], int(np_points), "noarg")

print(f"Bz min: {min(Bz):.6e} T")
print(f"Bz max: {max(Bz):.6e} T")
print(f"Bz at center: {Bz[np_points//2]:.6e} T")

# Set field computation precision
rad.FldCmpPrc('PrcB->1e-10, PrcCoord->1e-12')

# Calculate particle trajectory
print("\nCalculating particle trajectory...")
print("  Energy: 0.002 GeV")
print("  Initial position: [10.0, 0.0, 1.0, 0.0]")
print("  Parameter range: [0.0, 10.0]")
print("  Number of points: 50")

r = rad.FldPtcTrj(g, 0.002, [10.0, 0.0, 1.0, 0.0], [0.0, 10.0], 50)
r = array(r)

print(f"Trajectory calculated: {len(r)} points")
print(f"Start position: [{r[0][0]:.3f}, {r[0][1]:.3f}]")
print(f"End position: [{r[-1][0]:.3f}, {r[-1][1]:.3f}]")

# Sample a few points
print("\nTrajectory samples:")
for i in [0, len(r)//4, len(r)//2, 3*len(r)//4, len(r)-1]:
	print(f"  Point {i:2d}: x={r[i][0]:7.3f}, y={r[i][1]:7.3f}, z={r[i][2]:7.3f}")

print()
print("="*60)
print("SUCCESS: Practice 02 completed!")
print("="*60)
