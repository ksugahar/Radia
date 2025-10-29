#!/usr/bin/env python
"""
Test 2024_01_30_Radiaの練習_01.py
Modified version without OpenGL and GUI dependencies
"""

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))

import radia as rad
from numpy import *

print("="*60)
print("Practice 01: Arc current with rectangular magnets")
print("="*60)

# Create arc current
p = []
p.append(rad.ObjArcCur([0,0,0], [95, 105], [0,2*pi], 10, 36, 1e6/10/10, "man"))
g1 = rad.ObjCnt(p)
print(f"Arc current created: ID {g1}")

# Create rectangular magnet
M = [0,0,0]
g2 = rad.ObjRecMag([0,0,0], [10,10,5], [0,0,1])
mat2 = rad.MatLin([1,1], [0, 0, 1])
rad.MatApl(g2, mat2)
print(f"Rectangular magnet created: ID {g2}")

# Create multiple magnets
g3 = []
z = arange(5, 20, 5)
for n in range(len(z)):
	g3.append(rad.ObjRecMag([0,0,z[n]],[10,10,5],[0,0,1]))

print(f"Multiple magnets: {len(g3)} magnets created")
g3 = rad.ObjCnt(g3)

rad.TrfZerPerp(g3, [0,0,0], [0,0,1])
rad.ObjDrwAtr(g3, [0,1,1], 0.1)
mat3 = rad.MatLin([1,1], [0, 0, 1])
rad.MatApl(g3, mat3)

# Combine all objects
g = rad.ObjCnt([g1, g2, g3])
print(f"Container object: ID {g}")

# Solve
print("\nSolving...")
res = rad.Solve(g, 0.01, 1000, 4)
print(f"Solve result: {res}")

# Calculate field along z-axis
x1 =    0
y1 =    0
z1 = -100
x2 =    0
y2 =    0
z2 =  100
np =  100

print(f"\nCalculating field from z={z1} to z={z2} mm ({np} points)")
Bz = rad.FldLst(g,"Bz", [x1,y1,z1], [x2,y2,z2], int(np), "noarg")

print(f"Bz min: {min(Bz):.6e} T")
print(f"Bz max: {max(Bz):.6e} T")
print(f"Bz at origin: {Bz[np//2]:.6e} T")

# Analytical comparison
z_arr = linspace(z1,z2,np)
mu0 = 4*pi*1e-7*1e6
a = 0.1
I = 1
z_m = z_arr/1000
Bz_analytical = mu0*I*a**2/2/(a**2+z_m**2)**(3/2)

print(f"Analytical Bz at origin: {Bz_analytical[np//2]:.6e} T")

print()
print("="*60)
print("SUCCESS: Practice 01 completed!")
print("="*60)
