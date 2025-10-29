#!/usr/bin/env python
"""
Simplified test of 2024_02_01_振り分け電磁石.py
Tests only coil creation without Cubit dependency
"""

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))

import radia as rad
from numpy import *
from COIL import cCOIL

print("="*60)
print("Simplified test: Septum magnet coil creation")
print("="*60)

# Test with single current value
I_ = 130
COIL = []
I = 64*I_

print(f"\nTesting with current I = {I} A")

# Create coil segments based on original script
V = array([[1,0,0], [0,1,0], [0,0,1]]).T
W = 122
H = 122
L = 16.43186645*2
x0 = array([(48+170)/2,-L/2,-20-W/2])
R = []
phi = 0
tilt = 0

print("\n1. Creating straight section...")
COIL.append(cCOIL('GCE',I,x0,V,W,H,R,phi,tilt,L))
print(f"   Position: {COIL[-1].c}")

x0 = COIL[-1].x1
W  = COIL[-1].W
H  = COIL[-1].H
V  = COIL[-1].V1
R = 60+61
phi = 64.59228189
tilt = 90
L = 0.0

print("\n2. Creating first arc...")
COIL.append(cCOIL('ARC',I,x0,V,W,H,R,phi,tilt,L))
print(f"   Position: {COIL[-1].c}, R={R}, phi={phi}")

x0 = COIL[-1].x1
W  = COIL[-1].W
H  = COIL[-1].H
V  = COIL[-1].V1
R = []
phi = 0
tilt = 90
L = 1018.51313197

print("\n3. Creating long straight section...")
COIL.append(cCOIL('GCE',I,x0,V,W,H,R,phi,tilt,L))
print(f"   Position: {COIL[-1].c}")

print(f"\nTotal coil segments created: {len(COIL)}")

# Create Radia objects from COIL segments
print("\nConverting to Radia objects...")
coils = []
for n in range(len(COIL)):
	if (COIL[n].type=='ARC'):
		phi1 = COIL[n].rot[0]/180*pi
		if (phi1 <= 0):
			phi1 = phi1 + 2*pi
		phi2 = (COIL[n].rot[0]+COIL[n].phi)/180*pi
		if (phi1 > phi2) or (phi2 <= 0):
			phi2 = phi2 + 2*pi
		coil = rad.ObjArcCur([0,0,0],[COIL[n].R-COIL[n].W/2,COIL[n].R+COIL[n].W/2], [phi1,phi2], COIL[n].H, 10, COIL[n].I/COIL[n].W/COIL[n].H, "auto")
		trf = rad.TrfRot([0,0,0], [0,0,1], COIL[n].rot[2]/180*pi)
		trf = rad.TrfCmbR(trf,rad.TrfRot([0,0,0], [1,0,0], COIL[n].rot[1]/180*pi))
		print(f"  Segment {n+1}: ARC - Radia object created")
	else:
		coil = rad.ObjRecCur([0,0,0],[COIL[n].W, COIL[n].L, COIL[n].H], [0,COIL[n].I/COIL[n].W/COIL[n].H,0])
		trf = rad.TrfRot([0,0,0], [0,0,1], COIL[n].rot[2]/180*pi)
		trf = rad.TrfCmbR(trf,rad.TrfRot([0,0,0], [1,0,0], COIL[n].rot[1]/180*pi))
		trf = rad.TrfCmbR(trf,rad.TrfRot([0,0,0], [0,0,1], COIL[n].rot[0]/180*pi))
		print(f"  Segment {n+1}: STRAIGHT - Radia object created")
	trf = rad.TrfCmbL(trf,rad.TrfTrsl([COIL[n].c[0],COIL[n].c[1],COIL[n].c[2]]))
	coil = rad.TrfOrnt(coil,trf)
	coils.append(coil)

coils = rad.ObjCnt(coils)
rad.ObjDrwAtr(coils, [1,0,0], 0.1)
rad.TrfZerPara(coils, [0,0,0], [1,0,0])

print(f"\nCoils container object: ID {coils}")

# Calculate field at a few points
test_points = [[0,0,0], [0,0,600], [100,0,0]]
print("\nField calculations:")
for pt in test_points:
	B = rad.Fld(coils, "BxByBz", pt)
	print(f"  Point {pt}: Bx={B[0]:.6e}, By={B[1]:.6e}, Bz={B[2]:.6e} T")

print()
print("="*60)
print("SUCCESS: Septum magnet coils created successfully!")
print("="*60)
print("\nNote: This is a simplified test without the yoke (iron core).")
print("The full simulation requires Cubit for mesh generation.")
