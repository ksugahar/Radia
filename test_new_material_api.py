"""
Test script for new material API:
- MatLin(ksi) - isotropic linear material
- MatLin([ksi_par, ksi_perp], [ex, ey, ez]) - anisotropic linear material
- MatPM(Br, Hc, [mx, my, mz]) - permanent magnet with demagnetization
"""

import sys
import os
sys.path.insert(0, os.path.join('S:/Radia/01_GitHub', 'build/Release'))

import radia as rad
import numpy as np

print("="*60)
print("Testing New Material API")
print("="*60)

# Test 1: Isotropic linear material
print("\n1. Testing MatLin(ksi) - Isotropic Linear Material")
print("-"*60)
try:
	ksi = 1000  # Soft iron
	mat_iso = rad.MatLin(ksi)
	print(f"[OK] Created isotropic material with ksi={ksi}")
	print(f"  Material handle: {mat_iso}")

	# Apply to a block
	block = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,0])
	rad.MatApl(block, mat_iso)
	print(f"[OK] Applied material to block {block}")
except Exception as e:
	print(f"[ERROR] Error: {e}")

# Test 2: Anisotropic linear material with easy axis
print("\n2. Testing MatLin([ksi_par, ksi_perp], [ex, ey, ez]) - Anisotropic")
print("-"*60)
try:
	ksi_par = 0.06
	ksi_perp = 0.17
	easy_axis = [0, 0, 1]  # z-direction
	mat_aniso = rad.MatLin([ksi_par, ksi_perp], easy_axis)
	print(f"[OK] Created anisotropic material")
	print(f"  ksi_par={ksi_par}, ksi_perp={ksi_perp}")
	print(f"  Easy axis: {easy_axis}")
	print(f"  Material handle: {mat_aniso}")

	# Apply to a block
	block2 = rad.ObjRecMag([20,0,0], [10,10,10], [0,0,0])
	rad.MatApl(block2, mat_aniso)
	print(f"[OK] Applied material to block {block2}")
except Exception as e:
	print(f"[ERROR] Error: {e}")

# Test 3: Permanent magnet with demagnetization curve
print("\n3. Testing MatPM(Br, Hc, [mx, my, mz]) - Permanent Magnet")
print("-"*60)
try:
	Br = 1.43  # Tesla (NdFeB N52)
	Hc = 876000  # A/m
	mag_axis = [0, 0, 1]  # z-direction
	mat_pm = rad.MatPM(Br, Hc, mag_axis)
	print(f"[OK] Created permanent magnet material (NdFeB N52)")
	print(f"  Br={Br} T, Hc={Hc} A/m")
	print(f"  Magnetization axis: {mag_axis}")
	print(f"  Material handle: {mat_pm}")
	print(f"  Recoil permeability: μ_rec = {Br/(1.25663706212e-6 * Hc):.4f}")

	# Apply to a block
	block3 = rad.ObjRecMag([40,0,0], [20,20,10], [0,0,0])
	rad.MatApl(block3, mat_pm)
	print(f"[OK] Applied material to block {block3}")

	# Solve and check field
	container = rad.ObjCnt([block3])
	rad.Solve(container, 0.0001, 1000)
	B = rad.Fld(container, 'b', [40,0,15])
	print(f"[OK] Magnetic field at [40,0,15]: B = {B} T")
	print(f"  |B| = {np.linalg.norm(B):.4f} T")
except Exception as e:
	print(f"[ERROR] Error: {e}")

# Test 4: Backward compatibility - old MatLin API
print("\n4. Testing Backward Compatibility - Old MatLin API")
print("-"*60)
try:
	# Old form: MatLin([ksi_par, ksi_perp], Mr_scalar)
	mat_old = rad.MatLin([0.06, 0.17], 1e6)
	print(f"[OK] Old API works: MatLin([0.06, 0.17], 1e6)")
	print(f"  Material handle: {mat_old}")
except Exception as e:
	print(f"[ERROR] Error: {e}")

print("\n" + "="*60)
print("All Tests Completed!")
print("="*60)
