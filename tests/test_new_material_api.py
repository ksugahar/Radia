"""
Test script for material API:
- MatLin(mu_r) - isotropic linear material (mu_r = relative permeability)
- MatLin([mu_r_par, mu_r_perp], [ex, ey, ez]) - anisotropic linear material
- MatPM(Br, Hc, [mx, my, mz]) - permanent magnet with demagnetization
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import radia as rad
import numpy as np

print("="*60)
print("Testing Material API")
print("="*60)

# Test 1: Isotropic linear material
print("\n1. Testing MatLin(mu_r) - Isotropic Linear Material")
print("-"*60)
try:
	mu_r = 1000  # Soft iron
	mat_iso = rad.MatLin(mu_r)
	print(f"[OK] Created isotropic material with mu_r={mu_r}")
	print(f"  Material handle: {mat_iso}")

	# Apply to a block
	block = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,0])
	rad.MatApl(block, mat_iso)
	print(f"[OK] Applied material to block {block}")
except Exception as e:
	print(f"[ERROR] Error: {e}")

# Test 2: Anisotropic linear material with easy axis
print("\n2. Testing MatLin([mu_r_par, mu_r_perp], [ex, ey, ez]) - Anisotropic")
print("-"*60)
try:
	mu_r_par = 1.06
	mu_r_perp = 1.17
	easy_axis = [0, 0, 1]  # z-direction
	mat_aniso = rad.MatLin([mu_r_par, mu_r_perp], easy_axis)
	print(f"[OK] Created anisotropic material")
	print(f"  mu_r_par={mu_r_par}, mu_r_perp={mu_r_perp}")
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
	print(f"  Recoil permeability: mu_rec = {Br/(1.25663706212e-6 * Hc):.4f}")

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

# Test 4: mu_r validation (must be >= 1)
print("\n4. Testing mu_r >= 1 validation")
print("-"*60)
try:
	mat_bad = rad.MatLin(0.5)  # mu_r < 1 should fail
	print(f"[ERROR] Should have raised error for mu_r < 1")
except Exception as e:
	print(f"[OK] Correctly rejected mu_r < 1: {e}")

print("\n" + "="*60)
print("All Tests Completed!")
print("="*60)
