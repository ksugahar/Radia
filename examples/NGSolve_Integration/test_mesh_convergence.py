#!/usr/bin/env python
"""
Mesh convergence test for curl(A) = B verification

Tests different mesh resolutions to verify that error decreases
as mesh becomes finer.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "build", "Release"))

import radia as rad
try:
	from ngsolve import *
	from netgen.occ import *
	import rad_ngsolve
except ImportError:
	print("ERROR: NGSolve not available")
	sys.exit(1)

import numpy as np

print("=" * 80)
print("MESH CONVERGENCE TEST: curl(A) = B")
print("=" * 80)

# Create Radia magnet
rad.UtiDelAll()
rad.FldUnits('m')  # Set units to meters
magnet = rad.ObjRecMag([0, 0, 0], [0.040, 0.040, 0.060], [0, 0, 1.2])

# Create background field
def radia_field_with_A(coords):
	x, y, z = coords
	B = rad.Fld(magnet, 'b', [x, y, z])
	A = rad.Fld(magnet, 'a', [x, y, z])
	return {'B': list(B), 'A': list(A)}

bg_field = rad.ObjBckgCF(radia_field_with_A)

# Get CoefficientFunctions
A_cf = rad_ngsolve.RadiaField(bg_field, 'a')
B_cf = rad_ngsolve.RadiaField(bg_field, 'b')

# Test points (in meters)
test_points = [
	(0.030, 0.020, 0.030),
	(0.030, 0.020, 0.040),
	(0.030, 0.020, 0.050),
	(0.020, 0.020, 0.040),
	(0.040, 0.040, 0.050),
]

# Different mesh sizes to test
mesh_sizes = [0.016, 0.008, 0.004, 0.002, 0.001]  # mm: 16mm, 8mm, 4mm, 2mm, 1mm

print("\nTesting different mesh resolutions:")
print("-" * 80)
print(f"{'Mesh size (m)':<15} {'# Elements':<12} {'# DOFs':<12} {'Mean Error (T)':<18} {'Max Error (T)':<18}")
print("-" * 80)

results = []

for maxh in mesh_sizes:
	# Create mesh
	box = Box((0.01, 0.01, 0.02), (0.06, 0.06, 0.08))
	geo = OCCGeometry(box)
	mesh = Mesh(geo.GenerateMesh(maxh=maxh))

	# Create finite element spaces
	fes_hcurl = HCurl(mesh, order=2)

	# Project A and compute curl
	A_gf = GridFunction(fes_hcurl)
	A_gf.Set(A_cf)
	curl_A_gf = curl(A_gf)

	# Evaluate at test points
	errors = []
	for point in test_points:
		try:
			mip = mesh(*point)
			curl_A_val = np.array(curl_A_gf(mip))
			B_val = np.array(B_cf(mip))
			error_norm = np.linalg.norm(curl_A_val - B_val)
			errors.append(error_norm)
		except:
			pass

	if len(errors) > 0:
		mean_error = np.mean(errors)
		max_error = np.max(errors)
		results.append({
			'maxh': maxh,
			'ne': mesh.ne,
			'ndof': fes_hcurl.ndof,
			'mean_error': mean_error,
			'max_error': max_error
		})

		print(f"{maxh:<15.6f} {mesh.ne:<12} {fes_hcurl.ndof:<12} {mean_error:<18.6e} {max_error:<18.6e}")

print("-" * 80)

# Analyze convergence
print("\nConvergence Analysis:")
print("-" * 80)

if len(results) >= 2:
	print("\nError reduction rates:")
	for i in range(1, len(results)):
		ratio_h = results[i-1]['maxh'] / results[i]['maxh']
		ratio_mean = results[i-1]['mean_error'] / results[i]['mean_error']
		ratio_max = results[i-1]['max_error'] / results[i]['max_error']

		print(f"  Mesh refined by {ratio_h:.1f}x:")
		print(f"    Mean error reduced by {ratio_mean:.2f}x")
		print(f"    Max error reduced by {ratio_max:.2f}x")

		if ratio_mean > 1.1:
			print(f"    [OK] Mean error decreasing")
		else:
			print(f"    [!] Mean error NOT decreasing significantly")

		if ratio_max > 1.1:
			print(f"    [OK] Max error decreasing")
		else:
			print(f"    [!] Max error NOT decreasing significantly")
		print()

	# Overall verdict
	all_decreasing_mean = all(
		results[i]['mean_error'] < results[i-1]['mean_error'] * 0.9
		for i in range(1, len(results))
	)
	all_decreasing_max = all(
		results[i]['max_error'] < results[i-1]['max_error'] * 0.9
		for i in range(1, len(results))
	)

	print("=" * 80)
	if all_decreasing_mean and all_decreasing_max:
		print("RESULT: [CONVERGING]")
		print("\nErrors decrease consistently as mesh is refined.")
		print("This confirms that curl(A) = B is correctly implemented.")
	elif all_decreasing_mean or all_decreasing_max:
		print("RESULT: [PARTIALLY CONVERGING]")
		print("\nSome errors decrease with mesh refinement.")
	else:
		print("RESULT: [NOT CONVERGING]")
		print("\nErrors do not decrease with mesh refinement.")
		print("This may indicate an implementation error.")
	print("=" * 80)

rad.UtiDelAll()
