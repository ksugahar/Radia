#!/usr/bin/env python3
"""
Test: GridFunction.Set() vs .Interpolate()

Compare performance and behavior of Set() vs Interpolate() methods
for evaluating Radia field on NGSolve mesh.
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

import time
import numpy as np

print("=" * 80)
print("Test: GridFunction.Set() vs .Interpolate()")
print("=" * 80)

# Create magnet
rad.UtiDelAll()
n = 5
cube_size = 100.0  # mm
elem_size = cube_size / n
mag_value = 1.2  # T

print(f"\n[Setup] Creating magnet: {n}x{n}x{n} = {n**3} elements")
elements = []
for i in range(n):
	for j in range(n):
		for k in range(n):
			x = (i - n/2 + 0.5) * elem_size
			y = (j - n/2 + 0.5) * elem_size
			z = (k - n/2 + 0.5) * elem_size
			elem = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size],
			                      [0, 0, mag_value])
			elements.append(elem)

magnet = rad.ObjCnt(elements)

# Create NGSolve mesh
print("\n[Setup] Creating NGSolve mesh")
box_size_m = cube_size * 2.0 / 1000  # mm to m
geo = Box(
	Pnt(-box_size_m/2, -box_size_m/2, -box_size_m/2),
	Pnt( box_size_m/2,  box_size_m/2,  box_size_m/2)
)

# Test different mesh sizes
mesh_configs = [
	{"h": 0.05, "desc": "Coarse (h=0.05m)"},
	{"h": 0.025, "desc": "Medium (h=0.025m)"},
	{"h": 0.0125, "desc": "Fine (h=0.0125m)"},
]

print("\n" + "=" * 80)
print("PERFORMANCE COMPARISON: Set() vs Interpolate()")
print("=" * 80)

B_cf = rad_ngsolve.RadiaField(magnet, 'b')

for config in mesh_configs:
	h = config['h']
	desc = config['desc']

	print(f"\n{desc}")
	print("-" * 80)

	mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=h))
	print(f"  Mesh: {mesh.nv} vertices, {mesh.ne} elements")

	# Test with HCurl space
	print("\n  [HCurl Space]")
	fes_hcurl = HCurl(mesh, order=1)
	gf_hcurl_set = GridFunction(fes_hcurl)
	gf_hcurl_interp = GridFunction(fes_hcurl)

	# Method 1: Set()
	print(f"    Set() method:")
	t_start = time.perf_counter()
	gf_hcurl_set.Set(B_cf)
	t_set = time.perf_counter() - t_start
	print(f"      Time: {t_set*1000:.2f} ms ({t_set*1e6/mesh.nv:.2f} us/vertex)")

	# Method 2: Interpolate() - may not be available for HCurl
	try:
		print(f"    Interpolate() method:")
		t_start = time.perf_counter()
		gf_hcurl_interp.Interpolate(B_cf)
		t_interp = time.perf_counter() - t_start
		print(f"      Time: {t_interp*1000:.2f} ms ({t_interp*1e6/mesh.nv:.2f} us/vertex)")
		print(f"      Speedup: {t_set/t_interp:.2f}x")
	except Exception as e:
		print(f"      Not available for HCurl: {e}")

	# Test with H1 space (if possible)
	print("\n  [H1^3 Space (component-wise)]")
	try:
		fes_h1 = VectorH1(mesh, order=1)
		gf_h1_set = GridFunction(fes_h1)
		gf_h1_interp = GridFunction(fes_h1)

		# Method 1: Set()
		print(f"    Set() method:")
		t_start = time.perf_counter()
		gf_h1_set.Set(B_cf)
		t_set_h1 = time.perf_counter() - t_start
		print(f"      Time: {t_set_h1*1000:.2f} ms ({t_set_h1*1e6/mesh.nv:.2f} us/vertex)")

		# Method 2: Interpolate()
		print(f"    Interpolate() method:")
		t_start = time.perf_counter()
		gf_h1_interp.Interpolate(B_cf)
		t_interp_h1 = time.perf_counter() - t_start
		print(f"      Time: {t_interp_h1*1000:.2f} ms ({t_interp_h1*1e6/mesh.nv:.2f} us/vertex)")
		print(f"      Speedup: {t_set_h1/t_interp_h1:.2f}x")

		# Compare values at test point
		test_point = mesh(0.01, 0.0, 0.0)
		val_set = gf_h1_set(test_point)
		val_interp = gf_h1_interp(test_point)
		error = np.linalg.norm(np.array(val_set) - np.array(val_interp))
		print(f"\n    Value comparison at (0.01, 0, 0):")
		print(f"      Set():        [{val_set[0]:.6e}, {val_set[1]:.6e}, {val_set[2]:.6e}]")
		print(f"      Interpolate(): [{val_interp[0]:.6e}, {val_interp[1]:.6e}, {val_interp[2]:.6e}]")
		print(f"      Difference: {error:.6e}")

	except Exception as e:
		print(f"    H1 space test failed: {e}")

print("\n" + "=" * 80)
print("CONCLUSION")
print("=" * 80)

print("""
If Interpolate() shows significant speedup:
  - Interpolate() evaluates at vertices (fewer points than integration points)
  - This allows better utilization of batch evaluation
  - However, Interpolate() may give different results than Set() for HCurl

Key observations:
  - Set() evaluates at integration points (many per element)
  - Interpolate() evaluates at interpolation nodes (typically vertices)
  - For H1 spaces, Interpolate() is often faster
  - For HCurl/HDiv, Interpolate() may not be appropriate
""")

print("=" * 80)
