#!/usr/bin/env python3
"""
Test batch evaluation in radia_ngsolve.RadiaField

Tests whether the new batch evaluation method improves GridFunction.Set() performance.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "build", "Release"))

import radia as rad
try:
	from ngsolve import *
	from netgen.occ import *
	import radia_ngsolve
	NGSOLVE_AVAILABLE = True
except ImportError:
	print("ERROR: NGSolve not available")
	sys.exit(1)

import time
import numpy as np

print("=" * 80)
print("Batch Evaluation Test for radia_ngsolve.RadiaField")
print("=" * 80)

# Create magnet
rad.UtiDelAll()
rad.FldUnits('m')  # Set units to meters
n = 5  # 5x5x5 = 125 elements
cube_size = 0.100  # meters
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
box_size_m = cube_size * 2.0  # meters
geo = Box(
	Pnt(-box_size_m/2, -box_size_m/2, -box_size_m/2),
	Pnt( box_size_m/2,  box_size_m/2,  box_size_m/2)
)

# Test different mesh sizes
mesh_configs = [
	{"h": 0.05, "desc": "Coarse (h=0.05m)"},
	{"h": 0.025, "desc": "Medium (h=0.025m)"},
]

print("\n" + "=" * 80)
print("PERFORMANCE TEST: GridFunction.Set() with Batch Evaluation")
print("=" * 80)

for config in mesh_configs:
	h = config['h']
	desc = config['desc']

	print(f"\n{desc}")
	print("-" * 80)

	mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=h))
	print(f"  Mesh: {mesh.nv} vertices, {mesh.ne} elements")

	# Create CoefficientFunction
	B_cf = radia_ngsolve.RadiaField(magnet, 'b')

	# Create GridFunction
	fes = HCurl(mesh, order=1)
	gf = GridFunction(fes)

	# Measure time
	print(f"  Testing GridFunction.Set()...")
	t_start = time.perf_counter()
	gf.Set(B_cf)
	t_set = time.perf_counter() - t_start

	print(f"  Time: {t_set*1000:.2f} ms ({t_set*1e6/mesh.nv:.2f} us/vertex)")

	# Verify field value at a test point
	test_point = mesh(0.01, 0.0, 0.0)
	B_val = B_cf(test_point)
	print(f"  Sample field at (0.01, 0, 0) m: B = [{B_val[0]:.6f}, {B_val[1]:.6f}, {B_val[2]:.6f}] T")

print("\n" + "=" * 80)
print("Test all field types")
print("=" * 80)

mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=0.05))
field_types = [
	('b', 'Magnetic flux density B'),
	('h', 'Magnetic field H'),
	('a', 'Vector potential A'),
	('m', 'Magnetization M'),
]

for field_type, field_name in field_types:
	print(f"\n{field_name} ('{field_type}')")
	print("-" * 80)

	cf = radia_ngsolve.RadiaField(magnet, field_type)
	gf = GridFunction(HCurl(mesh))

	t_start = time.perf_counter()
	gf.Set(cf)
	t_set = time.perf_counter() - t_start

	print(f"  Time: {t_set*1000:.2f} ms")

	# Sample value
	test_point = mesh(0.01, 0.0, 0.0)
	val = cf(test_point)
	print(f"  Sample at (0.01, 0, 0) m: [{val[0]:.6e}, {val[1]:.6e}, {val[2]:.6e}]")

print("\n" + "=" * 80)
print("CONCLUSION")
print("=" * 80)

print("""
If batch evaluation is working:
  - GridFunction.Set() should be faster than before
  - Time per vertex should be reduced significantly

Compare with previous benchmark results:
  - Old implementation: ~400-500 ms for 826 vertices (125 elements)
  - Expected with batch: ~60-100 ms (6-8x speedup)
""")

print("=" * 80)
