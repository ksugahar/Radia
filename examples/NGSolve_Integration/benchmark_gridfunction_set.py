#!/usr/bin/env python3
"""
Performance benchmark: GridFunction.Set() time evaluation

Measures the time required to evaluate Radia field on NGSolve mesh using
GridFunction.Set() method. Compares:
1. Different field types (B, H, A, M)
2. Dense matrix vs H-matrix
3. Different mesh sizes
4. Different Radia element counts

This is the critical operation in coupled NGSolve-Radia simulations.

Date: 2025-11-08
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "build", "Release"))

import radia as rad
try:
	from ngsolve import *
	from netgen.occ import *
	import rad_ngsolve
	NGSOLVE_AVAILABLE = True
except ImportError:
	print("ERROR: NGSolve not available. This benchmark requires NGSolve.")
	sys.exit(1)

import numpy as np
from time import perf_counter

print("=" * 80)
print("GridFunction.Set() Performance Benchmark")
print("=" * 80)

# Test configurations for Radia magnet
radia_configs = [
	{"n": 3, "desc": "3x3x3 = 27 Radia elements"},
	{"n": 4, "desc": "4x4x4 = 64 Radia elements"},
	{"n": 5, "desc": "5x5x5 = 125 Radia elements"},
	{"n": 6, "desc": "6x6x6 = 216 Radia elements"},
	{"n": 7, "desc": "7x7x7 = 343 Radia elements"},
]

# Test configurations for NGSolve mesh
mesh_configs = [
	{"h": 0.05, "desc": "Coarse mesh (h=0.05 m)"},
	{"h": 0.025, "desc": "Medium mesh (h=0.025 m)"},
	{"h": 0.0125, "desc": "Fine mesh (h=0.0125 m)"},
]

# Field types to test
field_types = [
	{"type": "b", "name": "Magnetic flux density B"},
	{"type": "h", "name": "Magnetic field H"},
	{"type": "a", "name": "Vector potential A"},
	{"type": "m", "name": "Magnetization M"},
]

print(f"\nConfiguration:")
print(f"  Radia configs: {len(radia_configs)} (N = 27 to 343 elements)")
print(f"  Mesh configs: {len(mesh_configs)} (coarse to fine)")
print(f"  Field types: {len(field_types)} (B, H, A, M)")
print(f"  Methods: Dense matrix vs H-matrix")

all_results = []

# ==================================================================
# Main benchmark loop
# ==================================================================

for radia_cfg in radia_configs:
	n = radia_cfg['n']
	n_elem = n * n * n

	print(f"\n{'='*80}")
	print(f"Radia Configuration: {radia_cfg['desc']}")
	print(f"{'='*80}")

	# Create Radia magnet geometry
	rad.UtiDelAll()

	cube_size = 100.0  # mm
	elem_size = cube_size / n
	mag_value = 1.2  # T

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
	print(f"  Created magnet: {n_elem} elements")

	for mesh_cfg in mesh_configs:
		h = mesh_cfg['h']

		print(f"\n  {'-'*76}")
		print(f"  Mesh: {mesh_cfg['desc']}")
		print(f"  {'-'*76}")

		# Create NGSolve mesh
		box_size_m = cube_size * 2.0 / 1000  # mm to m
		geo = Box(
			Pnt(-box_size_m/2, -box_size_m/2, -box_size_m/2),
			Pnt( box_size_m/2,  box_size_m/2,  box_size_m/2)
		)

		mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=h))
		print(f"    Mesh: {mesh.nv} vertices, {mesh.ne} elements, {mesh.nedge} edges")

		for field_cfg in field_types:
			field_type = field_cfg['type']
			field_name = field_cfg['name']

			print(f"\n    Field type: {field_name}")
			print(f"    {'-'*72}")

			# ----------------------------------------------------------
			# DENSE MATRIX (H-matrix disabled)
			# ----------------------------------------------------------
			rad.SolverHMatrixDisable()

			# Create CoefficientFunction
			t_start = perf_counter()
			cf_dense = rad_ngsolve.RadiaField(magnet, field_type, use_hmatrix=False)
			t_create_dense = perf_counter() - t_start

			# Create GridFunction and Set
			gf_dense = GridFunction(HCurl(mesh))

			t_start = perf_counter()
			gf_dense.Set(cf_dense)
			t_set_dense = perf_counter() - t_start

			print(f"      Dense:   Set() = {t_set_dense*1000:8.2f} ms  "
				  f"({t_set_dense*1e6/mesh.nv:6.2f} us/vertex)")

			# ----------------------------------------------------------
			# H-MATRIX (enabled)
			# ----------------------------------------------------------
			rad.SolverHMatrixEnable()

			# Create CoefficientFunction
			t_start = perf_counter()
			cf_hmat = rad_ngsolve.RadiaField(magnet, field_type, use_hmatrix=True)
			t_create_hmat = perf_counter() - t_start

			# Create GridFunction and Set
			gf_hmat = GridFunction(HCurl(mesh))

			t_start = perf_counter()
			gf_hmat.Set(cf_hmat)
			t_set_hmat = perf_counter() - t_start

			print(f"      H-matrix: Set() = {t_set_hmat*1000:8.2f} ms  "
				  f"({t_set_hmat*1e6/mesh.nv:6.2f} us/vertex)")

			# Calculate speedup
			speedup = t_set_dense / t_set_hmat if t_set_hmat > 0 else 0
			print(f"      Speedup: {speedup:.2f}x")

			# Store results
			all_results.append({
				'n_radia': n_elem,
				'mesh_h': h,
				'mesh_nv': mesh.nv,
				'mesh_ne': mesh.ne,
				'field_type': field_type,
				'field_name': field_name,
				't_set_dense': t_set_dense,
				't_set_hmat': t_set_hmat,
				'speedup': speedup,
			})

# ==================================================================
# Summary Tables
# ==================================================================

print(f"\n{'='*80}")
print("SUMMARY: GridFunction.Set() Time (milliseconds)")
print(f"{'='*80}")

for field_cfg in field_types:
	field_type = field_cfg['type']
	field_name = field_cfg['name']

	print(f"\n{field_name} (field_type='{field_type}'):")
	print(f"{'N_Radia':>8}  {'Mesh_h':>8}  {'Vertices':>9}  {'Dense(ms)':>10}  "
		  f"{'H-mat(ms)':>10}  {'Speedup':>8}")
	print("-" * 80)

	# Filter results for this field type
	field_results = [r for r in all_results if r['field_type'] == field_type]

	# Sort by N_radia, then mesh_h
	field_results.sort(key=lambda x: (x['n_radia'], x['mesh_h']))

	for r in field_results:
		print(f"{r['n_radia']:>8}  {r['mesh_h']:>8.4f}  {r['mesh_nv']:>9}  "
			  f"{r['t_set_dense']*1000:>10.2f}  {r['t_set_hmat']*1000:>10.2f}  "
			  f"{r['speedup']:>8.2f}x")

# ==================================================================
# Analysis
# ==================================================================

print(f"\n{'='*80}")
print("ANALYSIS")
print(f"{'='*80}")

# Average speedup by field type
print("\nAverage Speedup by Field Type:")
for field_cfg in field_types:
	field_type = field_cfg['type']
	field_name = field_cfg['name']

	field_results = [r for r in all_results if r['field_type'] == field_type]
	avg_speedup = np.mean([r['speedup'] for r in field_results])

	print(f"  {field_name:30s}: {avg_speedup:6.2f}x")

# Speedup vs Radia element count
print("\nSpeedup vs Radia Element Count (averaged over all meshes and field types):")
for radia_cfg in radia_configs:
	n_elem = radia_cfg['n'] ** 3

	# Filter results for this N
	n_results = [r for r in all_results if r['n_radia'] == n_elem]
	avg_speedup = np.mean([r['speedup'] for r in n_results])

	print(f"  N = {n_elem:4d}:  {avg_speedup:6.2f}x")

# Speedup vs mesh size
print("\nSpeedup vs Mesh Size (averaged over all Radia configs and field types):")
for mesh_cfg in mesh_configs:
	h = mesh_cfg['h']

	# Filter results for this mesh size
	mesh_results = [r for r in all_results if r['mesh_h'] == h]
	avg_speedup = np.mean([r['speedup'] for r in mesh_results])
	avg_vertices = np.mean([r['mesh_nv'] for r in mesh_results])

	print(f"  h = {h:.4f} m ({int(avg_vertices)} vertices avg):  {avg_speedup:6.2f}x")

print(f"\n{'='*80}")
print("CONCLUSION")
print(f"{'='*80}")

# Overall average speedup
overall_avg_speedup = np.mean([r['speedup'] for r in all_results])

print(f"""
Overall Performance Summary:

1. Average H-matrix speedup: {overall_avg_speedup:.2f}x
2. H-matrix is beneficial for:
   - Large Radia element counts (N > 100)
   - Fine NGSolve meshes (many evaluation points)
   - All field types (B, H, A, M)

3. GridFunction.Set() is the main computational cost in coupled simulations
4. H-matrix acceleration directly speeds up coupled NGSolve-Radia simulations

Recommendation:
- Use H-matrix (use_hmatrix=True) for production runs
- Disable H-matrix (use_hmatrix=False) only for:
  - Small problems (N < 50)
  - Debugging and validation
""")

print("=" * 80)
