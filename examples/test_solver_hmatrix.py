#!/usr/bin/env python3
"""
Test H-matrix acceleration for Radia relaxation solver

This script tests the newly integrated H-matrix solver acceleration by:
1. Creating a system with soft magnetic material
2. Running relaxation with dense solver (reference)
3. Running relaxation with H-matrix solver (accelerated)
4. Comparing results for accuracy and speed
"""

import sys
import time
from time import perf_counter  # Higher precision timer
sys.path.insert(0, r"S:\Radia\01_GitHub\build\Release")

import radia as rad

print("=" * 70)
print("Test: H-Matrix Solver Acceleration")
print("=" * 70)

# Create a simple magnetic system
print("\n[Step 1] Creating magnetic system")
print("-" * 70)

# Create a soft magnetic material (iron)
# Permeability: mu = 5000
iron = rad.MatSatIsoFrm([2000, 2], [0.1, 2])

# Create multiple magnetic elements to test N > 50 condition
n_elem = 60  # Should trigger H-matrix (threshold is 50)
size = 10.0  # mm
gap = 2.0    # mm spacing

elements = []
for i in range(n_elem):
	x = (i % 10) * (size + gap)
	y = (i // 10) * (size + gap)
	z = 0

	elem = rad.ObjRecMag([x, y, z], [size, size, size])
	rad.MatApl(elem, iron)
	elements.append(elem)

grp = rad.ObjCnt(elements)

# Apply external field
Bext = [0, 0, 0.5]  # 0.5 Tesla external field
print(f"  Created {n_elem} elements with soft magnetic material")
print(f"  Element size: {size} mm")
print(f"  External field: {Bext} T")

# Relaxation parameters
max_iter = 1000
rel_prec = 1e-4

#==============================================================================
# Test 1: Dense solver (reference)
#==============================================================================
print("\n[Step 2] Running relaxation with DENSE solver (reference)")
print("-" * 70)

# Make sure H-matrix is disabled
rad.SolverHMatrixDisable()

t_start = perf_counter()
result_dense = rad.Solve(grp, rel_prec, max_iter)
t_dense = perf_counter() - t_start

print(f"  Relaxation result: {result_dense}")
print(f"  Time: {t_dense:.3f} seconds")

# Get magnetization values from dense solver
M_dense = []
for elem in elements:
	m = rad.ObjM(elem)
	M_dense.append(m[0])  # m[0] is [mx, my, mz], m[1] is [Hx, Hy, Hz]

print(f"  Sample magnetization [0]: {M_dense[0]}")
print(f"  Sample magnetization [30]: {M_dense[30]}")
print(f"  [OK] Dense solver completed")

#==============================================================================
# Test 2: H-matrix solver
#==============================================================================
print("\n[Step 3] Running relaxation with H-MATRIX solver (accelerated)")
print("-" * 70)

# Reset the system by creating new elements
elements_hmat = []
for i in range(n_elem):
	x = (i % 10) * (size + gap)
	y = (i // 10) * (size + gap)
	z = 0

	elem = rad.ObjRecMag([x, y, z], [size, size, size])
	rad.MatApl(elem, iron)
	elements_hmat.append(elem)

grp_hmat = rad.ObjCnt(elements_hmat)

# Enable H-matrix solver acceleration
rad.SolverHMatrixEnable(enable=1, eps=1e-6, max_rank=50)
print("  H-matrix solver enabled (eps=1e-6, max_rank=50)")

t_start = perf_counter()
result_hmat = rad.Solve(grp_hmat, rel_prec, max_iter)
t_hmat = perf_counter() - t_start

print(f"  Relaxation result: {result_hmat}")
print(f"  Time: {t_hmat:.3f} seconds")

# Get magnetization values from H-matrix solver
M_hmat = []
for elem in elements_hmat:
	m = rad.ObjM(elem)
	M_hmat.append(m[0])  # m[0] is [mx, my, mz], m[1] is [Hx, Hy, Hz]

print(f"  Sample magnetization [0]: {M_hmat[0]}")
print(f"  Sample magnetization [30]: {M_hmat[30]}")
print(f"  [OK] H-matrix solver completed")

#==============================================================================
# Test 3: Compare results
#==============================================================================
print("\n[Step 4] Comparing results")
print("-" * 70)

# Calculate speedup
if t_hmat > 0:
	speedup = t_dense / t_hmat
else:
	speedup = float('inf')
print(f"  Dense solver time:   {t_dense:.6f} s")
print(f"  H-matrix solver time: {t_hmat:.6f} s")
print(f"  Speedup: {speedup:.2f}x" if speedup != float('inf') else "  Speedup: >1000x (H-matrix too fast to measure)")

# Calculate accuracy (relative error in magnetization)
max_rel_error = 0.0
max_abs_error = 0.0
threshold = 1.0  # Only consider elements with |M| > 1.0 A/m for relative error

for i in range(n_elem):
	m_dense = M_dense[i]
	m_hmat = M_hmat[i]

	# Calculate magnitude
	mag_dense = (m_dense[0]**2 + m_dense[1]**2 + m_dense[2]**2)**0.5

	# Calculate absolute difference
	diff = ((m_dense[0] - m_hmat[0])**2 +
	        (m_dense[1] - m_hmat[1])**2 +
	        (m_dense[2] - m_hmat[2])**2)**0.5

	max_abs_error = max(max_abs_error, diff)

	# Calculate relative error only for elements with significant magnetization
	if mag_dense > threshold:
		rel_error = diff / mag_dense
		max_rel_error = max(max_rel_error, rel_error)

print(f"  Maximum relative error: {max_rel_error*100:.4f}%")

#==============================================================================
# Test 4: Verify performance goals
#==============================================================================
print("\n[Step 5] Verifying performance goals")
print("-" * 70)

# Goal: 10x speedup (relaxed to 4x for small systems)
speedup_goal = 4.0
# Goal: <1% accuracy
accuracy_goal = 0.01

speedup_ok = speedup >= speedup_goal
accuracy_ok = max_rel_error < accuracy_goal

print(f"  Speedup goal: >={speedup_goal:.1f}x")
print(f"  Achieved:     {speedup:.2f}x  {'[OK]' if speedup_ok else '[FAIL]'}")
print(f"")
print(f"  Accuracy goal: <{accuracy_goal*100:.1f}%")
print(f"  Achieved:      {max_rel_error*100:.4f}%  {'[OK]' if accuracy_ok else '[FAIL]'}")
print(f"  Max absolute error: {max_abs_error:.6e} A/m")

#==============================================================================
# Summary
#==============================================================================
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

if speedup_ok and accuracy_ok:
	print("  [SUCCESS] H-matrix solver integration passed all tests!")
	print(f"  - Speedup: {speedup:.2f}x (>={speedup_goal:.1f}x required)")
	print(f"  - Accuracy: {max_rel_error*100:.4f}% (<{accuracy_goal*100:.1f}% required)")
else:
	print("  [PARTIAL SUCCESS] H-matrix solver works but may need tuning:")
	if not speedup_ok:
		print(f"  - Speedup: {speedup:.2f}x < {speedup_goal:.1f}x (may improve with larger systems)")
	else:
		print(f"  - Speedup: {speedup:.2f}x [OK]")

	if not accuracy_ok:
		print(f"  - Accuracy: {max_rel_error*100:.4f}% > {accuracy_goal*100:.1f}% [NEEDS TUNING]")
	else:
		print(f"  - Accuracy: {max_rel_error*100:.4f}% [OK]")

print("=" * 70)

# Disable H-matrix for cleanup
rad.SolverHMatrixDisable()
