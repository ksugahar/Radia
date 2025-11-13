#!/usr/bin/env python
"""
Step 2: Load H-matrix from disk (instant startup)
This simulates restarting the program - geometry recreated, H-matrix loaded from disk
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))

import radia as rad
import time

print("="*70)
print("STEP 2: Load H-matrix from disk (instant startup)")
print("="*70)

# Re-create EXACT SAME geometry (16x16 = 256 elements)
print("\nRe-creating 16x16 = 256 element geometry...")
print("(Same as step 1 - should have matching geometry hash)")
g = rad.ObjCnt([])
for i in range(16):
    for j in range(16):
        x, y = i * 5.0, j * 5.0
        elem = rad.ObjRecMag([x, y, 0], [4, 4, 4], [0, 0, 100])
        mat = rad.MatSatIsoFrm([20000, 2], [0.1, 0.1], [0.1, 0.1])
        rad.MatApl(elem, mat)
        rad.ObjAddToCnt(g, [elem])

print("Created 256 elements with nonlinear material")

# Enable full serialization
print("\nEnabling full H-matrix disk serialization...")
rad.SolverHMatrixCacheFull(1)
rad.SolverHMatrixEnable(1, 1e-4, 30)

# Load H-matrix from disk (should be instant)
print("\nLoading H-matrix from disk (should be instant)...")
t_start = time.time()
rad.RlxPre(g, 1)
t_load = time.time() - t_start

print(f"\n[OK] H-matrix loaded in {t_load:.6f}s")
print("\nIf this is much faster than step 1, disk serialization is working!")
