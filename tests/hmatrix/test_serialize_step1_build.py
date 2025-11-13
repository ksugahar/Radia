#!/usr/bin/env python
"""
Step 1: Build H-matrix and save to disk
This simulates the first time a user runs their simulation
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))

import radia as rad
import time

print("="*70)
print("STEP 1: Build H-matrix and save to disk")
print("="*70)

# Create geometry (16x16 = 256 elements)
print("\nCreating 16x16 = 256 element geometry...")
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

# Build H-matrix (will save to disk)
print("\nBuilding H-matrix (first time, will save to disk)...")
t_start = time.time()
rad.RlxPre(g, 1)
t_build = time.time() - t_start

print(f"\n[OK] H-matrix built and saved to disk in {t_build:.3f}s")
print("\nNext: Run test_serialize_step2_load.py to test loading from disk")
