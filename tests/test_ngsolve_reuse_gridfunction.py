#!/usr/bin/env python
"""
Test if reusing GridFunction reduces memory leak in NGSolve integration
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

import gc
import tracemalloc
import numpy as np
import radia as rad
from ngsolve import *
import rad_ngsolve

print("=" * 80)
print("NGSolve GridFunction Reuse Test")
print("=" * 80)
print()

rad.FldUnits('m')

# Create simple mesh
from netgen.occ import *
box = Box((0,0,0), (0.1,0.1,0.1))
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=0.02))
print(f"Mesh: {mesh.nv} vertices, {mesh.ne} elements")

# Create FE space
fes = HCurl(mesh, order=2)
print(f"FE Space: {fes.ndof} DOFs")
print()

# Create REUSABLE GridFunction ONCE
gf_B = GridFunction(fes)
print("[Setup] Created GridFunction once (will reuse)")
print()

tracemalloc.start()
gc.collect()
mem_start = tracemalloc.get_traced_memory()[0]

NUM_STEPS = 100

for step in range(NUM_STEPS):
    x_pos = step * 0.001
    new_pos = [x_pos, 0.0, 0.0]

    # Update magnet position
    rad.UtiDelAll()
    magnet = rad.ObjRecMag(new_pos, [0.01, 0.01, 0.01], [0, 0, 1.0])

    # Create NEW CoefficientFunction each time
    B_cf = rad_ngsolve.RadiaField(magnet, 'b')

    # REUSE GridFunction - just call Set() again
    gf_B.Set(B_cf)

    # Test evaluation
    test_point = mesh(0.02, 0.0, 0.0)
    B_val = gf_B(test_point)

    # Clean up CF (but keep gf_B)
    del B_cf

    if step % 10 == 0:
        gc.collect()
        mem_current = tracemalloc.get_traced_memory()[0]
        mem_mb = mem_current / 1024 / 1024
        print(f"Step {step:3d}: Memory={mem_mb:.2f} MB")

gc.collect()
mem_end = tracemalloc.get_traced_memory()[0]
growth = (mem_end - mem_start) / 1024 / 1024

print()
print("=" * 80)
print(f"Memory growth: {growth:.2f} MB for {NUM_STEPS} steps")
print(f"Per step: {growth / NUM_STEPS * 1024:.2f} KB/step")
print("=" * 80)

if growth < 50:
    print("\n[PASS] Memory leak significantly reduced!")
else:
    print(f"\n[FAIL] Still leaking {growth:.2f} MB")

tracemalloc.stop()
rad.UtiDelAll()
