"""Minimal test to debug segfault."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

print("Step 1: Importing radia...")
import radia as rad
print("  OK")

print("Step 2: rad.UtiDelAll()...")
rad.UtiDelAll()
print("  OK")

print("Step 3: rad.FldUnits('mm')...")
rad.FldUnits('mm')
print("  OK")

print("Step 4: Creating simple hex...")
vertices = [
    [0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0],
    [0, 0, 10], [10, 0, 10], [10, 10, 10], [0, 10, 10]
]
hex_obj = rad.ObjHexahedron(vertices, [0, 0, 0])
print(f"  OK: {hex_obj}")

print("Step 5: Loading BH curve...")
import numpy as np
bh_data = np.loadtxt("BH.txt", comments='#').tolist()
print(f"  OK: {len(bh_data)} points")

print("Step 6: Creating material...")
mat = rad.MatSatIsoTab(bh_data)
print(f"  OK: {mat}")

print("Step 7: Applying material...")
rad.MatApl(hex_obj, mat)
print("  OK")

print("Step 8: Creating container...")
yoke = rad.ObjCnt([hex_obj])
print(f"  OK: {yoke}")

print("Step 9: Solving...")
result = rad.Solve(yoke, 0.0001, 100, 0)
print(f"  OK: {result}")

print("\nAll tests passed!")
