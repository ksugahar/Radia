#!/usr/bin/env python
"""Compare Method 4, 9, and 10 results."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))
import radia as rad
import time

B0 = 1.0
def uniform_B(pos):
    return [0, 0, B0]

def run_test(method, max_iter):
    rad.UtiDelAll()
    rad.FldUnits('m')

    cube = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 0])
    rad.ObjDivMag(cube, [4, 4, 4])

    chi = 99  # mu_r = 100
    mat = rad.MatLin(chi)
    rad.MatApl(cube, mat)

    bg = rad.ObjBckgCF(uniform_B)
    system = rad.ObjCnt([cube, bg])

    t0 = time.time()
    result = rad.Solve(system, 0.001, max_iter, method)
    elapsed = time.time() - t0

    B = rad.Fld(cube, 'b', [0.1, 0, 0])
    return result, elapsed, B

print("Method Comparison Test")
print("=" * 60)

# Method 4
result4, t4, B4 = run_test(4, 20000)
print(f"Method 4: iter={int(result4[3])}, time={t4:.4f}s, Bz={B4[2]:.6f}")

# Method 9
result9, t9, B9 = run_test(9, 1)
print(f"Method 9: iter={int(result9[3])}, time={t9:.4f}s, Bz={B9[2]:.6f}")

# Method 10
result10, t10, B10 = run_test(10, 100)
print(f"Method 10: iter={int(result10[3])}, time={t10:.4f}s, Bz={B10[2]:.6f}")

print()
print("Bz Differences:")
print(f"  M4 vs M9:  {abs(B4[2] - B9[2]):.2e}")
print(f"  M4 vs M10: {abs(B4[2] - B10[2]):.2e}")
print(f"  M9 vs M10: {abs(B9[2] - B10[2]):.2e}")
