"""Test HACApK MatVec timing"""
import sys
import os
_src_path = os.path.join(os.path.dirname(__file__), '../../../src/radia')
sys.path.insert(0, _src_path)
import radia as rad
import time
import numpy as np

# Use meters
rad.FldUnits('m')

# Parameters
N = 10
cube_size = 1.0
H_ext = 200000.0  # A/m
eps = 1e-4

print(f"Setting up N={N} cube ({N**3} elements)...")

# Create hexahedral mesh
cell_size = cube_size / N
half_size = cube_size / 2.0

# Create container and elements
container = rad.ObjCnt([])
for iz in range(N):
    for iy in range(N):
        for ix in range(N):
            cx = -half_size + (ix + 0.5) * cell_size
            cy = -half_size + (iy + 0.5) * cell_size
            cz = -half_size + (iz + 0.5) * cell_size
            elem = rad.ObjRecMag([cx, cy, cz], [cell_size]*3, [0, 0, 0])
            rad.ObjAddToCnt(container, [elem])

# Apply nonlinear material
steel_bh = [
    [0.0, 0.0],
    [100.0, 0.5],
    [500.0, 1.2],
    [1000.0, 1.5],
    [5000.0, 1.8],
    [10000.0, 1.95],
    [50000.0, 2.1],
    [100000.0, 2.2],
    [200000.0, 2.3],
]
mat = rad.MatSatIsoTab(steel_bh)
rad.MatApl(container, mat)

# Apply external field
rad.ObjDrwAtr(container, [0, 0, 0], 0.001)
rad.FldExtSrc(container, [0, 0, H_ext])

# Set HACApK parameters
print(f"Setting HACApK eps={eps}")
rad.SetHMatrixEpsilon(eps)

# Solve with HACApK
print("Solving with HACApK...")
t0 = time.time()
rad.Solve(container, 0.0001, 1000, 2)  # Method 2 = HACApK
t_solve = time.time() - t0

print(f"Total solve time: {t_solve:.3f} s")

# Get statistics
info = rad.UtiMag(container)
print(f"M_avg_z: {info[2]/1000:.0f} kA/m")
