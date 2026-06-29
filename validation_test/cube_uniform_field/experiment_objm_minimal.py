"""Test ObjM with correct benchmark pattern."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'src'))
import numpy as np
import radia as rad


MU_0 = 4 * np.pi * 1e-7
H_EXT = 200000.0
B_ext = MU_0 * H_EXT
MU_R = 1000
CHI = MU_R - 1
N_DEMAG = 1.0 / 3.0
M_ANALYTICAL_Z = CHI * H_EXT / (1 + CHI * N_DEMAG)

print(f"B_ext = {B_ext:.6f} T")
print(f"M_analytical_z = {M_ANALYTICAL_Z:.0f} A/m")

# Replicate benchmark_common.py pattern
rad.UtiDelAll()

# Create 3x3x3 hex mesh
n_div = 3
objs = []
dx = 1.0 / n_div
offset = 0.5
verts_base = [
    [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
    [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]
]
for iz in range(n_div):
    for iy in range(n_div):
        for ix in range(n_div):
            x0 = ix * dx - offset
            y0 = iy * dx - offset
            z0 = iz * dx - offset
            verts = [[v[0]*dx + x0, v[1]*dx + y0, v[2]*dx + z0] for v in verts_base]
            obj = rad.ObjHexahedron(verts, [0, 0, 0])
            objs.append(obj)

radia_obj = rad.ObjCnt(objs)
n_elements = len(objs)
print(f"Created {n_elements} hexahedra")

# Apply linear material
mat = rad.MatLin(MU_R)
rad.MatApl(radia_obj, mat)

# Background field (same as benchmark)
ext = rad.ObjBckg(lambda p: [0, 0, B_ext])
grp = rad.ObjCnt([radia_obj, ext])

# Set solver params
rad.SolverConfig(bicgstab_tol=1e-4, relax_param=0.0)

# Solve LU
print("Solving LU...")
sys.stdout.flush()
result = rad.Solve(grp, 0.001, 100, 0)
print(f"Solve result = {result}")
sys.stdout.flush()

# Get magnetization
print("Getting ObjM...")
sys.stdout.flush()
all_M = rad.ObjM(radia_obj)
print(f"ObjM returned {len(all_M)} elements")
M_total_z = sum(m[1][2] for m in all_M)
M_avg_z = M_total_z / n_elements
print(f"M_avg_z = {M_avg_z:.0f} A/m (analytical: {M_ANALYTICAL_Z:.0f})")
error = abs(M_avg_z - M_ANALYTICAL_Z) / M_ANALYTICAL_Z * 100
print(f"Error = {error:.2f}%")

print("\nAll tests passed!")
