import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/python'))

from numpy import *
from ngsolve import *
from netgen.occ import *
import radia as rad
import radia_ngsolve


# ========================================================================
# radiaモデル
# ========================================================================

magnet_size = 1.0
magnetization = 1.0

magnet_base = rad.ObjRecMag( [0, 2, 0], [magnet_size, magnet_size, magnet_size], [0, magnetization, 0])

test_points = [
    ([0, 0, 0], "原点"),
    ([0, 1, 0], "磁石から+Y方向1mm"),
    ([0, 2, 0], "磁石中心 (0, 2, 0)"),
    ([0, 3, 0], "磁石から+Y方向1mm"),
    ([1, 0, 0], "磁石から+X方向1mm"),
    ([0, 0, 1], "磁石から+Z方向1mm"),
]

print("="*60)
print("rad.Fld Results (Radia Direct)")
print("="*60)
for point, description in test_points:
	B_radia = rad.Fld(magnet_base, 'b', point)
	print(f"{description}: {B_radia}")

mesh_domain = 6.0e-3
air_region = Box((-mesh_domain, -mesh_domain, -mesh_domain), (mesh_domain, mesh_domain, mesh_domain)).mat("air")
mesh_maxh = 1.0e-3
mesh = air_region.GenerateMesh(maxh=mesh_maxh)

# FIXED: Remove dim=3 - it was creating a 3x3 tensor field instead of 3D vector field
fes = VectorH1(mesh, order=2)
B_cf = radia_ngsolve.RadiaField(magnet_base, 'b')
gf_B = GridFunction(fes)
gf_B.Set(B_cf)

print("\n" + "="*60)
print("rad_ngsolve Results (NGSolve via GridFunction)")
print("="*60)
for point, description in test_points:
	# Convert mm to m
	point_m = (point[0]/1000, point[1]/1000, point[2]/1000)
	B_ngsolve = gf_B(mesh(*point_m))
	print(f"{description}: {B_ngsolve}")

print("\n" + "="*60)
print("Comparison: rad.Fld vs rad_ngsolve")
print("="*60)
print(f"{'Point':<30s} {'rad.Fld By':>15s} {'NGSolve By':>15s} {'Error %':>12s}")
print("-"*60)
for point, description in test_points:
	B_radia = rad.Fld(magnet_base, 'b', point)
	point_m = (point[0]/1000, point[1]/1000, point[2]/1000)
	B_ngsolve = gf_B(mesh(*point_m))
	if abs(B_radia[1]) > 1e-6:
		rel_error = abs(B_radia[1] - B_ngsolve[1]) / abs(B_radia[1]) * 100
	else:
		rel_error = abs(B_radia[1] - B_ngsolve[1]) * 100
	print(f"{description:<30s} {B_radia[1]:15.6f} {B_ngsolve[1]:15.6f} {rel_error:11.2f}%")
