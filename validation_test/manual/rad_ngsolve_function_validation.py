import sys, os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "src", "radia"))

from numpy import *
from ngsolve import *
from ngsolve import TaskManager
from netgen.occ import *
import radia as rad


# ========================================================================
# Radia model
# ========================================================================

magnet_size = 1.0
magnetization = 1.0

magnet_base = rad.ObjRecMag( [0, 2, 0], [magnet_size, magnet_size, magnet_size], [0, magnetization, 0])

test_points = [
    ([0, 0, 0], "origin"),
    ([0, 1, 0], "1 mm from magnet in +Y"),
    ([0, 2, 0], "magnet center (0, 2, 0)"),
    ([0, 3, 0], "1 mm from magnet in +Y"),
    ([1, 0, 0], "1 mm from magnet in +X"),
    ([0, 0, 1], "1 mm from magnet in +Z"),
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
B_cf = rad.RadiaField(magnet_base, 'b')
gf_B = GridFunction(fes)
with TaskManager():
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
