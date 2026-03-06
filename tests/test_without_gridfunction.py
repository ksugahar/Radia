import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/python'))

from numpy import *
from ngsolve import *
from netgen.occ import *
import radia as rad
import radia_ngsolve

# Create magnet
magnet_size = 1.0  # mm
magnetization = 1.0  # T
magnet_base = rad.ObjRecMag([0, 2, 0], [magnet_size, magnet_size, magnet_size], [0, magnetization, 0])

test_points = [
    ([0, 0, 0], "Origin"),
    ([0, 1, 0], "Below magnet -1mm"),
    ([0, 2, 0], "Magnet center"),
    ([0, 3, 0], "Above magnet +1mm"),
    ([1, 0, 0], "Side +X 1mm"),
    ([0, 0, 1], "Side +Z 1mm"),
]

# Create mesh
mesh_domain = 6.0e-3
air_region = Box((-mesh_domain, -mesh_domain, -mesh_domain), 
                 (mesh_domain, mesh_domain, mesh_domain)).mat("air")
mesh = air_region.GenerateMesh(maxh=1.0e-3)

# Create CoefficientFunction (NO GridFunction!)
B_cf = radia_ngsolve.RadiaField(magnet_base, 'b')

print("="*70)
print("Comparison: rad.Fld vs rad_ngsolve (Direct CF Evaluation)")
print("="*70)
print(f"{'Point':<25s} {'rad.Fld By':>15s} {'rad_ngsolve By':>15s} {'Error %':>12s}")
print("-"*70)

for point, description in test_points:
    B_radia = rad.Fld(magnet_base, 'b', point)
    point_m = (point[0]/1000, point[1]/1000, point[2]/1000)
    
    # Direct CF evaluation (no GridFunction!)
    B_cf_eval = B_cf(mesh(*point_m))
    
    if abs(B_radia[1]) > 1e-6:
        rel_error = abs(B_radia[1] - B_cf_eval[1]) / abs(B_radia[1]) * 100
    else:
        rel_error = abs(B_radia[1] - B_cf_eval[1]) * 100
    
    print(f"{description:<25s} {B_radia[1]:15.6f} {B_cf_eval[1]:15.6f} {rel_error:11.4f}%")
    
print("\n" + "="*70)
print("Full Vector Comparison")
print("="*70)
for point, description in test_points:
    B_radia = rad.Fld(magnet_base, 'b', point)
    point_m = (point[0]/1000, point[1]/1000, point[2]/1000)
    B_cf_eval = B_cf(mesh(*point_m))
    
    print(f"\n{description}:")
    print(f"  rad.Fld:     [{B_radia[0]:10.6f}, {B_radia[1]:10.6f}, {B_radia[2]:10.6f}]")
    print(f"  rad_ngsolve: [{B_cf_eval[0]:10.6f}, {B_cf_eval[1]:10.6f}, {B_cf_eval[2]:10.6f}]")
