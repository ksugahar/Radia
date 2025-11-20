import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/python'))

from numpy import *
from ngsolve import *
from netgen.occ import *
import radia as rad
import radia_ngsolve

# Set Radia to use meters (required for NGSolve integration)
rad.FldUnits('m')

# ========================================================================
# Radia model
# ========================================================================

magnet_size = 0.001  # m (1 mm)
magnetization = 1.0  # T

# Magnet center at (0, 0.002, 0) m, size 0.001x0.001x0.001 m
# So magnet occupies: x=[-0.0005, 0.0005], y=[0.0015, 0.0025], z=[-0.0005, 0.0005]
magnet_base = rad.ObjRecMag([0, 0.002, 0], [magnet_size, magnet_size, magnet_size], [0, magnetization, 0])

test_points = [
    ([0, 0, 0], "Origin (0, 0, 0)"),
    ([0, 0.001, 0], "Below magnet -1mm"),
    ([0, 0.002, 0], "Magnet center"),
    ([0, 0.003, 0], "Above magnet +1mm"),
    ([0.001, 0, 0], "Side +X 1mm"),
    ([0, 0, 0.001], "Side +Z 1mm"),
]

print("="*70)
print("Radia Field Evaluation (rad.Fld) - Direct")
print("="*70)
for point, description in test_points:
	B_radia = rad.Fld(magnet_base, 'b', point)
	print(f"{description:25s}: Bx={B_radia[0]:12.6f} By={B_radia[1]:12.6f} Bz={B_radia[2]:12.6f}")

# Create NGSolve mesh
mesh_domain = 6.0e-3  # 6mm in meters
air_region = Box((-mesh_domain, -mesh_domain, -mesh_domain), 
                 (mesh_domain, mesh_domain, mesh_domain)).mat("air")
mesh_maxh = 1.0e-3  # 1mm in meters
mesh = air_region.GenerateMesh(maxh=mesh_maxh)

# Create CoefficientFunction and GridFunction
fes = VectorH1(mesh, order=2, dim=3)
B_cf = radia_ngsolve.RadiaField(magnet_base, 'b')
gf_B = GridFunction(fes)
gf_B.Set(B_cf)

print("\n" + "="*70)
print("NGSolve GridFunction Evaluation")
print("="*70)
for point, description in test_points:
	try:
		B_ngsolve = gf_B(mesh(*point))
		print(f"{description:25s}: Bx={B_ngsolve[0]:12.6f} By={B_ngsolve[1]:12.6f} Bz={B_ngsolve[2]:12.6f}")
	except Exception as e:
		print(f"{description:25s}: ERROR - {e}")

print("\n" + "="*70)
print("Comparison Table")
print("="*70)
print(f"{'Point':<25s} {'rad.Fld By':>15s} {'NGSolve By':>15s} {'Error %':>15s}")
print("-"*70)
for point, description in test_points:
	B_radia = rad.Fld(magnet_base, 'b', point)
	try:
		B_ngsolve = gf_B(mesh(*point))
		if abs(B_radia[1]) > 1e-6:
			rel_error = abs(B_radia[1] - B_ngsolve[1]) / abs(B_radia[1]) * 100
		else:
			rel_error = abs(B_radia[1] - B_ngsolve[1]) * 100
		print(f"{description:<25s} {B_radia[1]:15.6f} {B_ngsolve[1]:15.6f} {rel_error:14.2f}%")
	except Exception as e:
		print(f"{description:<25s} {'ERROR':<30s}")
