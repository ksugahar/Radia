import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/python'))

from numpy import *
from ngsolve import *
from netgen.occ import *
import radia as rad
import rad_ngsolve

# Set Radia to use meters (required for NGSolve integration)
rad.FldUnits('m')

# Create magnet
magnet_size = 0.001  # m (1 mm)
magnetization = 1.0  # T
magnet_base = rad.ObjRecMag([0, 0.002, 0], [magnet_size, magnet_size, magnet_size], [0, magnetization, 0])

test_points = [
    ([0, 0, 0], "Origin"),
    ([0, 0.002, 0], "Magnet center"),
    ([0, 0.003, 0], "Above magnet +1mm"),
]

print("="*70)
print("Test: VectorH1 with order=1")
print("="*70)

# Create mesh
mesh_domain = 6.0e-3
air_region = Box((-mesh_domain, -mesh_domain, -mesh_domain), 
                 (mesh_domain, mesh_domain, mesh_domain)).mat("air")
mesh = air_region.GenerateMesh(maxh=1.0e-3)

# Try order=1 (linear elements)
fes = VectorH1(mesh, order=1)
B_cf = rad_ngsolve.RadiaField(magnet_base, 'b')
gf_B = GridFunction(fes)
gf_B.Set(B_cf)

print(f"\nFinite Element Space: {fes}")
print(f"Number of DOFs: {fes.ndof}")

for point, description in test_points:
	B_radia = rad.Fld(magnet_base, 'b', point)
	B_ngsolve = gf_B(mesh(*point))

	print(f"\n{description} ({point[0]}, {point[1]}, {point[2]}) m:")
	print(f"  rad.Fld:    Bx={B_radia[0]:10.6f}, By={B_radia[1]:10.6f}, Bz={B_radia[2]:10.6f}")
	print(f"  rad_ngsolve: Bx={B_ngsolve[0]:10.6f}, By={B_ngsolve[1]:10.6f}, Bz={B_ngsolve[2]:10.6f}")

	# Direct CF evaluation (no GridFunction)
	mip = mesh(*point)
	B_cf_direct = B_cf(mip)
	print(f"  B_cf direct: Bx={B_cf_direct[0]:10.6f}, By={B_cf_direct[1]:10.6f}, Bz={B_cf_direct[2]:10.6f}")
