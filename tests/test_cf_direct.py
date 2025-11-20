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

# Test point: magnet center (0, 0.002, 0) m
point = [0, 0.002, 0]

print("Test Point: (0, 0.002, 0) m - Magnet Center")
print("="*60)

# 1. Radia direct evaluation
B_radia = rad.Fld(magnet_base, 'b', point)
print(f"rad.Fld result:          {B_radia}")

# 2. Create CoefficientFunction
B_cf = rad_ngsolve.RadiaField(magnet_base, 'b')
print(f"CoefficientFunction created: {B_cf}")

# 3. Try to evaluate CF at a point (this might not work directly)
print("\nAttempting direct CF evaluation:")
try:
    # CoefficientFunction cannot be evaluated without a mesh
    # We need to create a minimal mesh first
    mesh_domain = 6.0e-3
    air_region = Box((-mesh_domain, -mesh_domain, -mesh_domain),
                     (mesh_domain, mesh_domain, mesh_domain)).mat("air")
    mesh = air_region.GenerateMesh(maxh=1.0e-3)

    # Create a mesh point
    mip = mesh(*point)
    print(f"Mesh point created: {mip}")
    
    # Evaluate CF directly at mesh point
    B_cf_result = B_cf(mip)
    print(f"B_cf(mip) result:        {B_cf_result}")
    
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()

# 4. Use GridFunction approach
print("\n" + "="*60)
print("GridFunction approach:")
try:
    fes = VectorH1(mesh, order=2, dim=3)
    gf_B = GridFunction(fes)
    gf_B.Set(B_cf)
    
    B_gf_result = gf_B(mip)
    print(f"gf_B(mip) result:        {B_gf_result}")
    
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()
