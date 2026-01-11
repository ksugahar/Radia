"""Check MMM element count and geometry."""

import sys
sys.path.insert(0, '../../src/radia')
import radia as rad

rad.UtiDelAll()
rad.FldUnits('m')

# Create hexahedron
vertices = [
    [-0.015, -0.015, -0.015],
    [0.015, -0.015, -0.015],
    [0.015, 0.015, -0.015],
    [-0.015, 0.015, -0.015],
    [-0.015, -0.015, 0.015],
    [0.015, -0.015, 0.015],
    [0.015, 0.015, 0.015],
    [-0.015, 0.015, 0.015],
]
core = rad.ObjHexahedron(vertices, [0, 0, 0])
print(f"Hexahedron handle: {core}")

# Apply material
mat = rad.MatLin(1000)
rad.MatApl(core, mat)

# Try to subdivide
try:
    # ObjDivMag divides magnetic objects
    rad.ObjDivMag(core, [2, 2, 2])
    print("Subdivided into 2x2x2 = 8 elements")
except Exception as e:
    print(f"ObjDivMag failed: {e}")

# Try alternative subdivision
try:
    rad.ObjDivMagPln(core, [0, 0, 1], [0, 0, 0])
    print("Subdivided by plane")
except Exception as e:
    print(f"ObjDivMagPln failed: {e}")

# Get geometry info
try:
    geo = rad.ObjGeoVol(core)
    print(f"Volume: {geo} m^3")
except Exception as e:
    print(f"ObjGeoVol failed: {e}")

# Check element count via container
try:
    cnt = rad.ObjCnt([core])
    elements = rad.ObjCntSize(cnt)
    print(f"Container size: {elements}")
except Exception as e:
    print(f"Container check failed: {e}")

# Check with interaction
print("\n--- Using PreRelax to create interaction ---")
# PreRelax creates the interaction matrix
result = rad.Solve(core, 0.0001, 10, 0)
print(f"Solve result: {result}")

# Try to get magnetization
try:
    M = rad.Fld(core, 'm', [0, 0, 0])
    print(f"Magnetization at center: {M}")
except Exception as e:
    print(f"Fld failed: {e}")
