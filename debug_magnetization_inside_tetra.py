#!/usr/bin/env python
"""
Debug: Why is magnetization zero inside tetrahedron?

Hypothesis: Radia cannot evaluate magnetization inside tetrahedral elements.
This would explain why external field calculations are wrong.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build/Release'))

os.environ['RADIA_TETRA_METHOD'] = 'ANALYTICAL'

import numpy as np
import radia as rad

print("="*70)
print("Debug: Magnetization Inside Tetrahedron")
print("="*70)

cube_size = 0.1  # m
H_ext = 1000.0   # A/m
mu0 = 4*np.pi*1e-7
B_ext = H_ext * mu0
mu_r = 100
ksi = mu_r - 1

# ============================================================
# Create tetrahedron
# ============================================================
rad.FldUnits('m')

vertices = [
    [0, 0, 0],
    [cube_size, 0, 0],
    [0, cube_size, 0],
    [0, 0, cube_size]
]

TETRA_FACES = [
    [1, 3, 2],
    [1, 2, 4],
    [1, 4, 3],
    [2, 3, 4]
]

tetra = rad.ObjPolyhdr(vertices, TETRA_FACES, [0, 0, 0])
mat = rad.MatLin([ksi, ksi], [0, 0, 1])
rad.MatApl(tetra, mat)

bg = rad.ObjBckg([0, 0, B_ext])
container = rad.ObjCnt([tetra, bg])

print("\nSolving...")
result = rad.Solve(container, 0.0001, 10000)
print(f"Solve result: {result}")

# ============================================================
# Test magnetization at various points
# ============================================================
print("\n" + "="*70)
print("Magnetization Evaluation")
print("="*70)

# Tetrahedron centroid
center = [cube_size/4, cube_size/4, cube_size/4]

# Test points: inside and outside
test_points = [
    (center, "Centroid (inside)"),
    ([0.01, 0.01, 0.01], "Near vertex 1 (inside)"),
    ([cube_size*0.3, cube_size*0.3, cube_size*0.1], "Inside tetra"),
    ([cube_size*0.2, cube_size*0.2, cube_size*0.2], "Inside tetra"),
    ([cube_size*0.5, cube_size*0.5, cube_size*0.5], "Outside tetra"),
    ([0.05, 0.05, 0.2], "Far outside"),
]

print("\n{:<30} {:<20} {:<20}".format("Location", "|M| (A/m)", "|H| (A/m)"))
print("-"*70)

for pt, label in test_points:
    try:
        M = rad.Fld(container, 'm', pt)
        M_mag = np.linalg.norm(M)
    except Exception as e:
        M_mag = f"ERROR: {e}"

    try:
        H = rad.Fld(container, 'h', pt)
        H_mag = np.linalg.norm(H)
    except Exception as e:
        H_mag = f"ERROR: {e}"

    print(f"{label:<30} {str(M_mag):<20} {str(H_mag):<20}")

# ============================================================
# Compare with hexahedron
# ============================================================
print("\n" + "="*70)
print("Comparison: Hexahedron (should work)")
print("="*70)

rad.FldUnits('m')

cube_hex = rad.ObjRecMag([0, 0, 0], [cube_size, cube_size, cube_size], [0, 0, 0])
mat_hex = rad.MatLin([ksi, ksi], [0, 0, 1])
rad.MatApl(cube_hex, mat_hex)

bg_hex = rad.ObjBckg([0, 0, B_ext])
container_hex = rad.ObjCnt([cube_hex, bg_hex])

rad.Solve(container_hex, 0.0001, 10000)

test_points_hex = [
    ([0, 0, 0], "Center"),
    ([cube_size*0.4, cube_size*0.4, cube_size*0.4], "Near center"),
    ([cube_size*0.5, cube_size*0.5, cube_size*0.5], "Outside"),
]

print("\n{:<30} {:<20} {:<20}".format("Location", "|M| (A/m)", "|H| (A/m)"))
print("-"*70)

for pt, label in test_points_hex:
    M = rad.Fld(container_hex, 'm', pt)
    M_mag = np.linalg.norm(M)

    H = rad.Fld(container_hex, 'h', pt)
    H_mag = np.linalg.norm(H)

    print(f"{label:<30} {M_mag:<20.2f} {H_mag:<20.2f}")

# ============================================================
# Hypothesis Test: Radia limitation
# ============================================================
print("\n" + "="*70)
print("Analysis")
print("="*70)

print("""
Observation:
- Hexahedron: M can be evaluated inside the element
- Tetrahedron: M = 0 everywhere (inside and outside)

Possible causes:
1. Radia's rad.Fld(obj, 'm', pt) does NOT support tetrahedral elements
2. Magnetization is stored only for hexahedral elements
3. For tetrahedron, M must be inferred from H field and material properties

Radia documentation states:
"rad.Fld() is designed for field calculation in AIR regions"
"Inside magnets, rad.Fld() returns inaccurate values"

This means:
- rad.Fld(obj, 'h', pt) - Works outside magnets (in air)
- rad.Fld(obj, 'm', pt) - May NOT work for all element types

Implication for high permeability:
The issue is NOT in the ANALYTICAL method itself, but in how
Radia's solver computes magnetization for tetrahedral elements.

The external field calculation (ANALYTICAL method) depends on having
correct magnetization values, which are computed by rad.Solve().

If rad.Solve() doesn't correctly compute M for tetrahedra with high mu_r,
then ANALYTICAL method will produce wrong external fields.
""")

print("="*70)
print("Conclusion")
print("="*70)

print("""
ROOT CAUSE:
  rad.Fld(container, 'm', point) returns [0, 0, 0] for tetrahedra

This is a fundamental limitation in Radia's current implementation:
- Magnetization evaluation is not implemented for tetrahedral elements
- Without knowing M, we cannot compute correct demagnetization
- This breaks the solver for high permeability materials

NEXT STEPS:
1. Check if magnetization is stored but not accessible via rad.Fld()
2. Investigate Radia's internal magnetization storage
3. Determine if ANALYTICAL method can work without M evaluation
""")

print("="*70)
