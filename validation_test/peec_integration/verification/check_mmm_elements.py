"""Check MMM element count and geometry."""

import sys
sys.path.insert(0, '../../../src/radia')
import radia as rad


def hex_vertices(cx, cy, cz, dx, dy, dz):
    """Generate hexahedron vertices from center and dimensions."""
    hx, hy, hz = dx/2, dy/2, dz/2
    return [
        [cx-hx, cy-hy, cz-hz], [cx+hx, cy-hy, cz-hz],
        [cx+hx, cy+hy, cz-hz], [cx-hx, cy+hy, cz-hz],
        [cx-hx, cy-hy, cz+hz], [cx+hx, cy-hy, cz+hz],
        [cx+hx, cy+hy, cz+hz], [cx-hx, cy+hy, cz+hz]
    ]


rad.UtiDelAll()

# Create hexahedral mesh manually (2x2x2 = 8 elements)
cube_half = 0.015  # m
n_div = 2
elem_size = 2 * cube_half / n_div
n_elements = n_div ** 3

print(f"Creating {n_div}x{n_div}x{n_div} = {n_elements} hexahedral elements...")

elements = []
for i in range(n_div):
    for j in range(n_div):
        for k in range(n_div):
            cx = -cube_half + (i + 0.5) * elem_size
            cy = -cube_half + (j + 0.5) * elem_size
            cz = -cube_half + (k + 0.5) * elem_size
            vertices = hex_vertices(cx, cy, cz, elem_size, elem_size, elem_size)
            elem = rad.ObjHexahedron(vertices, [0, 0, 0])
            elements.append(elem)
            print(f"  Element {len(elements)}: center=({cx:.4f}, {cy:.4f}, {cz:.4f})")

core = rad.ObjCnt(elements)
print(f"Container handle: {core}")

# Apply material
mat = rad.MatLin(1000)
rad.MatApl(core, mat)
print(f"Applied material with mu_r=1000")

# Get geometry info
try:
    geo = rad.ObjGeoVol(core)
    print(f"Volume: {geo} m^3")
except Exception as e:
    print(f"ObjGeoVol failed: {e}")

# Check element count via container
try:
    cnt_size = rad.ObjCntSize(core)
    print(f"Container size: {cnt_size}")
except Exception as e:
    print(f"Container check failed: {e}")

# Solve
print("\n--- Solving ---")
result = rad.Solve(core, 0.0001, 10, 0)
print(f"Solve result: {result}")

# Try to get magnetization
try:
    M = rad.Fld(core, 'm', [0, 0, 0])
    print(f"Magnetization at center: {M}")
except Exception as e:
    print(f"Fld failed: {e}")
