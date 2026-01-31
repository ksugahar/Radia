"""Test creating one hexahedron at a time."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import numpy as np
import radia as rad
from meg_to_vol import parse_meg_file

print("=" * 60)
print("Test: Create hexahedra one by one")
print("=" * 60)

rad.UtiDelAll()
rad.FldUnits('mm')

meg_file = r"S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/nonlinear_20000AT/ELF_MMB8T_EIEM2_1x1x1/ELF_magic.meg"
nodes, elements, scale = parse_meg_file(meg_file)

# Convert to mm
for nid in nodes:
    nodes[nid] = [c / scale for c in nodes[nid]]

yoke_elements = [e for e in elements if e['type'] == 'MMB8T']
print(f"Yoke elements: {len(yoke_elements)}")

# Load BH curve
bh_data = np.loadtxt("BH.txt", comments='#').tolist()
mat = rad.MatSatIsoTab(bh_data)

# Different orderings to try
orderings = {
    'elf_default': [0, 1, 3, 2, 4, 5, 7, 6],
    'direct': [0, 1, 2, 3, 4, 5, 6, 7],
    'swapped': [0, 3, 2, 1, 4, 7, 6, 5],
}

# Test just the first element with different orderings
elem = yoke_elements[0]
print(f"\nTesting element {elem['id']}")
node_ids = elem['nodes']

for name, reorder in orderings.items():
    print(f"\n  Ordering '{name}': {reorder}")
    rad.UtiDelAll()
    mat = rad.MatSatIsoTab(bh_data)

    vertices = [list(nodes[node_ids[i]]) for i in reorder]
    print(f"    Vertices:")
    for i, v in enumerate(vertices):
        print(f"      V{i}: [{v[0]:.1f}, {v[1]:.1f}, {v[2]:.1f}]")

    try:
        hex_obj = rad.ObjHexahedron(vertices, [0, 0, 0])
        rad.MatApl(hex_obj, mat)
        print(f"    Created: {hex_obj}")

        # Try to compute field
        centroid = [sum(v[i] for v in vertices) / 8 for i in range(3)]
        B = rad.Fld(hex_obj, 'b', centroid)
        print(f"    Field at centroid: {B}")
    except Exception as e:
        print(f"    FAILED: {e}")

print("\nNow testing all 13 elements with direct ordering...")
rad.UtiDelAll()
mat = rad.MatSatIsoTab(bh_data)

success_count = 0
for elem in yoke_elements:
    elem_id = elem['id']
    node_ids = elem['nodes']

    # Use direct ordering (no reorder)
    vertices = [list(nodes[nid]) for nid in node_ids]

    try:
        hex_obj = rad.ObjHexahedron(vertices, [0, 0, 0])
        rad.MatApl(hex_obj, mat)
        success_count += 1
        print(f"  Element {elem_id}: OK ({hex_obj})")
    except Exception as e:
        print(f"  Element {elem_id}: FAILED - {e}")

print(f"\nSuccess: {success_count}/{len(yoke_elements)}")

print("\nDone!")
