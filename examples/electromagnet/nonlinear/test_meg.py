"""Test MEG file parsing and element creation."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import numpy as np
import radia as rad
from meg_to_vol import parse_meg_file

print("=" * 60)
print("MEG File Parsing Test")
print("=" * 60)

rad.UtiDelAll()
rad.FldUnits('mm')

meg_file = r"S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/nonlinear_20000AT/ELF_MMB8T_EIEM2_1x1x1/ELF_magic.meg"

print(f"\nParsing: {meg_file}")
nodes, elements, scale = parse_meg_file(meg_file)

print(f"  Nodes: {len(nodes)}")
print(f"  Elements: {len(elements)}")
print(f"  Scale: {scale}")

# Convert back to mm
for nid in nodes:
    nodes[nid] = [c / scale for c in nodes[nid]]

# Get yoke elements only
yoke_elements = [e for e in elements if e['type'] == 'MMB8T']
print(f"  Yoke elements: {len(yoke_elements)}")

# Load BH curve
print("\nLoading BH curve...")
bh_data = np.loadtxt("BH.txt", comments='#').tolist()
mat = rad.MatSatIsoTab(bh_data)
print(f"  Material: {mat}")

# Create hex elements one by one
print("\nCreating hex elements:")
yoke_objs = []

def get_centroid(nodes, node_ids):
    coords = np.array([nodes[nid] for nid in node_ids])
    return coords.mean(axis=0)

# ELF to Radia vertex reorder
reorder = [0, 1, 3, 2, 4, 5, 7, 6]

for elem in yoke_elements:
    elem_id = elem['id']
    node_ids = elem['nodes']

    # Get vertices
    vertices = [list(nodes[node_ids[i]]) for i in reorder]
    centroid = get_centroid(nodes, node_ids)

    print(f"\n  Element {elem_id}:")
    print(f"    Node IDs: {node_ids}")
    print(f"    Centroid: [{centroid[0]:.1f}, {centroid[1]:.1f}, {centroid[2]:.1f}]")

    # Print vertices
    for i, v in enumerate(vertices):
        print(f"    V{i}: [{v[0]:.1f}, {v[1]:.1f}, {v[2]:.1f}]")

    try:
        hex_obj = rad.ObjHexahedron(vertices, [0, 0, 0])
        rad.MatApl(hex_obj, mat)
        yoke_objs.append(hex_obj)
        print(f"    Created: {hex_obj}")
    except Exception as e:
        print(f"    FAILED: {e}")

print(f"\n\nCreated {len(yoke_objs)} elements")

if len(yoke_objs) > 0:
    print("\nCreating container...")
    yoke = rad.ObjCnt(yoke_objs)
    print(f"  Container: {yoke}")

    print("\nApplying symmetry...")
    print("  X mirror...")
    rad.TrfMlt(yoke, rad.TrfPlSym([0, 0, 0], [1, 0, 0]), 2)
    print("  Z mirror...")
    rad.TrfMlt(yoke, rad.TrfPlSym([0, 0, 0], [0, 0, 1]), 2)
    print("  Done")

    print("\nSolving...")
    result = rad.Solve(yoke, 0.0001, 100, 0)
    print(f"  Result: {result}")

print("\nDone!")
