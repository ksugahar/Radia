"""Test MEG file parsing only - no Radia."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import numpy as np
from meg_to_vol import parse_meg_file

print("=" * 60)
print("MEG File Parsing Test (no Radia)")
print("=" * 60)

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

def get_centroid(nodes, node_ids):
    coords = np.array([nodes[nid] for nid in node_ids])
    return coords.mean(axis=0)

# ELF to Radia vertex reorder
reorder = [0, 1, 3, 2, 4, 5, 7, 6]

print("\nYoke element details:")
for elem in yoke_elements:
    elem_id = elem['id']
    node_ids = elem['nodes']
    vertices = [list(nodes[node_ids[i]]) for i in reorder]
    centroid = get_centroid(nodes, node_ids)

    print(f"\n  Element {elem_id}:")
    print(f"    Centroid: [{centroid[0]:.1f}, {centroid[1]:.1f}, {centroid[2]:.1f}]")
    for i, v in enumerate(vertices):
        print(f"    V{i}: [{v[0]:.1f}, {v[1]:.1f}, {v[2]:.1f}]")

print("\nDone!")
