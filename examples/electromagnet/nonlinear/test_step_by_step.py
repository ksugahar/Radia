"""Step by step test."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))
sys.stdout.flush()

print("Step 1: Import numpy", flush=True)
import numpy as np
print("  OK", flush=True)

print("Step 2: Import radia", flush=True)
import radia as rad
print("  OK", flush=True)

print("Step 3: Import meg_to_vol", flush=True)
from meg_to_vol import parse_meg_file
print("  OK", flush=True)

print("Step 4: UtiDelAll", flush=True)
rad.UtiDelAll()
print("  OK", flush=True)

print("Step 5: FldUnits", flush=True)
rad.FldUnits('mm')
print("  OK", flush=True)

print("Step 6: Parse MEG file", flush=True)
meg_file = r"S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/nonlinear_20000AT/ELF_MMB8T_EIEM2_1x1x1/ELF_magic.meg"
nodes, elements, scale = parse_meg_file(meg_file)
print(f"  OK: {len(nodes)} nodes, {len(elements)} elements", flush=True)

print("Step 7: Convert to mm", flush=True)
for nid in nodes:
    nodes[nid] = [c / scale for c in nodes[nid]]
print("  OK", flush=True)

print("Step 8: Get first yoke element", flush=True)
yoke_elements = [e for e in elements if e['type'] == 'MMB8T']
elem = yoke_elements[0]
node_ids = elem['nodes']
print(f"  Element {elem['id']}, nodes: {node_ids}", flush=True)

print("Step 9: Get vertices", flush=True)
vertices = [list(nodes[nid]) for nid in node_ids]
print(f"  V0: {vertices[0]}", flush=True)
print(f"  V1: {vertices[1]}", flush=True)
print(f"  V2: {vertices[2]}", flush=True)
print(f"  V3: {vertices[3]}", flush=True)
print(f"  V4: {vertices[4]}", flush=True)
print(f"  V5: {vertices[5]}", flush=True)
print(f"  V6: {vertices[6]}", flush=True)
print(f"  V7: {vertices[7]}", flush=True)

print("Step 10: Load BH curve", flush=True)
bh_data = np.loadtxt("BH.txt", comments='#').tolist()
print(f"  OK: {len(bh_data)} points", flush=True)

print("Step 11: Create material", flush=True)
mat = rad.MatSatIsoTab(bh_data)
print(f"  OK: {mat}", flush=True)

print("Step 12: Create hexahedron", flush=True)
try:
    hex_obj = rad.ObjHexahedron(vertices, [0, 0, 0])
    print(f"  OK: {hex_obj}", flush=True)
except Exception as e:
    print(f"  FAILED: {e}", flush=True)
    sys.exit(1)

print("Step 13: Apply material", flush=True)
rad.MatApl(hex_obj, mat)
print("  OK", flush=True)

print("Step 14: Solve", flush=True)
try:
    result = rad.Solve(hex_obj, 0.0001, 100, 0)
    print(f"  OK: {result}", flush=True)
except Exception as e:
    print(f"  FAILED: {e}", flush=True)

print("\nAll steps completed!", flush=True)
