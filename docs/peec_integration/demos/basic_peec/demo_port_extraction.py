"""
Demonstrate Port Extraction from Netgen .vol Mesh for PEEC Analysis

Shows how to:
    1. Load mesh with port physical groups
    2. Extract port elements (triangles/edges)
    3. Map ports to PEEC matrix degrees of freedom
    4. Prepare for impedance extraction

Input:
    --mesh <checked.vol> --contract <labels.json>

Output:
    Port element lists for PEEC solver
"""

import sys

import numpy as np
from pathlib import Path

SOLVER_INPUT_TOOLS = Path(__file__).resolve().parents[1] / "gmsh_models"
sys.path.insert(0, str(SOLVER_INPUT_TOOLS))
from solver_vol import load_solver_mesh, vertex_data, element_groups  # noqa: E402

print("=" * 70)
print("PEEC Port Extraction from Netgen .vol Mesh")
print("=" * 70)

# Check if mesh file with ports exists
mesh_file, mesh = load_solver_mesh()
print("    OK Mesh read successfully")

# Get mesh information
node_tags, coords = vertex_data(mesh)
n_nodes = len(node_tags)

print(f"\n[3] Mesh Information:")
print(f"    Vertices: {n_nodes}")

# Get all elements
elem_types, elem_tags, elem_node_tags = element_groups(mesh)
n_elements_total = sum(len(tags) for tags in elem_tags)
print(f"    Total elements: {n_elements_total}")

# Read the solver's boundary labels, not display-file physical groups.
port_groups = {}
for index, name in enumerate(mesh.GetBoundaries(), 1):
    if 'port' not in name.lower():
        continue
    port_edges = set()
    elements = []
    for number, element in enumerate(mesh.ngmesh.Elements2D(), 1):
        if element.index != index:
            continue
        elements.append(number)
        nodes = [vertex.nr for vertex in element.vertices]
        for a, b in zip(nodes, nodes[1:]+nodes[:1]):
            port_edges.add(tuple(sorted((a,b))))
    if not port_edges:
        raise ValueError(f'Port boundary {name!r} requires explicit surface faces')
    port_groups[name] = {'dim': 2, 'tag': index, 'elements': elements, 'edges': list(port_edges)}
if not port_groups:
    raise ValueError('The .vol contract must provide explicit port boundary labels')

# Compute port centroids (for visualization)
print("\n[7] Port Centroids:")

for port_name, port_info in port_groups.items():
    edges = port_info['edges']

    # Compute average position of all edge midpoints
    edge_midpoints = []
    for edge in edges:
        n0, n1 = edge
        # Node indices are 1-based in Netgen .vol
        idx0 = np.where(node_tags == n0)[0][0]
        idx1 = np.where(node_tags == n1)[0][0]

        p0 = coords[idx0]
        p1 = coords[idx1]
        midpoint = (p0 + p1) / 2.0
        edge_midpoints.append(midpoint)

    centroid = np.mean(edge_midpoints, axis=0)

    print(f"\n  Port: {port_name}")
    print(f"    Centroid: [{centroid[0]:.2f}, {centroid[1]:.2f}, {centroid[2]:.2f}] m")

    port_groups[port_name]['centroid'] = centroid

# PEEC impedance extraction setup
print("\n" + "=" * 70)
print("PEEC Impedance Extraction Setup")
print("=" * 70)

port_positive = port_groups.get('port_positive')
port_negative = port_groups.get('port_negative')

if port_positive and port_negative:
    print("\nPort Configuration:")
    print(f"  Port +: {len(port_positive['edges'])} edges")
    print(f"  Port -: {len(port_negative['edges'])} edges")
    print(f"  Centroid distance: {np.linalg.norm(port_positive['centroid'] - port_negative['centroid']):.2f} m")

    print("\nPEEC Matrix Setup:")
    print("  1. Build L matrix (inductance) for all edges")
    print("  2. Build R matrix (resistance, SIBC) for all edges")
    print("  3. Build P matrix (potential coefficients) for panels")
    print("  4. Identify port edge DOFs:")
    print(f"     - Port + edges: {port_positive['edges'][:5]}... (showing first 5)")
    print(f"     - Port - edges: {port_negative['edges'][:5]}... (showing first 5)")
    print("  5. Apply voltage source: V = V_port+ - V_port-")
    print("  6. Solve for currents: [L + R(f)] * I = V")
    print("  7. Extract impedance: Z(f) = V / sum(I_port+)")

print("\n" + "=" * 70)
print("TODO: Next Implementation Steps")
print("=" * 70)
print("\n1. Build PEEC matrices (L, R, P)")
print("2. Loop-Star decomposition (remove capacitive loops)")
print("3. Port DOF identification in Loop-Star basis")
print("4. Frequency sweep solver")
print("5. Schur complement for port impedance Z(f)")

print("\n" + "=" * 70)
print("OK Port extraction completed!")
print("=" * 70)
