"""
Demonstrate Port Extraction from GMSH Mesh for PEEC Analysis

Shows how to:
    1. Load mesh with port physical groups
    2. Extract port elements (triangles/edges)
    3. Map ports to PEEC matrix degrees of freedom
    4. Prepare for impedance extraction

Input:
    gmsh_models/circular_coil_with_ports.msh

Output:
    Port element lists for PEEC solver
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import numpy as np
from pathlib import Path

GMSH_MODELS = Path(__file__).resolve().parents[1] / "gmsh_models"
sys.path.insert(0, str(GMSH_MODELS))
from gmsh_ascii_reader import read_gmsh_ascii  # noqa: E402

print("=" * 70)
print("PEEC Port Extraction from GMSH Mesh")
print("=" * 70)

# Check if mesh file with ports exists
mesh_file = GMSH_MODELS / "circular_coil_with_ports.msh"

if not mesh_file.exists():
    print(f"\nERROR: Mesh file not found: {mesh_file}")
    print("\n" + "=" * 70)
    print("Please generate mesh with ports first:")
    print("=" * 70)
    print("\nStep 1: Navigate to gmsh_models directory")
    print("  cd gmsh_models")
    print("\nStep 2: Run Cubit mesh generation with ports")
    print('  "<Coreform Cubit 2025.8+>/bin/python3/python.exe" generate_coil_with_ports.py')
    print("\nStep 3: Come back and run this script")
    print("  cd ..")
    print("  python demo_port_extraction.py")
    sys.exit(1)

print(f"\n[1] Loading mesh file: {mesh_file.name}")
print(f"    Path: {mesh_file}")

print("\n[2] Reading GMSH ASCII fixture directly...")
mesh = read_gmsh_ascii(mesh_file)
print("    OK Mesh read successfully")

# Get mesh information
node_tags = mesh.node_tags
n_nodes = len(node_tags)
coords = mesh.coords

print(f"\n[3] Mesh Information:")
print(f"    Vertices: {n_nodes}")

# Get all elements
elem_types, elem_tags, elem_node_tags = mesh.grouped_elements()
n_elements_total = sum(len(tags) for tags in elem_tags)
print(f"    Total elements: {n_elements_total}")

# Extract physical groups
print("\n[4] Physical Groups (Ports):")
phys_groups = mesh.physical_groups()

port_groups = {}
for dim, tag in phys_groups:
    name = mesh.physical_name(dim, tag)
    print(f"    Dimension {dim}, Tag {tag}: {name}")

    # Store port groups
    if 'port' in name.lower():
        port_groups[name] = {'dim': dim, 'tag': tag}

# Extract port elements
print("\n[5] Extracting Port Elements:")

for port_name, port_info in port_groups.items():
    dim = port_info['dim']
    tag = port_info['tag']

    print(f"\n  Port: {port_name}")
    print(f"    Physical tag: {tag}")

    port_elements = [
        elem.tag
        for elem in mesh.elements_for_physical_group(dim, tag)
        if elem.element_type == 2
    ]

    print(f"    Elements in port: {len(port_elements)}")

    # Store for PEEC solver
    port_groups[port_name]['elements'] = port_elements

# Extract port edges (for PEEC DOF mapping)
print("\n[6] Port Edge Extraction (for PEEC DOF):")
print("    PEEC uses EDGES as degrees of freedom")

for port_name, port_info in port_groups.items():
    port_elements = port_info['elements']

    # Collect all edges from port triangles
    port_edges = set()

    for elem_tag in port_elements:
        # Get nodes for this element
        elem = mesh.element(elem_tag)
        elem_nodes = elem.nodes

        # Triangle has 3 edges: (n0,n1), (n1,n2), (n2,n0)
        n0, n1, n2 = elem_nodes

        # Create edges (sorted for uniqueness)
        edge1 = tuple(sorted([n0, n1]))
        edge2 = tuple(sorted([n1, n2]))
        edge3 = tuple(sorted([n2, n0]))

        port_edges.add(edge1)
        port_edges.add(edge2)
        port_edges.add(edge3)

    print(f"\n  Port: {port_name}")
    print(f"    Unique edges: {len(port_edges)}")

    # Store edge list
    port_groups[port_name]['edges'] = list(port_edges)

# Compute port centroids (for visualization)
print("\n[7] Port Centroids:")

for port_name, port_info in port_groups.items():
    edges = port_info['edges']

    # Compute average position of all edge midpoints
    edge_midpoints = []
    for edge in edges:
        n0, n1 = edge
        # Node indices are 1-based in GMSH
        idx0 = np.where(node_tags == n0)[0][0]
        idx1 = np.where(node_tags == n1)[0][0]

        p0 = coords[idx0]
        p1 = coords[idx1]
        midpoint = (p0 + p1) / 2.0
        edge_midpoints.append(midpoint)

    centroid = np.mean(edge_midpoints, axis=0)

    print(f"\n  Port: {port_name}")
    print(f"    Centroid: [{centroid[0]:.2f}, {centroid[1]:.2f}, {centroid[2]:.2f}] mm")

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
    print(f"  Centroid distance: {np.linalg.norm(port_positive['centroid'] - port_negative['centroid']):.2f} mm")

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
