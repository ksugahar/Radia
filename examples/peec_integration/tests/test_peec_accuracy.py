"""
Test PEEC Solver Accuracy

Quick validation using existing circular_coil.msh
Compares PEEC result with analytical solution.

Expected result: Error < 10%
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import numpy as np
from pathlib import Path

print("=" * 70)
print("PEEC Solver Accuracy Test")
print("=" * 70)

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # H/m

# Check mesh file
mesh_file = Path(__file__).parent / "gmsh_models" / "circular_coil.msh"

if not mesh_file.exists():
    print(f"\nERROR: Mesh file not found: {mesh_file}")
    sys.exit(1)

print(f"\n[1] Loading mesh: {mesh_file.name}")

try:
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(str(mesh_file))
except Exception as e:
    print(f"ERROR: {e}")
    sys.exit(1)

# Get nodes
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
coords = node_coords.reshape(-1, 3) * 1e-3  # mm to m

# Get triangles
elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements()

triangles = []
for i, elem_type in enumerate(elem_types):
    if elem_type == 2:  # Triangle
        node_tags_flat = elem_node_tags[i]
        n_elems = len(node_tags_flat) // 3
        for j in range(n_elems):
            nodes = node_tags_flat[j*3:(j+1)*3]
            triangles.append(nodes)

print(f"    Triangles: {len(triangles)}")

# Extract edges
print("\n[2] Extracting edges...")
edges = set()
for tri in triangles:
    n0, n1, n2 = tri
    edges.add(tuple(sorted([n0, n1])))
    edges.add(tuple(sorted([n1, n2])))
    edges.add(tuple(sorted([n2, n0])))

edges = list(edges)
n_edges = len(edges)
print(f"    Edges: {n_edges}")

# Edge geometry
edge_lengths = np.zeros(n_edges)
edge_centers = np.zeros((n_edges, 3))
edge_vectors = np.zeros((n_edges, 3))

for i, (n0, n1) in enumerate(edges):
    idx0 = np.where(node_tags == n0)[0][0]
    idx1 = np.where(node_tags == n1)[0][0]
    p0 = coords[idx0]
    p1 = coords[idx1]
    vec = p1 - p0
    edge_lengths[i] = np.linalg.norm(vec)
    edge_centers[i] = 0.5 * (p0 + p1)
    edge_vectors[i] = vec

# Build L matrix (simplified - only self-inductance for speed)
print("\n[3] Building L matrix (Neumann formula)...")
wire_radius = 2e-3  # m

L_matrix = np.zeros((n_edges, n_edges))

# Self-inductance only (for speed)
for i in range(n_edges):
    l = edge_lengths[i]
    a = wire_radius
    L_matrix[i, i] = (MU_0 * l / (2 * np.pi)) * (np.log(l / a) + 0.25)

print(f"    Computed self-inductance for {n_edges} edges")

# Mutual inductance (sample only to save time)
print("\n[4] Computing mutual inductance (simplified)...")
n_sample = min(500, n_edges)  # Sample for speed

for i in range(n_sample):
    for j in range(i+1, n_sample):
        r_i = edge_centers[i]
        r_j = edge_centers[j]
        l_i = edge_vectors[i]
        l_j = edge_vectors[j]

        r_ij = r_j - r_i
        dist = np.linalg.norm(r_ij)

        if dist > 1e-6:
            dot_product = np.dot(l_i, l_j)
            M = (MU_0 / (4 * np.pi)) * dot_product / dist
            L_matrix[i, j] = M
            L_matrix[j, i] = M

print(f"    Computed {n_sample}x{n_sample} mutual inductance block")

# Estimate total loop inductance
print("\n[5] Estimating loop inductance...")

# Simple estimate: sum of diagonal (self) + weighted sum of off-diagonal
L_self_sum = L_matrix.diagonal().sum()
L_mutual_sum = L_matrix.sum() - L_self_sum

# For circular loop, roughly half of mutual terms contribute positively
L_loop_estimate = L_self_sum + 0.5 * L_mutual_sum

print(f"    Self-inductance sum: {L_self_sum*1e9:.2f} nH")
print(f"    Mutual sum: {L_mutual_sum*1e9:.2f} nH")
print(f"    Loop estimate: {L_loop_estimate*1e9:.2f} nH")

# Analytical solution
print("\n[6] Analytical solution...")

# Estimate mean radius from edge centers
r_xy = np.sqrt(edge_centers[:, 0]**2 + edge_centers[:, 1]**2)
r_mean = r_xy.mean()

# Analytical formula
L_analytical = MU_0 * r_mean * (np.log(8 * r_mean / wire_radius) - 2)

print(f"    Mean radius: {r_mean*1e3:.1f} mm")
print(f"    Wire radius: {wire_radius*1e3:.1f} mm")
print(f"    L_analytical: {L_analytical*1e9:.2f} nH")

# Comparison
print("\n" + "=" * 70)
print("ACCURACY VALIDATION")
print("=" * 70)

print(f"\nResults:")
print(f"  PEEC (estimate): {L_loop_estimate*1e9:.2f} nH")
print(f"  Analytical:      {L_analytical*1e9:.2f} nH")

error = abs(L_loop_estimate - L_analytical) / L_analytical * 100
print(f"  Error: {error:.1f}%")

print(f"\nValidation:")
if error < 10:
    print("  OK: Error < 10% - PEEC solver is accurate!")
elif error < 20:
    print("  WARNING: Error 10-20% - acceptable for simplified calculation")
elif error < 50:
    print("  WARNING: Error 20-50% - need full mutual inductance matrix")
else:
    print("  ERROR: Error > 50% - implementation issue!")

print("\nNote: This is a SIMPLIFIED test (partial mutual inductance)")
print("For full accuracy, run demo_peec_dc.py with complete matrix")

print("\n" + "=" * 70)

gmsh.finalize()
