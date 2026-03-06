"""
Debug Edge Lengths to Find Self-Inductance Issue
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import numpy as np
from pathlib import Path

mesh_file = Path(__file__).parent / "gmsh_models" / "circular_coil.msh"

import gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.open(str(mesh_file))

node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
coords = node_coords.reshape(-1, 3) * 1e-3  # mm to m

elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements()

triangles = []
for i, elem_type in enumerate(elem_types):
    if elem_type == 2:
        node_tags_flat = elem_node_tags[i]
        n_elems = len(node_tags_flat) // 3
        for j in range(n_elems):
            nodes = node_tags_flat[j*3:(j+1)*3]
            triangles.append(nodes)

edges = set()
for tri in triangles:
    n0, n1, n2 = tri
    edges.add(tuple(sorted([n0, n1])))
    edges.add(tuple(sorted([n1, n2])))
    edges.add(tuple(sorted([n2, n0])))

edges = list(edges)
n_edges = len(edges)

edge_lengths = []
for n0, n1 in edges:
    idx0 = np.where(node_tags == n0)[0][0]
    idx1 = np.where(node_tags == n1)[0][0]
    p0 = coords[idx0]
    p1 = coords[idx1]
    length = np.linalg.norm(p1 - p0)
    edge_lengths.append(length)

edge_lengths = np.array(edge_lengths)

print("=" * 70)
print("Edge Length Analysis")
print("=" * 70)

print(f"\nNumber of edges: {n_edges}")
print(f"\nEdge lengths (meters):")
print(f"  Min:    {edge_lengths.min()*1e3:.4f} mm")
print(f"  Max:    {edge_lengths.max()*1e3:.4f} mm")
print(f"  Mean:   {edge_lengths.mean()*1e3:.4f} mm")
print(f"  Median: {np.median(edge_lengths)*1e3:.4f} mm")

# Histogram
bins = [0, 0.5e-3, 1e-3, 1.5e-3, 2e-3, 3e-3, 5e-3, 10e-3, 100e-3]
hist, _ = np.histogram(edge_lengths, bins=bins)

print(f"\nLength distribution:")
for i in range(len(bins)-1):
    print(f"  {bins[i]*1e3:.1f} - {bins[i+1]*1e3:.1f} mm: {hist[i]:6d} edges ({hist[i]/n_edges*100:.1f}%)")

# Self-inductance formula test
MU_0 = 4 * np.pi * 1e-7

print(f"\n" + "=" * 70)
print("Self-Inductance Formula Test")
print("=" * 70)

wire_radius = 2e-3  # m (assumed)
print(f"\nAssumed wire radius: {wire_radius*1e3:.1f} mm")

# Test formula for different edge lengths
test_lengths = [0.5e-3, 1.0e-3, 1.5e-3, 2.0e-3, 3.0e-3]
print(f"\nL_self = (mu_0 * l / 2pi) * [ln(l/a) + 0.25]")
print(f"\nwhere a = {wire_radius*1e3:.1f} mm\n")

for l in test_lengths:
    log_term = np.log(l / wire_radius)
    L_self = (MU_0 * l / (2 * np.pi)) * (log_term + 0.25)
    print(f"  l = {l*1e3:.1f} mm: ln(l/a) = {log_term:.3f}, L_self = {L_self*1e9:.2f} nH")

print(f"\n" + "=" * 70)
print("PROBLEM IDENTIFIED")
print("=" * 70)

n_negative = np.sum(edge_lengths < wire_radius)
print(f"\nEdges shorter than wire radius ({wire_radius*1e3:.1f} mm): {n_negative} ({n_negative/n_edges*100:.1f}%)")
print(f"\nFor these edges, ln(l/a) < 0, causing NEGATIVE self-inductance!")

print(f"\n" + "=" * 70)
print("SOLUTION")
print("=" * 70)

# Calculate effective wire radius from mesh
avg_edge_length = edge_lengths.mean()
effective_radius = avg_edge_length / 10  # Rule of thumb: radius ~ 1/10 of mesh size

print(f"\nAverage edge length: {avg_edge_length*1e3:.4f} mm")
print(f"Effective wire radius (l_avg/10): {effective_radius*1e3:.4f} mm")

# Re-test with effective radius
print(f"\nSelf-inductance with effective radius a = {effective_radius*1e3:.4f} mm:\n")

L_self_total = 0
for l in edge_lengths[:10]:  # Test first 10 edges
    log_term = np.log(l / effective_radius)
    L_self = (MU_0 * l / (2 * np.pi)) * (log_term + 0.25)
    L_self_total += L_self
    if L_self < 0:
        print(f"  l = {l*1e3:.4f} mm: L_self = {L_self*1e9:.2f} nH (STILL NEGATIVE!)")

n_negative_new = np.sum(edge_lengths < effective_radius)
print(f"\nEdges shorter than effective radius: {n_negative_new} ({n_negative_new/n_edges*100:.1f}%)")

print(f"\n" + "=" * 70)
print("RECOMMENDED FIX")
print("=" * 70)

print(f"\n1. Use minimum radius constraint:")
print(f"   a_eff = max(l/e, a_min)")
print(f"   where e = 2.71828 (Euler's number)")
print(f"   This ensures ln(l/a) >= -1\n")

print(f"2. Use alternative formula for short segments:")
print(f"   For l < critical_length, use constant approximation\n")

print(f"3. Increase mesh size (current: ~{avg_edge_length*1e3:.2f} mm)")
print(f"   Recommended: > 2mm for 2mm wire radius\n")

gmsh.finalize()
