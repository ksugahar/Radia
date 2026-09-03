"""
PEEC Analysis from 1D Line Mesh (GMSH)

Workflow:
    1. Load 1D line mesh (circular_coil_1d.msh)
    2. Extract edge elements (start_node, end_node pairs)
    3. Convert edges to PEEC segments (center, direction, length)
    4. Build PEEC matrices with specified cross-section
    5. Compute port-to-port impedance

1D Mesh Structure:
    - Edge elements: 2-node line segments along conductor centerline
    - Ports: Nodesets defining voltage application points
    - Cross-section: Specified separately (not in mesh)

Author: Radia Development Team
Date: 2026-02-12
"""

import sys
from pathlib import Path

import numpy as np

GMSH_MODELS = Path(__file__).resolve().parents[1] / "gmsh_models"
sys.path.insert(0, str(GMSH_MODELS))
from gmsh_ascii_reader import read_gmsh_ascii  # noqa: E402

print("=" * 70)
print("PEEC Analysis from 1D Line Mesh")
print("Using Radia PEEC Solver (FastMaxwell-based)")
print("=" * 70)

# Import PEEC module
try:
    from radia.peec_matrices import PEECBuilder
    print("\nUsing C++ PEEC implementation")
except ImportError:
    print("\nERROR: PEEC module not available")
    print("Build with: pwsh -NoProfile -ExecutionPolicy Bypass -File .\\Build.ps1")
    sys.exit(1)

# Check mesh file
mesh_file = GMSH_MODELS / "circular_coil_1d.msh"

if not mesh_file.exists():
    print(f"\nERROR: Mesh file not found: {mesh_file}")
    print("\nGenerate with:")
    print("  cd gmsh_models")
    print("  python generate_1d_coil_mesh.py   # drives Coreform Cubit headless")
    sys.exit(1)

print(f"\n[1] Reading 1D line mesh fixture: {mesh_file.name}")

mesh = read_gmsh_ascii(mesh_file)

node_tags = mesh.node_tags
n_nodes = len(node_tags)
coords = mesh.coords * 1e-3  # mm to m

print(f"    Nodes: {n_nodes}")

# Get edge elements (type 1 = 2-node line)
elem_types, elem_tags, elem_node_tags = mesh.grouped_elements()

edges = []
for i, elem_type in enumerate(elem_types):
    if elem_type == 1:  # Line element (2 nodes)
        node_tags_flat = elem_node_tags[i]
        n_elems = len(node_tags_flat) // 2
        for j in range(n_elems):
            n0 = node_tags_flat[j*2]
            n1 = node_tags_flat[j*2 + 1]
            edges.append((n0, n1))

n_edges = len(edges)
print(f"    Edge elements: {n_edges}")

if n_edges == 0:
    print("\nERROR: No edge elements found!")
    print("Check mesh generation - should create 1D edge mesh")
    sys.exit(1)

# Define port nodes by coordinate-based search
print("\n[2] Defining port nodes (coordinate-based)...")

# Estimate mean radius from node coordinates
r_xy = np.sqrt(coords[:, 0]**2 + coords[:, 1]**2)
r_mean = r_xy.mean()

# Port targets: phi=0 (+X) and phi=180 (-X)
port_positive_target = np.array([r_mean, 0, 0])
port_negative_target = np.array([-r_mean, 0, 0])

# Find closest nodes to port targets
port_positive_node = None
port_negative_node = None
min_dist_pos = float('inf')
min_dist_neg = float('inf')

for i, tag in enumerate(node_tags):
    pos = coords[i]

    # Distance to positive port target
    dist_pos = np.linalg.norm(pos - port_positive_target)
    if dist_pos < min_dist_pos:
        min_dist_pos = dist_pos
        port_positive_node = tag
        port_positive_idx = i

    # Distance to negative port target
    dist_neg = np.linalg.norm(pos - port_negative_target)
    if dist_neg < min_dist_neg:
        min_dist_neg = dist_neg
        port_negative_node = tag
        port_negative_idx = i

if port_positive_node:
    pos = coords[port_positive_idx]
    print(f"    Port positive: node {port_positive_node} at "
          f"({pos[0]:.4f}, {pos[1]:.4f}, {pos[2]:.4f}) [dist={min_dist_pos*1e3:.2f}mm]")

if port_negative_node:
    pos = coords[port_negative_idx]
    print(f"    Port negative: node {port_negative_node} at "
          f"({pos[0]:.4f}, {pos[1]:.4f}, {pos[2]:.4f}) [dist={min_dist_neg*1e3:.2f}mm]")

# Conductor cross-section parameters
# (NOT in mesh - specified separately for PEEC)
width = 4e-3   # 4mm wire width (radial direction)
height = 4e-3  # 4mm wire height (axial direction)
sigma = 5.8e7  # S/m (copper)

print(f"\n[3] Conductor parameters:")
print(f"    Cross-section: {width*1e3:.1f} x {height*1e3:.1f} mm")
print(f"    Area: {width*height*1e6:.2f} mm^2")
print(f"    Conductivity: {sigma:.2e} S/m (copper)")

# Build PEEC model
print(f"\n[4] Building PEEC model...")
print(f"    Adding {n_edges} segments from 1D mesh edges...")

builder = PEECBuilder()

# Convert each edge to PEEC segment
# Edge: (node0, node1) -> Segment: (center, direction, length, width, height)
for n0, n1 in edges:
    # Get node positions
    idx0 = np.where(node_tags == n0)[0][0]
    idx1 = np.where(node_tags == n1)[0][0]
    p0 = coords[idx0]
    p1 = coords[idx1]

    # Create wire segment
    # PEECBuilder.create_wire() takes (start, end, width, height, n_segments=1, sigma)
    builder.create_wire(p0, p1, width, height, 1, sigma)

print(f"    Building matrices...")
L, R, P, M_LS = builder.build()

n_loop = builder.n_loop
n_star = builder.n_star

print(f"\n[5] PEEC matrices built:")
print(f"    L matrix: {L.shape}")
print(f"    R matrix: {len(R)}")
print(f"    P matrix: {P.shape}")
print(f"    M_LS: {M_LS.shape}")
print(f"    Loop DOFs: {n_loop}")
print(f"    Star DOFs: {n_star}")

# DC parameters
print(f"\n[6] DC parameters:")

R_dc = np.sum(R)
print(f"    Total DC resistance: {R_dc*1e3:.4f} mOhm")

# Verify against analytical formula
circumference = 0  # Will compute from mesh
for n0, n1 in edges:
    idx0 = np.where(node_tags == n0)[0][0]
    idx1 = np.where(node_tags == n1)[0][0]
    circumference += np.linalg.norm(coords[idx1] - coords[idx0])

R_dc_analytical = circumference / (sigma * width * height)
error_R = abs(R_dc - R_dc_analytical) / R_dc_analytical * 100

print(f"    Analytical R_dc: {R_dc_analytical*1e3:.4f} mOhm")
print(f"    Error: {error_R:.2f}%")

# Total inductance
L_total = np.sum(L)
print(f"\n[7] Inductance:")
print(f"    PEEC L (sum): {L_total*1e6:.3f} uH")

# Analytical formula - estimate mean radius from coordinates
r_xy = np.sqrt(coords[:, 0]**2 + coords[:, 1]**2)
r_mean = r_xy.mean()

a = np.sqrt(width * height / np.pi)  # Equivalent wire radius
L_analytical = 4 * np.pi * 1e-7 * r_mean * (np.log(8 * r_mean / a) - 2)

print(f"    Analytical L: {L_analytical*1e6:.3f} uH")
print(f"    Ratio (PEEC/Analytical): {L_total/L_analytical:.3f}\"")

error_L = abs(L_total - L_analytical) / L_analytical * 100
print(f"    Error: {error_L:.1f}%")

if error_L < 10:
    print(f"    OK Error < 10% - Excellent!")
elif error_L < 20:
    print(f"    OK Error < 20% - Good")
else:
    print(f"    X Error > 20% - Check mesh quality")

# Port-to-port impedance analysis
print(f"\n[8] Port-to-port impedance analysis:")

# Create port excitation vector
# Find edges connected to port nodes
port_positive_edges = []
port_negative_edges = []

for i, (n0, n1) in enumerate(edges):
    if n0 == port_positive_node or n1 == port_positive_node:
        port_positive_edges.append(i)
    if n0 == port_negative_node or n1 == port_negative_node:
        port_negative_edges.append(i)

print(f"    Port positive: {len(port_positive_edges)} connected edges")
print(f"    Port negative: {len(port_negative_edges)} connected edges")

# Port current excitation: +1 at port_positive, -1 at port_negative
I_port = np.zeros(n_loop)
for idx in port_positive_edges:
    I_port[idx] = 1.0 / len(port_positive_edges)  # Normalized
for idx in port_negative_edges:
    I_port[idx] = -1.0 / len(port_negative_edges)  # Normalized

test_freqs = [1e3, 10e3, 100e3, 1e6]

for f in test_freqs:
    # Port impedance: Z = V_port / I_total
    Z_port = builder.compute_impedance(f, I_port)

    mag = np.abs(Z_port)
    angle = np.angle(Z_port, deg=True)

    print(f"    {f/1e3:8.1f} kHz: |Z| = {mag:.4e} Ohm, angle={angle:.1f}°")

print(f"\n" + "=" * 70)
print("PEEC Analysis from 1D Mesh Completed!")
print("=" * 70)

print(f"\nSummary:")
print(f"  Mesh: {mesh_file.name}")
print(f"  Nodes: {n_nodes}, Edges: {n_edges}")
print(f"  R_dc: {R_dc*1e3:.4f} mOhm (error {error_R:.1f}%)")
print(f"  L: {L_total*1e6:.3f} uH (error {error_L:.1f}%)")
print(f"  Method: 1D line mesh -> PEEC segments")
