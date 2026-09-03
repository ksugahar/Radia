"""
PEEC Analysis with Port-Defined GMSH Mesh

Workflow:
    1. Load GMSH mesh with port definitions (circular_coil_with_ports.msh)
    2. Extract port regions (port_positive, port_negative)
    3. Build PEEC model from mesh edges
    4. Create port excitation vector
    5. Compute port-to-port impedance

Port Configuration:
    - port_positive (phi=0, +X): Voltage applied
    - port_negative (phi=180, -X): Ground
    - Z_port = V_applied / I_total (two-port impedance)

FastMaxwell Integration:
    - Mesh edges -> PEEC segments
    - Port regions -> Excitation vector
    - Loop-Star decomposition

Author: Radia Development Team
Date: 2026-02-12
"""

import sys
from pathlib import Path

# Repo-local import path (works before `pip install radia`): put <repo>/src on
# sys.path so `radia.peec_matrices` resolves from the source tree.
_repo = next(p for p in Path(__file__).resolve().parents
             if (p / "src" / "radia").exists())
sys.path.insert(0, str(_repo / "src"))

import numpy as np

GMSH_MODELS = Path(__file__).resolve().parent
sys.path.insert(0, str(GMSH_MODELS))
from gmsh_ascii_reader import read_gmsh_ascii  # noqa: E402

print("=" * 70)
print("PEEC Analysis with Port-Defined GMSH Mesh")
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
mesh_file = GMSH_MODELS / "circular_coil_with_ports.msh"

if not mesh_file.exists():
    print(f"\nERROR: Mesh file not found: {mesh_file}")
    print("\nGenerate with:")
    print("  cd gmsh_models")
    print("  python generate_coil_with_ports.py   # drives Coreform Cubit headless")
    sys.exit(1)

print(f"\n[1] Reading GMSH mesh fixture with ports: {mesh_file.name}")

mesh = read_gmsh_ascii(mesh_file)

node_tags = mesh.node_tags
n_nodes = len(node_tags)
coords = mesh.coords * 1e-3  # mm to m

print(f"    Nodes: {n_nodes}")

# Get triangles
elem_types, elem_tags, elem_node_tags = mesh.grouped_elements()

triangles = []
for i, elem_type in enumerate(elem_types):
    if elem_type == 2:  # Triangle
        node_tags_flat = elem_node_tags[i]
        n_elems = len(node_tags_flat) // 3
        for j in range(n_elems):
            nodes = node_tags_flat[j*3:(j+1)*3]
            triangles.append(tuple(nodes))

print(f"    Triangles: {len(triangles)}")

# Define port regions using coordinate-based selection
# (Sideset not exported by export gmsh, use coordinates instead)
print("\n[2] Defining port regions (coordinate-based)...")

# Triangle ID -> nodes mapping
tri_id_to_nodes = {}
for i, elem_type in enumerate(elem_types):
    if elem_type == 2:  # Triangle
        node_tags_flat = elem_node_tags[i]
        tri_tags_flat = elem_tags[i]
        n_elems = len(node_tags_flat) // 3
        for j in range(n_elems):
            tri_id = tri_tags_flat[j]
            nodes = tuple(node_tags_flat[j*3:(j+1)*3])
            tri_id_to_nodes[tri_id] = nodes

# Select port triangles by centroid position
port_positive_tris = set()
port_negative_tris = set()

for tri_id, nodes in tri_id_to_nodes.items():
    # Compute triangle centroid
    n0, n1, n2 = nodes
    idx0 = np.where(node_tags == n0)[0][0]
    idx1 = np.where(node_tags == n1)[0][0]
    idx2 = np.where(node_tags == n2)[0][0]

    centroid = (coords[idx0] + coords[idx1] + coords[idx2]) / 3.0

    # Port positive: phi=0 (+X), centroid x > 45mm
    if centroid[0] > 0.045:  # 45mm in meters
        port_positive_tris.add(tri_id)

    # Port negative: phi=180 (-X), centroid x < -45mm
    elif centroid[0] < -0.045:
        port_negative_tris.add(tri_id)

print(f"    Port positive triangles: {len(port_positive_tris)} (x > 45mm)")
print(f"    Port negative triangles: {len(port_negative_tris)} (x < -45mm)")

if not port_positive_tris or not port_negative_tris:
    print("\nERROR: Port regions not found!")
    print("Check coordinate thresholds.")
    sys.exit(1)

# Extract edges and identify port edges
print("\n[3] Extracting edges...")

# Extract unique edges and track which triangles they belong to
edge_to_tris = {}  # edge -> set of triangle IDs
for tri_id, nodes in tri_id_to_nodes.items():
    n0, n1, n2 = nodes
    edges_tri = [
        tuple(sorted([n0, n1])),
        tuple(sorted([n1, n2])),
        tuple(sorted([n2, n0]))
    ]
    for edge in edges_tri:
        if edge not in edge_to_tris:
            edge_to_tris[edge] = set()
        edge_to_tris[edge].add(tri_id)

# Extract BOUNDARY edges only (connected to exactly 1 triangle)
# Interior edges (connected to 2 triangles) are excluded
# For surface mesh of a conductor, boundary edges represent the conductor outline
boundary_edges = []
interior_edges = []

for edge, tri_ids in edge_to_tris.items():
    if len(tri_ids) == 1:
        boundary_edges.append(edge)
    else:
        interior_edges.append(edge)

print(f"    Total edges: {len(edge_to_tris)}")
print(f"    Boundary edges: {len(boundary_edges)} (used for PEEC)")
print(f"    Interior edges: {len(interior_edges)} (excluded)")

# Use only boundary edges for PEEC
edges = boundary_edges
n_edges = len(edges)

# Classify edges by port membership
port_positive_edges = set()
port_negative_edges = set()

for edge, tri_ids in edge_to_tris.items():
    # Check if edge belongs to port_positive (any connected triangle in port_positive)
    if tri_ids & port_positive_tris:
        port_positive_edges.add(edge)

    # Check if edge belongs to port_negative
    if tri_ids & port_negative_tris:
        port_negative_edges.add(edge)

print(f"    Port positive edges: {len(port_positive_edges)}")
print(f"    Port negative edges: {len(port_negative_edges)}")

# Estimate conductor parameters
edge_lengths = []
for n0, n1 in list(edges)[:1000]:  # Sample
    idx0 = np.where(node_tags == n0)[0][0]
    idx1 = np.where(node_tags == n1)[0][0]
    p0 = coords[idx0]
    p1 = coords[idx1]
    edge_lengths.append(np.linalg.norm(p1 - p0))

avg_edge_length = np.mean(edge_lengths)
effective_width = 0.277e-3  # From previous analysis
width = effective_width
height = effective_width
sigma = 5.8e7  # Copper

print(f"\n[4] Conductor parameters:")
print(f"    Average edge length: {avg_edge_length*1e3:.3f} mm")
print(f"    Cross-section: {width*1e3:.3f} x {height*1e3:.3f} mm")
print(f"    Conductivity: {sigma:.2e} S/m (copper)")

# Build PEEC model
print(f"\n[5] Building PEEC model...")
print(f"    Adding {n_edges} segments from mesh edges...")

builder = PEECBuilder()

# Create edge ID to index mapping
edge_to_idx = {}
for idx, edge in enumerate(edges):
    # Get node positions
    n0, n1 = edge
    idx0 = np.where(node_tags == n0)[0][0]
    idx1 = np.where(node_tags == n1)[0][0]
    p0 = coords[idx0]
    p1 = coords[idx1]

    # Create wire segment
    builder.create_wire(p0, p1, width, height, 1, sigma)
    edge_to_idx[edge] = idx

print(f"    Building matrices...")
L, R, P, M_LS = builder.build()

n_loop = builder.n_loop
n_star = builder.n_star

print(f"\n[6] PEEC matrices built:")
print(f"    L matrix: {L.shape}")
print(f"    R matrix: {len(R)}")
print(f"    P matrix: {P.shape}")
print(f"    M_LS: {M_LS.shape}")
print(f"    Loop DOFs: {n_loop}")
print(f"    Star DOFs: {n_star}")

# Create port excitation vector
print(f"\n[7] Creating port excitation vector...")

# Port current vector: +1 for port_positive edges, -1 for port_negative edges
I_port = np.zeros(n_loop)

n_pos = 0
n_neg = 0
for edge, idx in edge_to_idx.items():
    if edge in port_positive_edges:
        I_port[idx] = 1.0
        n_pos += 1
    elif edge in port_negative_edges:
        I_port[idx] = -1.0
        n_neg += 1

# Normalize by number of port edges
if n_pos > 0:
    for edge, idx in edge_to_idx.items():
        if edge in port_positive_edges:
            I_port[idx] = 1.0 / n_pos

if n_neg > 0:
    for edge, idx in edge_to_idx.items():
        if edge in port_negative_edges:
            I_port[idx] = -1.0 / n_neg

print(f"    Port excitation created:")
print(f"      Positive port: {n_pos} edges (normalized)")
print(f"      Negative port: {n_neg} edges (normalized)")

# Frequency response with port impedance
print(f"\n[8] Computing port impedance:")

test_frequencies = [1e3, 10e3, 100e3, 1e6]  # Hz

for f in test_frequencies:
    # Use PEECBuilder's built-in impedance calculation
    Z_port = builder.compute_impedance(f, I_port)

    mag = np.abs(Z_port)
    angle = np.angle(Z_port, deg=True)

    print(f"    {f/1e3:8.1f} kHz: |Z| = {mag:.4e} Ohm, angle = {angle:.1f} deg")

# DC parameters
print(f"\n[9] DC parameters:")

R_dc = np.sum(R)
L_total = np.sum(L)

print(f"    Total DC resistance: {R_dc*1e3:.4f} mOhm")
print(f"    Total inductance (sum): {L_total*1e6:.4f} uH")

# Analytical comparison
r_xy = np.sqrt(coords[:, 0]**2 + coords[:, 1]**2)
r_mean = r_xy.mean()

a = np.sqrt(width * height / np.pi)
L_analytical = 4*np.pi*1e-7 * r_mean * (np.log(8*r_mean/a) - 2)

print(f"\n[10] Analytical comparison:")
print(f"    Mean radius: {r_mean*1e3:.1f} mm")
print(f"    Analytical inductance: {L_analytical*1e6:.4f} uH")
print(f"    PEEC inductance: {L_total*1e6:.4f} uH")
print(f"    Ratio (PEEC/Analytical): {L_total/L_analytical:.3f}")

error = abs(L_total - L_analytical) / L_analytical * 100
print(f"    Error: {error:.1f}%")

print(f"\n" + "=" * 70)
print("PEEC Port Analysis Completed!")
print("=" * 70)

print(f"\nSummary:")
print(f"  Mesh: {mesh_file.name}")
print(f"  Nodes: {n_nodes}, Edges: {n_edges}")
print(f"  Port positive: {len(port_positive_edges)} edges")
print(f"  Port negative: {len(port_negative_edges)} edges")
print(f"  L_PEEC: {L_total*1e6:.4f} uH")
print(f"  L_analytical: {L_analytical*1e6:.4f} uH")
print(f"  Error: {error:.1f}%")
