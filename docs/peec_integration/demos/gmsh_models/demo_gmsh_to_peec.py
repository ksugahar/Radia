"""
PEEC Analysis from GMSH Mesh

Workflow:
    1. Load GMSH mesh (circular_coil.msh)
    2. Extract triangular edges
    3. Convert edges to PEEC segments (FastMaxwell approach)
    4. Build PEEC matrices and compute inductance

FastMaxwell Integration:
    - Each mesh edge becomes a PEEC segment
    - Uses internal CreateWireSegments() approach:
      edge (start, end) -> segment (center, direction, length)
    - Loop-Star decomposition for AC analysis
    - DC analysis first (resistance + inductance verification)

Author: Radia Development Team
Date: 2026-02-12
"""

import sys
from pathlib import Path

import numpy as np

GMSH_MODELS = Path(__file__).resolve().parent
sys.path.insert(0, str(GMSH_MODELS))
from gmsh_ascii_reader import read_gmsh_ascii  # noqa: E402

print("=" * 70)
print("PEEC Analysis from GMSH Mesh")
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
mesh_file = GMSH_MODELS / "circular_coil.msh"

if not mesh_file.exists():
    print(f"\nERROR: Mesh file not found: {mesh_file}")
    print("\nGenerate with:")
    print("  cd gmsh_models")
    print("  python generate_coil_cubit.py   # drives Coreform Cubit headless")
    sys.exit(1)

print(f"\n[1] Reading GMSH mesh fixture: {mesh_file.name}")

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

# Extract unique edges
print("\n[2] Extracting edges...")
edges_set = set()
for tri in triangles:
    n0, n1, n2 = tri
    edges_set.add(tuple(sorted([n0, n1])))
    edges_set.add(tuple(sorted([n1, n2])))
    edges_set.add(tuple(sorted([n2, n0])))

edges = list(edges_set)
n_edges = len(edges)

print(f"    Unique edges: {n_edges}")

# Estimate conductor parameters from mesh
edge_lengths = []
for n0, n1 in edges[:1000]:  # Sample
    idx0 = np.where(node_tags == n0)[0][0]
    idx1 = np.where(node_tags == n1)[0][0]
    p0 = coords[idx0]
    p1 = coords[idx1]
    edge_lengths.append(np.linalg.norm(p1 - p0))

avg_edge_length = np.mean(edge_lengths)

print(f"\n[3] Mesh properties:")
print(f"    Average edge length: {avg_edge_length*1e3:.3f} mm")

# Estimate wire cross-section from coil geometry
# Total surface area
total_area = sum(
    0.5 * np.linalg.norm(np.cross(
        coords[np.where(node_tags == tri[1])[0][0]] - coords[np.where(node_tags == tri[0])[0][0]],
        coords[np.where(node_tags == tri[2])[0][0]] - coords[np.where(node_tags == tri[0])[0][0]]
    ))
    for tri in triangles[:100]  # Sample
) * (len(triangles) / 100)

# Estimate perimeter (sum of edge lengths)
perimeter = sum(edge_lengths) * (n_edges / len(edge_lengths))

# Effective wire width
effective_width = total_area / perimeter

print(f"    Total surface area: {total_area*1e6:.2f} mm^2")
print(f"    Estimated perimeter: {perimeter*1e3:.2f} mm")
print(f"    Effective width: {effective_width*1e3:.3f} mm")

# Conductor parameters
width = effective_width
height = effective_width  # Assume square cross-section
sigma = 5.8e7  # Copper

print(f"\n[4] Conductor parameters:")
print(f"    Cross-section: {width*1e3:.3f} x {height*1e3:.3f} mm")
print(f"    Conductivity: {sigma:.2e} S/m (copper)")

# Build PEEC model
print(f"\n[5] Building PEEC model...")
print(f"    Adding {n_edges} segments from mesh edges...")

builder = PEECBuilder()

# Convert each edge to PEEC segment
# Using FastMaxwell approach: edge (start, end) -> segment (center, direction, length)
for n0, n1 in edges:
    # Get node positions
    idx0 = np.where(node_tags == n0)[0][0]
    idx1 = np.where(node_tags == n1)[0][0]
    p0 = coords[idx0]
    p1 = coords[idx1]

    # Create wire segment using create_wire (n_segments=1 for single edge)
    # This internally calls CreateWireSegments() which converts start/end to center/direction/length
    builder.create_wire(p0, p1, width, height, 1, sigma)

print(f"    Building matrices...")
L, R, P, M_LS = builder.build()

n_loop = builder.n_loop
n_star = builder.n_star

print(f"\n[6] PEEC matrices built:")
print(f"    L matrix: {L.shape} (loop inductances)")
print(f"    R matrix: {len(R)} (loop resistances)")
print(f"    P matrix: {P.shape} (star capacitances)")
print(f"    M_LS: {M_LS.shape} (loop-star coupling)")
print(f"    Loop DOFs: {n_loop}")
print(f"    Star DOFs: {n_star}")

# Physical parameters
print(f"\n[7] DC parameters:")

R_dc = np.sum(R)
print(f"    Total DC resistance: {R_dc*1e3:.4f} mOhm")

# Total inductance (sum of L matrix)
L_total = np.sum(L)
print(f"    Total inductance (sum): {L_total*1e9:.2f} nH")

# Analytical comparison
print(f"\n[8] Analytical comparison:")

# Estimate mean radius from coordinates
r_xy = np.sqrt(coords[:, 0]**2 + coords[:, 1]**2)
r_mean = r_xy.mean()

# Analytical formula for circular loop: L = mu_0 * R * (ln(8R/a) - 2)
a = np.sqrt(width * height / np.pi)  # Equivalent wire radius
L_analytical = 4*np.pi*1e-7 * r_mean * (np.log(8*r_mean/a) - 2)

print(f"    Mean radius: {r_mean*1e3:.1f} mm")
print(f"    Equivalent wire radius: {a*1e3:.3f} mm")
print(f"    Analytical inductance: {L_analytical*1e9:.2f} nH")
print(f"    PEEC inductance: {L_total*1e9:.2f} nH")
print(f"    Ratio (PEEC/Analytical): {L_total/L_analytical:.3f}")

error = abs(L_total - L_analytical) / L_analytical * 100
print(f"    Error: {error:.1f}%")

if error < 10:
    print(f"    OK Error < 10% - PEEC solver validated!")
elif error < 20:
    print(f"    WARNING Error 10-20% - Acceptable for mesh approximation")
else:
    print(f"    ERROR Error > 20% - Check mesh quality")

# Frequency response
print(f"\n[9] Frequency response:")

test_frequencies = [1e3, 10e3, 100e3, 1e6, 10e6]  # Hz

for f in test_frequencies:
    omega = 2 * np.pi * f

    # Impedance: Z = R + j*omega*L (simplified, no capacitance)
    Z_LL = np.diag(R) + 1j * omega * np.diag(np.diag(L))

    # Port impedance (all edges in series)
    Z_port = np.sum(Z_LL.diagonal())

    mag = np.abs(Z_port)
    angle = np.angle(Z_port, deg=True)

    print(f"    {f/1e3:8.1f} kHz: |Z| = {mag:.4e} Ohm, angle = {angle:.1f} deg")

print(f"\n" + "=" * 70)
print("PEEC Analysis from GMSH Completed!")
print("=" * 70)

print(f"\nSummary:")
print(f"  Mesh: {mesh_file.name}")
print(f"  Nodes: {n_nodes}, Edges: {n_edges}")
print(f"  L_PEEC: {L_total*1e9:.2f} nH")
print(f"  L_analytical: {L_analytical*1e9:.2f} nH")
print(f"  Error: {error:.1f}%")
print(f"\nMethod: FastMaxwell Loop-Star decomposition")
print(f"  - Mesh edges -> PEEC segments")
print(f"  - Darwin approximation G(r) = 1/(4*pi*r)")
