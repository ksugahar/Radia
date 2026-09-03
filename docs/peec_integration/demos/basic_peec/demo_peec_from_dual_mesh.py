"""
PEEC Analysis from Dual Mesh (Filament + Panel)

Loads GMSH mesh containing:
    - Block 1: 1D edge elements (Filaments) for Loop DOFs
    - Block 2: 2D surface elements (Panels) for Star DOFs

Builds complete Loop-Star PEEC matrices:
    - L matrix: Filament-Filament inductance
    - P matrix: Panel-Panel potential coefficient
    - M_LS matrix: Filament-Panel coupling
    - R matrix: Filament resistance

Author: Radia Development Team
Date: 2026-02-13
"""

import sys
from pathlib import Path

# Repo-local import path (works before `pip install radia`): put <repo>/src on
# sys.path so `radia.peec_matrices` resolves from the source tree.
_repo = next(p for p in Path(__file__).resolve().parents
             if (p / "src" / "radia").exists())
sys.path.insert(0, str(_repo / "src"))

import numpy as np

GMSH_MODELS = Path(__file__).resolve().parents[1] / "gmsh_models"
sys.path.insert(0, str(GMSH_MODELS))
from gmsh_ascii_reader import read_gmsh_ascii  # noqa: E402

print("=" * 70)
print("PEEC Analysis from Dual Mesh (Filament + Panel)")
print("Using FastImp-style Loop-Star Decomposition")
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
mesh_file = GMSH_MODELS / "dual_mesh_filament_panel.msh"

if not mesh_file.exists():
    print(f"\nERROR: Mesh file not found: {mesh_file}")
    print("\nGenerate with:")
    print("  cd gmsh_models")
    print("  python generate_dual_mesh_filament_panel.py"
          "   # drives Coreform Cubit headless")
    sys.exit(1)

print(f"\n[1] Reading dual mesh fixture: {mesh_file.name}")

mesh = read_gmsh_ascii(mesh_file)

node_tags = mesh.node_tags
n_nodes = len(node_tags)
coords = mesh.coords * 1e-3  # mm to m

print(f"    Nodes: {n_nodes}")

# Get elements by type
elem_types, elem_tags, elem_node_tags = mesh.grouped_elements()

# Extract Filament elements (1D edges, type 1)
filament_edges = []
for i, elem_type in enumerate(elem_types):
    if elem_type == 1:  # Line element (2 nodes)
        node_tags_flat = elem_node_tags[i]
        n_elems = len(node_tags_flat) // 2
        for j in range(n_elems):
            n0 = node_tags_flat[j*2]
            n1 = node_tags_flat[j*2 + 1]
            filament_edges.append((n0, n1))

# Extract Panel elements (2D surface, type 2=tri, type 3=quad)
panel_triangles = []
panel_quads = []

for i, elem_type in enumerate(elem_types):
    if elem_type == 2:  # Triangle (3 nodes)
        node_tags_flat = elem_node_tags[i]
        n_elems = len(node_tags_flat) // 3
        for j in range(n_elems):
            nodes = [node_tags_flat[j*3 + k] for k in range(3)]
            panel_triangles.append(nodes)

    elif elem_type == 3:  # Quad (4 nodes)
        node_tags_flat = elem_node_tags[i]
        n_elems = len(node_tags_flat) // 4
        for j in range(n_elems):
            nodes = [node_tags_flat[j*4 + k] for k in range(4)]
            panel_quads.append(nodes)

n_filaments = len(filament_edges)
n_panels = len(panel_triangles) + len(panel_quads)

print(f"    Filament elements (1D): {n_filaments}")
print(f"    Panel elements (2D): {n_panels} ({len(panel_triangles)} tri + {len(panel_quads)} quad)")

if n_filaments == 0:
    print("\nERROR: No filament elements found!")
    print("Check mesh - Block 1 should contain 1D edge elements")
    sys.exit(1)

if n_panels == 0:
    print("\nERROR: No panel elements found!")
    print("Check mesh - Block 2 should contain 2D surface elements")
    sys.exit(1)

# Conductor parameters
width = 4e-3   # 4mm wire width
height = 4e-3  # 4mm wire height
sigma = 5.8e7  # S/m (copper)

print(f"\n[2] Conductor parameters:")
print(f"    Cross-section: {width*1e3:.1f} x {height*1e3:.1f} mm")
print(f"    Conductivity: {sigma:.2e} S/m (copper)")

# ============================================================================
# Build PEEC model with Filaments and Panels
# ============================================================================

print(f"\n[3] Building PEEC model...")
print(f"    Adding {n_filaments} Filament segments...")

builder = PEECBuilder()

# Add Filament segments (Loop elements)
for n0, n1 in filament_edges:
    idx0 = np.where(node_tags == n0)[0][0]
    idx1 = np.where(node_tags == n1)[0][0]
    p0 = coords[idx0]
    p1 = coords[idx1]

    builder.create_wire(p0, p1, width, height, 1, sigma)

print(f"    Adding {n_panels} Panel elements (Star elements)...")

# Add Panel elements (Star elements) - TRUE 2D SURFACE INTEGRATION
# Convert node indices to vertex coordinates and add panels

# Add triangular panels
for node_indices in panel_triangles:
    vertices = []
    for node_idx in node_indices:
        # Find coordinate index (GMSH node tags start at 1)
        coord_idx = np.where(node_tags == node_idx)[0][0]
        vertices.append(coords[coord_idx].tolist())
    builder.add_panel(vertices)

# Add quadrilateral panels
for node_indices in panel_quads:
    vertices = []
    for node_idx in node_indices:
        coord_idx = np.where(node_tags == node_idx)[0][0]
        vertices.append(coords[coord_idx].tolist())
    builder.add_panel(vertices)

print(f"    Added {len(panel_triangles)} triangle panels + {len(panel_quads)} quad panels")
print(f"    Using TRUE 2D analytical integration (Wilton + Gauss quadrature)")

# Build matrices
print(f"\n[4] Building PEEC matrices...")
L, R, P, M_LS = builder.build()

n_loop = builder.n_loop
n_star = builder.n_star

print(f"\n[5] PEEC matrices built:")
print(f"    L matrix: {L.shape} (Filament-Filament inductance)")
print(f"    R matrix: {len(R)} (Filament resistance)")
print(f"    P matrix: {P.shape} (Panel-Panel potential coefficient)")
print(f"    M_LS: {M_LS.shape} (Filament-Panel coupling)")
print(f"    Loop DOFs (Filaments): {n_loop}")
print(f"    Star DOFs (Panels): {n_star}")

# DC parameters
print(f"\n[6] DC parameters:")

R_dc = np.sum(R)
print(f"    Total DC resistance: {R_dc*1e3:.4f} mOhm")

# Total inductance
L_total = np.sum(L)
print(f"\n[7] Inductance:")
print(f"    PEEC L (sum): {L_total*1e6:.3f} uH")

# Analytical formula - circular loop
r_xy = np.sqrt(coords[:, 0]**2 + coords[:, 1]**2)
r_mean = r_xy.mean()

a = np.sqrt(width * height / np.pi)  # Equivalent wire radius
L_analytical = 4 * np.pi * 1e-7 * r_mean * (np.log(8 * r_mean / a) - 2)

print(f"    Analytical L: {L_analytical*1e6:.3f} uH")
print(f"    Ratio (PEEC/Analytical): {L_total/L_analytical:.3f}")

error_L = abs(L_total - L_analytical) / L_analytical * 100
print(f"    Error: {error_L:.1f}%")

if error_L < 2:
    print(f"    OK Error < 2% - Excellent!")
elif error_L < 10:
    print(f"    OK Error < 10% - Good")
else:
    print(f"    X Error > 10% - Check implementation")

print(f"\n" + "=" * 70)
print("PEEC Analysis from Dual Mesh Completed!")
print("=" * 70)

print(f"\nSummary:")
print(f"  Mesh: {mesh_file.name}")
print(f"  Filaments (Loop): {n_filaments} elements")
print(f"  Panels (Star): {n_panels} elements (TRUE 2D analytical integration)")
print(f"  L: {L_total*1e6:.3f} uH (error {error_L:.1f}%)")
print(f"\nPanel Integration:")
print(f"  - Self-potential: Wilton analytical formula")
print(f"  - Near-field (d < 3*sqrt(A)): 3-point Gauss quadrature")
print(f"  - Far-field (d > 3*sqrt(A)): Centroid approximation")
