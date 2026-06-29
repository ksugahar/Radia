"""
Generate Dual Mesh for PEEC: Filament (1D) + Panel (2D)

Creates both filament elements (centerline) and panel elements (surface)
for FastImp-style Loop-Star PEEC analysis.

Workflow:
    1. Define conductor centerline (1D curve)
    2. Sweep/Extrude to create surface (2D shell)
    3. Mesh centerline with 1D edge elements (Filaments)
    4. Mesh surface with 2D surface elements (Panels)
    5. Export to GMSH with both element types

Author: Radia Development Team
Date: 2026-02-13
"""

import sys
import os
import numpy as np

# Add Cubit to path
# Auto-detect Cubit installation
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
from install_panels import find_cubit_bin as _fcb
_cubit_path = _fcb()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

import cubit
cubit.init(['cubit', '-nojournal', '-batch'])


print("=" * 70)
print("Dual Mesh Generation: Filament (1D) + Panel (2D)")
print("=" * 70)

# Coil parameters (in mm)
mean_radius = 50.0      # Mean radius of coil (mm)
wire_width = 4.0        # Wire width (mm)
wire_height = 4.0       # Wire height (mm)
n_segments_circ = 36    # Segments around circle
n_segments_cross = 8    # Segments around cross-section

print(f"\n[1] Geometry Parameters:")
print(f"    Mean radius: {mean_radius} mm")
print(f"    Wire cross-section: {wire_width} x {wire_height} mm")
print(f"    Circumferential segments: {n_segments_circ}")
print(f"    Cross-section segments: {n_segments_cross}")

# Reset Cubit
cubit.cmd("reset")

# ============================================================================
# Step 1: Create centerline curve (Filament path)
# ============================================================================
print("\n[2] Creating centerline curve (for Filaments)...")

# Create circle in XY plane (z=0)
cubit.cmd(f"create curve arc radius {mean_radius} center location 0 0 0 "
          f"normal 0 0 1 start angle 0 stop angle 360")

centerline_curve_id = cubit.get_last_id("curve")
print(f"    Centerline curve ID: {centerline_curve_id}")

# ============================================================================
# Step 2: Create cross-section profile
# ============================================================================
print("\n[3] Creating rectangular cross-section profile...")

# Create rectangular profile perpendicular to centerline
# Position at (r, 0, 0) with normal in radial direction
profile_x = mean_radius
profile_y = 0
profile_z = 0

# Create rectangle in YZ plane at x=r
cubit.cmd(f"create surface rectangle width {wire_height} height {wire_width} "
          f"xplane offset {profile_x} {profile_y} {profile_z}")

profile_surface_id = cubit.get_last_id("surface")
print(f"    Cross-section surface ID: {profile_surface_id}")

# ============================================================================
# Step 3: Sweep cross-section along centerline to create conductor surface
# ============================================================================
print("\n[4] Sweeping cross-section along centerline...")

cubit.cmd(f"sweep surface {profile_surface_id} along curve {centerline_curve_id}")

conductor_surface_id = cubit.get_last_id("surface")
print(f"    Conductor surface ID: {conductor_surface_id}")

# Note: After sweep, the centerline curve still exists and can be meshed separately

# ============================================================================
# Step 4: Mesh centerline with 1D edge elements (Filaments)
# ============================================================================
print("\n[5] Meshing centerline (1D Filament elements)...")

cubit.cmd(f"curve {centerline_curve_id} interval {n_segments_circ}")
cubit.cmd(f"curve {centerline_curve_id} scheme equal")
cubit.cmd(f"mesh curve {centerline_curve_id}")

n_nodes_1d = cubit.get_node_count()
n_edges_1d = cubit.get_edge_count()

print(f"    Filament mesh: {n_nodes_1d} nodes, {n_edges_1d} edge elements")

# ============================================================================
# Step 5: Mesh surface with 2D elements (Panels)
# ============================================================================
print("\n[6] Meshing conductor surface (2D Panel elements)...")

# Set mesh size for surface
cubit.cmd(f"surface {conductor_surface_id} size auto factor 5")
cubit.cmd(f"surface {conductor_surface_id} scheme pave")  # Quad/tri meshing
cubit.cmd(f"mesh surface {conductor_surface_id}")

n_tris = cubit.get_tri_count()
n_quads = cubit.get_quad_count()

print(f"    Panel mesh: {n_tris} triangles, {n_quads} quads")

# ============================================================================
# Step 6: Define physical groups (blocks)
# ============================================================================
print("\n[7] Defining physical groups...")

# Block 1: Filament elements (1D edges on centerline)
cubit.cmd(f"block 1 add curve {centerline_curve_id}")
cubit.cmd("block 1 name 'filament'")
print(f"    Block 1: 'filament' (1D edge elements)")

# Block 2: Panel elements (2D surface elements)
cubit.cmd(f"block 2 add surface {conductor_surface_id}")
cubit.cmd("block 2 name 'panel'")
print(f"    Block 2: 'panel' (2D surface elements)")

# ============================================================================
# Step 7: Export to GMSH v4.1 format
# ============================================================================
output_file = "dual_mesh_filament_panel.msh"
print(f"\n[8] Exporting to GMSH v4.1 format...")

cubit.cmd(f'export gmsh "{output_file}" overwrite')

print(f"    OK Created: {output_file}")
print(f"    Format: GMSH v4.1 (1D Filament + 2D Panel)")

# ============================================================================
# Verification
# ============================================================================
print("\n[9] Mesh verification:")
print(f"    1D edge elements: {n_edges_1d} (Filaments)")
print(f"    2D surface elements: {n_tris + n_quads} (Panels)")
print(f"    Total nodes: {cubit.get_node_count()}")

print("\n" + "=" * 70)
print("Dual Mesh Summary")
print("=" * 70)
print("\nPhysical Groups:")
print("  1. 'filament' (Block 1): 1D edge elements along centerline")
print("     - Loop DOFs: Current flow")
print("     - Used for inductance matrix (L)")
print("\n  2. 'panel' (Block 2): 2D surface elements on conductor surface")
print("     - Star DOFs: Surface charge distribution")
print("     - Used for potential coefficient matrix (P)")

print("\nLoop-Star Coupling:")
print("  M_LS matrix couples filament currents with panel charges")

print("\n" + "=" * 70)
print("Next Steps")
print("=" * 70)
print("\n1. View dual mesh in GMSH GUI:")
print(f"   gmsh {output_file}")
print("   - Should see both 1D edges (centerline) and 2D surface")
print("\n2. Load into Python and create PEEC model:")
print("   python demo_peec_from_dual_mesh.py")

print("\n" + "=" * 70)
print("OK Dual mesh (Filament + Panel) generated!")
print("=" * 70)
