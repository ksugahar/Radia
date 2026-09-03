"""
Generate 1D Line Mesh for PEEC Analysis using Coreform Cubit

Creates a circular coil as 1D curve mesh:
    - Centerline curve only (no volume, no surface)
    - 1D edge elements along the curve
    - Export to GMSH format for Python loading

PEEC Requirements:
    - 1D segments (centerline of conductor)
    - Rectangular cross-section specified separately (width, height)
    - Ports defined by node positions

Author: Radia Development Team
Date: 2026-02-12
"""

import sys

import numpy as np

# Add Cubit to path (auto-detect via radia's install_panels helper)
from radia.install_panels import find_cubit_bin as _fcb
_cubit_path = _fcb()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

import cubit
cubit.init(['cubit', '-nojournal', '-batch', '-nographics'])


print("=" * 70)
print("1D Circular Coil Mesh for PEEC (Cubit -> GMSH)")
print("=" * 70)

# Coil parameters (in mm)
mean_radius = 50.0      # Mean radius of coil (mm)
n_segments = 36         # Number of segments along circle

print(f"\n[1] Coil Parameters:")
print(f"    Mean radius: {mean_radius} mm")
print(f"    Segments: {n_segments}")

# Reset Cubit
cubit.cmd("reset")

# Create circular curve (centerline only)
print("\n[2] Creating circular centerline curve...")

# Create circle in XY plane (z=0)
cubit.cmd(f"create curve arc radius {mean_radius} center location 0 0 0 "
          f"normal 0 0 1 start angle 0 stop angle 360")

# Get curve ID
curve_id = cubit.get_last_id("curve")
print(f"    Created curve {curve_id}")

# Set mesh size
edge_length = 2 * np.pi * mean_radius / n_segments
cubit.cmd(f"curve {curve_id} interval {n_segments}")
cubit.cmd(f"curve {curve_id} scheme equal")

# Generate 1D mesh (edge elements only)
print("\n[3] Generating 1D edge mesh...")
cubit.cmd(f"mesh curve {curve_id}")

# Get mesh info
n_nodes = cubit.get_node_count()
n_edges = cubit.get_edge_count()

print(f"    Nodes: {n_nodes}")
print(f"    Edge elements: {n_edges}")

# Define physical group (block for conductor)
print("\n[4] Defining physical groups...")
cubit.cmd(f"block 1 add curve {curve_id}")
cubit.cmd("block 1 name 'conductor'")
# Note: element type specification not needed for 1D edge elements
print("    Block 1: 'conductor' (1D edge elements)")
print("    NOTE: Ports will be defined by coordinate-based search in Python")
print("          (nodeset/sideset not supported by export gmsh)")

# Export to GMSH v4.1 format
output_file = "circular_coil_1d.msh"
print(f"\n[6] Exporting to GMSH v4.1 format...")
cubit.cmd(f'export gmsh "{output_file}" overwrite')
print(f"    OK Created: {output_file}")
print(f"    Format: GMSH v4.1 (1D edge elements with port nodes)")

# Verify: should have NO surface or volume elements
if cubit.get_tri_count() > 0 or cubit.get_quad_count() > 0:
    print("    WARNING: Surface elements found - PEEC only needs 1D edges!")
if cubit.get_tet_count() > 0 or cubit.get_hex_count() > 0:
    print("    WARNING: Volume elements found - PEEC only needs 1D edges!")
else:
    print("    OK No surface/volume elements - correct for PEEC")

print("\n" + "=" * 70)
print("1D Mesh Summary")
print("=" * 70)
print("\nPhysical Groups:")
print("  1. 'conductor' (Block 1): 1D edge elements along centerline")
print("\nPort Definition Strategy:")
print("  - Ports defined by coordinate-based search in Python script")
print("  - Port positive: Node closest to (+R, 0, 0)")
print("  - Port negative: Node closest to (-R, 0, 0)")
print("  - No need for nodeset/sideset in Cubit mesh")
print("\nPEEC Segment Generation:")
print("  - Each edge element becomes a PEEC segment")
print("  - Segment: (start_node, end_node) -> (center, direction, length)")
print("  - Cross-section: specified separately (width, height)")

print("\n" + "=" * 70)
print("Next Steps")
print("=" * 70)
print("\n1. View 1D mesh in GMSH GUI:")
print(f"   gmsh {output_file}")
print("   - Should see curve with nodes only (no surface)")
print("\n2. Load into Python and create PEEC model:")
print("   cd ..")
print("   python demo_peec_from_1d_mesh.py")

print("\n" + "=" * 70)
print("OK 1D coil mesh for PEEC generated!")
print("=" * 70)
