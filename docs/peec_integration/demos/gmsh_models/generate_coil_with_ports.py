"""
Generate Circular Coil Mesh with PORT DEFINITIONS using Coreform Cubit

Creates a toroidal coil with:
    - Surface mesh for PEEC analysis
    - PORT 1: phi=0 (positive terminal)
    - PORT 2: phi=180 (negative terminal)

Output:
    circular_coil_with_ports.msh - Surface mesh with port physical groups

Requirements:
    - Coreform Cubit 2025.8+
    - cubit_mesh_export module (via radia package)
"""

import sys
from pathlib import Path

import numpy as np

# Add Cubit to path (auto-detect via radia's install_panels helper)
_repo = next(p for p in Path(__file__).resolve().parents
             if (p / "src" / "radia").exists())
sys.path.insert(0, str(_repo / "src" / "radia"))
from install_panels import find_cubit_bin as _fcb
_cubit_path = _fcb()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

import cubit
cubit.init(['cubit', '-nojournal', '-batch', '-nographics'])


print("=" * 70)
print("Circular Coil with PORT Definitions (Cubit -> GMSH)")
print("=" * 70)

# Coil parameters (in mm)
mean_radius = 50.0      # Mean radius of coil (mm)
wire_width = 4.0        # Wire width (radial direction) (mm)
wire_height = 4.0       # Wire height (axial direction) (mm)
mesh_size = 1.0         # Target element size (mm)

print(f"\n[1] Coil Parameters:")
print(f"    Mean radius: {mean_radius} mm")
print(f"    Wire width (radial): {wire_width} mm")
print(f"    Wire height (axial): {wire_height} mm")
print(f"    Target mesh size: {mesh_size} mm")

# Reset Cubit
cubit.cmd("reset")

# Create toroidal coil
print("\n[2] Creating toroidal coil geometry...")
major_radius = mean_radius
minor_radius_x = wire_width / 2.0
minor_radius_z = wire_height / 2.0

cubit.cmd(f"create torus major radius {major_radius} minor radius {minor_radius_x}")
cubit.cmd(f"volume 1 scale z {minor_radius_z / minor_radius_x}")

cubit.cmd("volume all scheme tetmesh")
cubit.cmd(f"volume all size {mesh_size}")

print("    Geometry created: Toroidal coil with rectangular cross-section")

# Generate SURFACE mesh only
print("\n[3] Generating surface mesh (PEEC requirement: surface only)...")
cubit.cmd("mesh surface all")

nodes = cubit.get_node_count()
tris = cubit.get_tri_count()
print(f"    Nodes: {nodes}")
print(f"    Surface elements (triangles): {tris}")

# Define physical groups
print("\n[4] Defining physical groups (conductor + ports)...")

# Main conductor block
cubit.cmd("block 1 add tri all")
cubit.cmd("block 1 name 'conductor'")
cubit.cmd("block 1 element type tri3")
print("    Block 1: 'conductor' (all surface triangles)")

# PORT regions: Use Python API to select triangles by centroid location
print("\n[5] Defining PORT regions (coordinate-based selection)...")

# Get all triangle IDs
all_tris = list(range(1, cubit.get_tri_count() + 1))

# PORT 1: Select triangles with centroid x > 45mm (phi ~ 0)
port1_tris = []
for tri_id in all_tris:
    # Get triangle nodes
    nodes = cubit.get_connectivity("tri", tri_id)
    # Get node coordinates
    coords = np.array([cubit.get_nodal_coordinates(n) for n in nodes])
    centroid = coords.mean(axis=0)

    # Check if centroid is in port region (x > 45)
    if centroid[0] > 45:
        port1_tris.append(tri_id)

if port1_tris:
    # Create sideset for port 1
    tri_list = " ".join(map(str, port1_tris))
    cubit.cmd(f"sideset 1 add tri {tri_list}")
    cubit.cmd("sideset 1 name 'port_positive'")
    print(f"    Sideset 1: 'port_positive' (phi=0, +X) - {len(port1_tris)} triangles")
else:
    print("    WARNING: No triangles found for port_positive!")

# PORT 2: Select triangles with centroid x < -45mm (phi ~ 180)
port2_tris = []
for tri_id in all_tris:
    nodes = cubit.get_connectivity("tri", tri_id)
    coords = np.array([cubit.get_nodal_coordinates(n) for n in nodes])
    centroid = coords.mean(axis=0)

    if centroid[0] < -45:
        port2_tris.append(tri_id)

if port2_tris:
    tri_list = " ".join(map(str, port2_tris))
    cubit.cmd(f"sideset 2 add tri {tri_list}")
    cubit.cmd("sideset 2 name 'port_negative'")
    print(f"    Sideset 2: 'port_negative' (phi=180, -X) - {len(port2_tris)} triangles")
else:
    print("    WARNING: No triangles found for port_negative!")

# Export to GMSH v4.1 format
output_file = "circular_coil_with_ports.msh"
print(f"\n[6] Exporting to GMSH v4.1 format...")
cubit.cmd(f'export gmsh "{output_file}" overwrite')
print(f"    OK Created: {output_file}")
print(f"    Format: GMSH v4.1 (with port physical groups)")

# Verify: should have NO volume elements
if cubit.get_tet_count() > 0 or cubit.get_hex_count() > 0:
    print("    WARNING: Volume elements found - PEEC only needs surface!")
else:
    print("    OK No volume elements - correct for PEEC")

print("\n" + "=" * 70)
print("Port Definition Summary")
print("=" * 70)
print("\nPhysical Groups:")
print("  1. 'conductor' (Block 1): All surface triangles")
print("  2. 'port_positive' (Sideset 1): phi=0 region (+X)")
print("  3. 'port_negative' (Sideset 2): phi=180 region (-X)")
print("\nPEEC Port Configuration:")
print("  - Voltage source: Between port_positive and port_negative")
print("  - Current flow: +X -> toroidal direction -> -X")
print("  - Impedance extraction: Z = V_applied / I_total")

print("\n" + "=" * 70)
print("Next Steps")
print("=" * 70)
print("\n1. View ports in GMSH GUI:")
print(f"   gmsh {output_file}")
print("   - Tools -> Visibility -> Physical groups")
print("   - Check 'port_positive' and 'port_negative'")
print("\n2. Load into Python and extract port elements:")
print("   cd ..")
print("   python demo_gmsh_to_peec.py")

print("\n" + "=" * 70)
print("OK Coil mesh with ports generated!")
print("=" * 70)
