"""
Generate Circular Coil Mesh using Coreform Cubit

Creates a toroidal coil (circular coil with rectangular cross-section)
and exports as GMSH v4.1 format for PEEC analysis.

Workflow:
    1. Create toroidal coil geometry in Cubit
    2. Generate SURFACE mesh only (PEEC requirement)
    3. Export to GMSH v4.1 format (.msh)
    4. View in GMSH GUI: gmsh circular_coil.msh

Output:
    circular_coil.msh - Surface mesh in GMSH v4.1 format

Requirements:
    - Coreform Cubit 2025.8+
    - cubit_mesh_export module (via radia package)
"""

import sys
from pathlib import Path

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
print("Circular Coil Surface Mesh Generation (Cubit -> GMSH)")
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

# Create toroidal coil using torus
print("\n[2] Creating toroidal coil geometry...")
major_radius = mean_radius
minor_radius_x = wire_width / 2.0
minor_radius_z = wire_height / 2.0

# Create torus (circular cross-section first, then modify)
cubit.cmd(f"create torus major radius {major_radius} minor radius {minor_radius_x}")

# Scale in Z direction to get rectangular cross-section
scale_z = minor_radius_z / minor_radius_x
cubit.cmd(f"volume 1 scale z {scale_z}")

# Synchronize geometry
cubit.cmd("volume all scheme tetmesh")
cubit.cmd(f"volume all size {mesh_size}")

print("    Geometry created: Toroidal coil with rectangular cross-section")

# Generate SURFACE mesh only (PEEC requirement)
print("\n[3] Generating surface mesh (PEEC requirement: surface only)...")
cubit.cmd("mesh surface all")

# Get mesh statistics
nodes = cubit.get_node_count()
tris = cubit.get_tri_count()
print(f"    Nodes: {nodes}")
print(f"    Surface elements (triangles): {tris}")

# Define physical groups (blocks) for export
print("\n[4] Defining physical groups (blocks)...")
cubit.cmd("block 1 add tri all")
cubit.cmd("block 1 name 'conductor'")
cubit.cmd("block 1 element type tri3")
print("    Block 1: 'conductor' (surface triangles)")

# Export to GMSH v4.1 format
output_file = "circular_coil.msh"
print(f"\n[5] Exporting to GMSH v4.1 format...")
cubit.cmd(f'export gmsh "{output_file}" overwrite')
print(f"    [OK] Created: {output_file}")
print(f"    Format: GMSH v4.1 (use GMSH Python API to load)")

# Verify: should have NO volume elements
if cubit.get_tet_count() > 0 or cubit.get_hex_count() > 0:
    print("    ⚠️  WARNING: Volume elements found - PEEC only needs surface!")
else:
    print("    [OK] No volume elements - correct for PEEC")

print("\n" + "=" * 70)
print("Next Steps")
print("=" * 70)
print("\n1. View in GMSH GUI:")
print(f"   gmsh {output_file}")
print("\n2. Load into Python (using GMSH Python API):")
print("   import gmsh")
print("   gmsh.initialize()")
print(f"   gmsh.open('{output_file}')")
print("\n3. Run PEEC analysis:")
print("   cd ..")
print("   python demo_gmsh_to_peec.py")

print("\n" + "=" * 70)
print("[OK] Coil mesh generation completed!")
print("=" * 70)
