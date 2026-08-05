"""
Cubit: Generate Circular Wire Centerline Mesh for PEEC Validation

Creates a straight wire centerline for validating PEEC SIBC implementation
against Bessel function analytical solution.

Author: Radia Development Team
Date: 2026-02-13
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
print("Cubit: Generate Circular Wire Centerline for PEEC SIBC Validation")
print("=" * 70)

# ============================================================================
# Wire Parameters
# ============================================================================

wire_length = 1.0    # 1 m
wire_radius = 1e-3   # 1 mm radius (2 mm diameter)
n_segments = 10      # Number of segments

print(f"\nWire parameters:")
print(f"  Length: {wire_length * 1e3:.0f} mm")
print(f"  Radius: {wire_radius * 1e3:.2f} mm (diameter: {2*wire_radius*1e3:.2f} mm)")
print(f"  Segments: {n_segments}")

# ============================================================================
# Create Centerline Geometry
# ============================================================================

print("\n[Step 1] Creating centerline curve...")

# Reset
cubit.cmd("reset")

# Create straight line along X-axis
cubit.cmd(f"create curve location 0 0 0 location {wire_length} 0 0")

centerline_curve_id = cubit.get_last_id("curve")
print(f"  Centerline curve ID: {centerline_curve_id}")

# ============================================================================
# Mesh Centerline
# ============================================================================

print("\n[Step 2] Meshing centerline...")

cubit.cmd(f"curve {centerline_curve_id} interval {n_segments}")
cubit.cmd(f"curve {centerline_curve_id} scheme equal")
cubit.cmd(f"mesh curve {centerline_curve_id}")

n_nodes = cubit.get_node_count()
n_edges = cubit.get_edge_count()

print(f"  Nodes: {n_nodes}")
print(f"  Edge elements: {n_edges}")

# ============================================================================
# Define Block (Physical Group)
# ============================================================================

print("\n[Step 3] Defining physical group...")

cubit.cmd(f"block 1 add curve {centerline_curve_id}")
cubit.cmd("block 1 name 'wire_centerline'")
# NOTE: keep the default BAR element type.  `element type BEAM2` makes the
# gmsh exporter skip the 1D elements (BEAM is not in its BAR mapping).

print(f"  Block 1: 'wire_centerline' (1D edge elements)")

# ============================================================================
# Export to GMSH v4.1
# ============================================================================

output_file = "circular_wire_centerline.msh"
print(f"\n[Step 4] Exporting to GMSH v4.1 format...")

cubit.cmd(f'export gmsh "{output_file}" overwrite')

print(f"  OK Created: {output_file}")

# ============================================================================
# Save Wire Parameters
# ============================================================================

param_file = "circular_wire_params.txt"
print(f"\n[Step 5] Saving wire parameters...")

with open(param_file, "w") as f:
    f.write(f"# Circular wire parameters for PEEC validation\n")
    f.write(f"radius {wire_radius}\n")
    f.write(f"length {wire_length}\n")
    f.write(f"sigma 5.8e7  # Copper conductivity [S/m]\n")
    f.write(f"n_segments {n_segments}\n")

print(f"  OK Created: {param_file}")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

print(f"\nGenerated mesh:")
print(f"  Wire length: {wire_length * 1e3:.0f} mm")
print(f"  Wire radius: {wire_radius * 1e3:.2f} mm")
print(f"  Segments: {n_segments}")
print(f"  Nodes: {n_nodes}")
print(f"  Elements: {n_edges}")

print(f"\nOutput files:")
print(f"  - {output_file} (GMSH mesh)")
print(f"  - {param_file} (Parameters)")

print("\nNext steps:")
print("  1. Run validate_circular_wire_sibc.py to compare with Bessel solution")
print("  2. Check error vs analytical solution")

print("\n" + "=" * 70)
print("OK Cubit mesh generation complete!")
print("=" * 70)
