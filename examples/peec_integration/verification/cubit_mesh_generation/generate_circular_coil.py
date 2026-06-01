"""
Cubit: Generate Circular Coil Centerline Mesh for PEEC Validation

Creates a single-turn circular coil centerline for validating PEEC SIBC
against analytical solutions (Neumann's formula, Grover's formula).

Author: Radia Development Team
Date: 2026-02-13
"""

import sys
import os
# Auto-detect Cubit installation
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), 'src', 'radia'))
from install_panels import find_cubit_bin as _fcb
_cubit_path = _fcb()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)
import cubit


cubit.init(['cubit', '-nojournal', '-batch'])

print("=" * 70)
print("Cubit: Generate Circular Coil Centerline for PEEC SIBC Validation")
print("=" * 70)

# ============================================================================
# Coil Parameters
# ============================================================================

coil_radius = 50e-3      # 50 mm radius (mean radius of coil)
wire_radius = 1e-3       # 1 mm wire radius (cross-section)
n_segments = 36          # Number of segments (10 degrees each)

print(f"\nCoil parameters:")
print(f"  Coil radius: {coil_radius * 1e3:.1f} mm")
print(f"  Wire radius: {wire_radius * 1e3:.2f} mm (diameter: {2*wire_radius*1e3:.2f} mm)")
print(f"  Segments: {n_segments}")

# ============================================================================
# Create Circular Centerline Geometry
# ============================================================================

print("\n[Step 1] Creating circular coil centerline...")

# Reset
cubit.cmd("reset")

# Create circle in XY plane (Z-axis normal)
# Use arc creation: create curve arc center location ... radius ... normal ...
cubit.cmd(f"create curve arc center location 0 0 0 radius {coil_radius} normal 0 0 1 full")

centerline_curve_id = cubit.get_last_id("curve")
print(f"  Centerline curve ID: {centerline_curve_id}")

# Verify circumference
circumference = 2.0 * 3.14159265359 * coil_radius
print(f"  Theoretical circumference: {circumference * 1e3:.2f} mm")

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
cubit.cmd("block 1 name 'coil_centerline'")
cubit.cmd("block 1 element type BEAM2")

print(f"  Block 1: 'coil_centerline' (1D edge elements)")

# ============================================================================
# Export to GMSH v4.1
# ============================================================================

output_file = "circular_coil_centerline.msh"
print(f"\n[Step 4] Exporting to GMSH v4.1 format...")

cubit.cmd(f'export gmsh "{output_file}" overwrite')

print(f"  OK Created: {output_file}")

# ============================================================================
# Save Coil Parameters
# ============================================================================

param_file = "circular_coil_params.txt"
print(f"\n[Step 5] Saving coil parameters...")

with open(param_file, "w") as f:
    f.write(f"# Circular coil parameters for PEEC validation\n")
    f.write(f"coil_radius {coil_radius}\n")
    f.write(f"wire_radius {wire_radius}\n")
    f.write(f"sigma 5.8e7  # Copper conductivity [S/m]\n")
    f.write(f"n_segments {n_segments}\n")
    f.write(f"n_turns 1  # Single turn\n")

print(f"  OK Created: {param_file}")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

print(f"\nGenerated mesh:")
print(f"  Coil radius: {coil_radius * 1e3:.1f} mm")
print(f"  Wire radius: {wire_radius * 1e3:.2f} mm")
print(f"  Segments: {n_segments}")
print(f"  Nodes: {n_nodes}")
print(f"  Elements: {n_edges}")

print(f"\nOutput files:")
print(f"  - {output_file} (GMSH mesh)")
print(f"  - {param_file} (Parameters)")

print("\nNext steps:")
print("  1. Run validate_circular_coil_sibc.py to compare with analytical solution")
print("  2. Check inductance vs Neumann/Grover formula")
print("  3. Check resistance vs Dowell/Bessel solution")

print("\n" + "=" * 70)
print("OK Cubit mesh generation complete!")
print("=" * 70)
