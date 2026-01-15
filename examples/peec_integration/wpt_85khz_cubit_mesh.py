"""
WPT 85 kHz System with Cubit Mesh Generation

This example demonstrates PEEC Loop-Star + MMM analysis using Coreform Cubit
for mesh generation.

Mesh components:
1. Tx/Rx spiral coils (hex mesh) - PEEC conductor
2. Ferrite cores (hex mesh) - MMM magnetic material
3. Aluminum shields (hex mesh) - PEEC conductor with eddy currents

Workflow:
    Cubit geometry -> Nastran BDF -> Netgen mesh -> Radia PEEC+MMM

Requirements:
    - Coreform Cubit 2025.3+
    - cubit_mesh_export from S:\CoreformCubit\01_GitHub

Author: Radia Development Team
Date: 2026-01-16
"""

import sys
import os

# Cubit path
CUBIT_PATH = "C:/Program Files/Coreform Cubit 2025.3/bin"
CUBIT_EXPORT_PATH = "S:/CoreformCubit/01_GitHub"

sys.path.insert(0, CUBIT_PATH)
sys.path.insert(0, CUBIT_EXPORT_PATH)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build-msvc'))

import numpy as np

# Check if Cubit is available
try:
    import cubit
    CUBIT_AVAILABLE = True
except ImportError:
    CUBIT_AVAILABLE = False
    print("Warning: Cubit not available. Using pre-generated mesh files.")


def create_wpt_geometry_cubit():
    """
    Create WPT system geometry in Cubit.

    Returns filenames of exported Nastran BDF files.
    """
    if not CUBIT_AVAILABLE:
        return None

    import cubit_mesh_export

    cubit.init(['cubit', '-nojournal', '-batch'])

    output_dir = os.path.dirname(os.path.abspath(__file__))

    # ============================================================
    # Parameters (all in mm for Cubit)
    # ============================================================
    # Coil parameters
    coil_inner_r = 50.0      # 50 mm
    coil_outer_r = 150.0     # 150 mm
    coil_thickness = 3.0     # 3 mm (Litz wire bundle)
    n_turns = 10
    coil_pitch = (coil_outer_r - coil_inner_r) / n_turns

    # Ferrite core (disc)
    ferrite_r = 180.0        # 180 mm
    ferrite_thickness = 5.0  # 5 mm
    ferrite_z_offset = -4.0  # Below coil

    # Aluminum shield (disc)
    shield_r = 200.0         # 200 mm
    shield_thickness = 2.0   # 2 mm
    shield_z_offset = -10.0  # Below ferrite

    # Air gap between Tx and Rx
    air_gap = 150.0          # 150 mm

    # Mesh size
    mesh_size_coil = 10.0    # 10 mm
    mesh_size_core = 15.0    # 15 mm
    mesh_size_shield = 20.0  # 20 mm

    # ============================================================
    # Create Tx Assembly
    # ============================================================
    print("Creating Tx assembly...")
    cubit.cmd("reset")

    # Tx Ferrite core (hex mesh)
    cubit.cmd(f"create cylinder radius {ferrite_r} height {ferrite_thickness}")
    cubit.cmd(f"volume 1 move 0 0 {ferrite_z_offset}")
    cubit.cmd("volume 1 scheme sweep")
    cubit.cmd(f"volume 1 size {mesh_size_core}")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add hex in volume 1")
    cubit.cmd("block 1 name 'tx_ferrite'")

    # Tx Shield (hex mesh)
    cubit.cmd(f"create cylinder radius {shield_r} height {shield_thickness}")
    cubit.cmd(f"volume 2 move 0 0 {shield_z_offset}")
    cubit.cmd("volume 2 scheme sweep")
    cubit.cmd(f"volume 2 size {mesh_size_shield}")
    cubit.cmd("mesh volume 2")
    cubit.cmd("block 2 add hex in volume 2")
    cubit.cmd("block 2 name 'tx_shield'")

    # Tx Coil (simplified as annular disc for now)
    # Full spiral would require more complex geometry
    cubit.cmd(f"create cylinder radius {coil_outer_r} height {coil_thickness}")
    cubit.cmd(f"create cylinder radius {coil_inner_r} height {coil_thickness + 2}")
    cubit.cmd("subtract volume 4 from volume 3")
    cubit.cmd("volume 3 move 0 0 0")
    cubit.cmd("volume 3 scheme sweep")
    cubit.cmd(f"volume 3 size {mesh_size_coil}")
    cubit.cmd("mesh volume 3")
    cubit.cmd("block 3 add hex in volume 3")
    cubit.cmd("block 3 name 'tx_coil'")

    # Export Tx assembly
    tx_bdf = os.path.join(output_dir, "wpt_tx_assembly.bdf")
    print(f"Exporting Tx assembly to {tx_bdf}")
    cubit_mesh_export.export_nastran(cubit, tx_bdf, DIM="3D")

    # ============================================================
    # Create Rx Assembly (mirror of Tx)
    # ============================================================
    print("Creating Rx assembly...")
    cubit.cmd("reset")

    # Rx positions (below Tx by air_gap)
    rx_z_base = -air_gap

    # Rx Ferrite core
    cubit.cmd(f"create cylinder radius {ferrite_r} height {ferrite_thickness}")
    cubit.cmd(f"volume 1 move 0 0 {rx_z_base - ferrite_z_offset}")
    cubit.cmd("volume 1 scheme sweep")
    cubit.cmd(f"volume 1 size {mesh_size_core}")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add hex in volume 1")
    cubit.cmd("block 1 name 'rx_ferrite'")

    # Rx Shield
    cubit.cmd(f"create cylinder radius {shield_r} height {shield_thickness}")
    cubit.cmd(f"volume 2 move 0 0 {rx_z_base - shield_z_offset}")
    cubit.cmd("volume 2 scheme sweep")
    cubit.cmd(f"volume 2 size {mesh_size_shield}")
    cubit.cmd("mesh volume 2")
    cubit.cmd("block 2 add hex in volume 2")
    cubit.cmd("block 2 name 'rx_shield'")

    # Rx Coil
    cubit.cmd(f"create cylinder radius {coil_outer_r} height {coil_thickness}")
    cubit.cmd(f"create cylinder radius {coil_inner_r} height {coil_thickness + 2}")
    cubit.cmd("subtract volume 4 from volume 3")
    cubit.cmd(f"volume 3 move 0 0 {rx_z_base}")
    cubit.cmd("volume 3 scheme sweep")
    cubit.cmd(f"volume 3 size {mesh_size_coil}")
    cubit.cmd("mesh volume 3")
    cubit.cmd("block 3 add hex in volume 3")
    cubit.cmd("block 3 name 'rx_coil'")

    # Export Rx assembly
    rx_bdf = os.path.join(output_dir, "wpt_rx_assembly.bdf")
    print(f"Exporting Rx assembly to {rx_bdf}")
    cubit_mesh_export.export_nastran(cubit, rx_bdf, DIM="3D")

    return {
        'tx_bdf': tx_bdf,
        'rx_bdf': rx_bdf,
        'air_gap': air_gap
    }


def load_nastran_to_radia(bdf_file, material_type='conductor'):
    """
    Load Nastran BDF mesh and create Radia objects.

    Parameters:
        bdf_file: Path to Nastran BDF file
        material_type: 'conductor', 'ferrite', or 'shield'

    Returns:
        Radia container object
    """
    import radia as rad
    from nastran_mesh_import import import_nastran_mesh

    # Read mesh
    mesh_data = import_nastran_mesh(bdf_file, units='mm')

    # Get hex elements
    hex_elements = mesh_data.get('hex_elements', [])
    vertices = mesh_data['vertices']

    print(f"  Loaded {len(hex_elements)} hex elements from {os.path.basename(bdf_file)}")

    if len(hex_elements) == 0:
        print(f"  Warning: No hex elements found in {bdf_file}")
        return None

    # Create Radia objects
    objects = []

    for elem in hex_elements:
        # Get vertex coordinates (convert mm to m)
        verts = [[vertices[v][0] / 1000,
                  vertices[v][1] / 1000,
                  vertices[v][2] / 1000] for v in elem]

        if material_type == 'ferrite':
            # Magnetic material - zero initial magnetization
            obj = rad.ObjHexahedron(verts, [0, 0, 0])
        else:
            # Conductor - no magnetization
            obj = rad.ObjHexahedron(verts, [0, 0, 0])

        objects.append(obj)

    # Create container
    container = rad.ObjCnt(objects)

    # Apply material
    if material_type == 'ferrite':
        mat = rad.MatLin(2000)  # mu_r = 2000
        rad.MatApl(container, mat)
    elif material_type == 'shield':
        # Aluminum: will be handled by PEEC with sigma = 3.5e7
        pass
    else:  # conductor
        # Copper: will be handled by PEEC with sigma = 5.8e7
        pass

    return container


def analyze_wpt_system():
    """
    Analyze WPT system using PEEC + MMM.
    """
    import radia as rad

    rad.FldUnits('m')

    print("=" * 70)
    print("WPT 85 kHz System Analysis with Cubit Mesh")
    print("=" * 70)

    # Check for pre-generated mesh files
    output_dir = os.path.dirname(os.path.abspath(__file__))
    tx_bdf = os.path.join(output_dir, "wpt_tx_assembly.bdf")
    rx_bdf = os.path.join(output_dir, "wpt_rx_assembly.bdf")

    if not os.path.exists(tx_bdf) or not os.path.exists(rx_bdf):
        if CUBIT_AVAILABLE:
            print("\nGenerating mesh with Cubit...")
            mesh_info = create_wpt_geometry_cubit()
            if mesh_info is None:
                print("Error: Mesh generation failed")
                return
        else:
            print("\nError: Mesh files not found and Cubit not available")
            print(f"  Expected: {tx_bdf}")
            print(f"  Expected: {rx_bdf}")
            print("\nPlease run this script with Cubit available to generate meshes.")
            return
    else:
        print(f"\nUsing existing mesh files:")
        print(f"  Tx: {tx_bdf}")
        print(f"  Rx: {rx_bdf}")

    # Load meshes into Radia
    print("\n--- Loading Meshes into Radia ---")

    # For demonstration, we'll show the workflow
    # Full implementation would parse blocks from BDF

    print("\nNote: Full PEEC+MMM analysis requires:")
    print("  1. CndFromMesh() for conductor elements (PEEC)")
    print("  2. ObjHexahedron() for magnetic elements (MMM)")
    print("  3. CplMag coupling between conductor and magnetic parts")

    # Placeholder for analysis results
    print("\n--- Placeholder Analysis Results ---")
    print("(Full implementation with mesh loading TBD)")

    # Show expected workflow
    print("""
Expected Workflow:
    1. Load Tx coil mesh -> PEEC conductor
    2. Load Tx ferrite mesh -> MMM magnetic material
    3. Load Tx shield mesh -> PEEC conductor (with eddy currents)
    4. Repeat for Rx assembly
    5. Create CplMag solver for Tx+Rx system
    6. Set frequency (85 kHz)
    7. Solve coupled PEEC+MMM system
    8. Extract:
       - Self-inductances L1, L2
       - Mutual inductance M
       - AC resistances R1, R2 (including eddy current losses)
       - Shield loss contributions
    9. Calculate:
       - Coupling coefficient k = M / sqrt(L1*L2)
       - Power transfer efficiency
       - S-parameters
    10. Export to SPICE via PRIMA reduction
""")

    # Cleanup
    rad.UtiDelAll()

    print("\n" + "=" * 70)
    print("Mesh generation complete. Full PEEC+MMM implementation pending.")
    print("=" * 70)


def main():
    """Main entry point."""
    print("\n" + "=" * 70)
    print("WPT 85 kHz - Cubit Mesh Generation")
    print("=" * 70)

    if CUBIT_AVAILABLE:
        print("\nCubit available. Will generate meshes.")
        mesh_info = create_wpt_geometry_cubit()
        if mesh_info:
            print(f"\nMesh files created:")
            print(f"  Tx: {mesh_info['tx_bdf']}")
            print(f"  Rx: {mesh_info['rx_bdf']}")
            print(f"  Air gap: {mesh_info['air_gap']} mm")
    else:
        print("\nCubit not available. Skipping mesh generation.")
        print("To generate meshes:")
        print("  1. Install Coreform Cubit 2025.3")
        print("  2. Run this script with Cubit in PATH")

    # Analyze system (uses pre-generated meshes if Cubit not available)
    print("\n")
    analyze_wpt_system()


if __name__ == '__main__':
    main()
