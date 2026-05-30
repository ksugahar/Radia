"""
WPT 85 kHz System with Cubit Mesh Generation

This example demonstrates PEEC Loop-Star + MMM analysis using Coreform Cubit
for mesh generation.

Mesh components:
1. Tx/Rx spiral coils (hex mesh) - PEEC conductor
2. Ferrite cores (hex mesh) - MMM magnetic material
3. Aluminum shields (hex mesh) - PEEC conductor with eddy currents

Workflow:
    Cubit geometry -> NGSolve mesh (direct via export_NGSolveCurvedMesh) -> Radia PEEC+MMM

Note: Netgen alone cannot create 3D hexahedral meshes.
      Cubit is required for hex mesh generation.
      Use cubit_mesh_export.extract_curved_mesh() for direct Cubit -> NGSolve conversion.

Requirements:
    - Coreform Cubit 2025.8+
    - cubit_mesh_export module (via radia package)

Author: Radia Development Team
Date: 2026-01-16
"""

import sys
import os

# Auto-detect Cubit installation
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), 'src', 'radia'))
from install_panels import find_cubit_bin
_cubit_path = find_cubit_bin()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.insert(0, _cubit_path)

# numpy not needed for mesh generation

# Check if Cubit is available
try:
    import cubit
    from cubit_mesh_export import extract_curved_mesh
    CUBIT_AVAILABLE = True
except ImportError:
    CUBIT_AVAILABLE = False
    print("Warning: Cubit not available. Using pre-generated mesh files.")


def create_wpt_geometry_cubit():
    """
    Create WPT system geometry in Cubit.

    Returns NGSolve Mesh objects directly (no intermediate file format).

    Uses cubit_mesh_export.extract_curved_mesh() for direct Cubit -> NGSolve conversion.
    """
    if not CUBIT_AVAILABLE:
        return None

    from ngsolve import Mesh
    from ngsolve import TaskManager
    cubit.init(['cubit', '-nojournal', '-batch'])

    # ============================================================
    # Parameters (all in mm for Cubit)
    # ============================================================
    # Coil parameters
    coil_inner_r = 50.0      # 50 mm
    coil_outer_r = 150.0     # 150 mm
    coil_thickness = 3.0     # 3 mm (Litz wire bundle)
    n_turns = 10  # Number of turns (for future spiral coil implementation)
    _ = n_turns  # Suppress unused warning

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

    # Export Tx assembly to NGSolve directly (no Nastran intermediate)
    print("Exporting Tx assembly to NGSolve mesh...")
    with TaskManager():
        tx_mesh = Mesh(extract_curved_mesh(cubit, order=1))

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

        # Export Rx assembly to NGSolve directly (no Nastran intermediate)
        print("Exporting Rx assembly to NGSolve mesh...")
        rx_mesh = Mesh(extract_curved_mesh(cubit, order=1))

        return {
            'tx_mesh': tx_mesh,
            'rx_mesh': rx_mesh,
            'air_gap': air_gap
        }


def load_mesh_to_radia(mesh, material_type='conductor'):
    """
    Load NGSolve mesh and create Radia objects.

    Parameters:
        mesh: NGSolve Mesh object (from cubit_mesh_export.extract_curved_mesh())
        material_type: 'conductor', 'ferrite', or 'shield'

    Returns:
        Radia container object
    """
    import radia as rad
    from netgen_mesh_import import netgen_mesh_to_radia

    # Use netgen_mesh_to_radia for direct conversion
    # This handles hex/tet elements automatically
    magnetization = [0, 0, 0]  # Zero initial magnetization

    container = netgen_mesh_to_radia(
        mesh,
        material={'magnetization': magnetization},
        units='m'  # Cubit exports in mm, but export_NGSolveCurvedMesh converts to m
    )

    if container is None:
        print(f"  Warning: No elements found in mesh")
        return None

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


    print("=" * 70)
    print("WPT 85 kHz System Analysis with Cubit Mesh")
    print("=" * 70)

    # Generate mesh with Cubit (direct Netgen export)
    if CUBIT_AVAILABLE:
        print("\nGenerating mesh with Cubit -> Netgen direct export...")
        mesh_info = create_wpt_geometry_cubit()
        if mesh_info is None:
            print("Error: Mesh generation failed")
            return
        print(f"  Air gap: {mesh_info['air_gap']} mm")
    else:
        print("\nError: Cubit not available")
        print("This example requires Coreform Cubit for hex mesh generation.")
        print("\nTo run this example:")
        print("  1. Install Coreform Cubit 2025.8+")
        print("  2. Ensure cubit_mesh_export module is available")
        return

    # Load meshes into Radia via Netgen
    print("\n--- Loading Netgen Meshes into Radia ---")

    # Convert Netgen meshes to Radia objects
    # Full implementation would use mesh_info['tx_mesh'] and mesh_info['rx_mesh']

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
    print("WPT 85 kHz - Cubit -> Netgen Direct Export")
    print("=" * 70)

    if CUBIT_AVAILABLE:
        print("\nCubit available. Will generate meshes directly to Netgen format.")
        mesh_info = create_wpt_geometry_cubit()
        if mesh_info:
            print(f"\nNGSolve meshes created in memory:")
            print(f"  Tx mesh: {type(mesh_info['tx_mesh'])}")
            print(f"  Rx mesh: {type(mesh_info['rx_mesh'])}")
            print(f"  Air gap: {mesh_info['air_gap']} mm")
    else:
        print("\nCubit not available. Cannot generate hex meshes.")
        print("To generate meshes:")
        print("  1. Install Coreform Cubit 2025.8+")
        print("  2. Ensure cubit_mesh_export module is available")
        print("\nNote: Netgen alone cannot create 3D hexahedral meshes.")

    # Analyze system
    print("\n")
    analyze_wpt_system()


if __name__ == '__main__':
    main()
