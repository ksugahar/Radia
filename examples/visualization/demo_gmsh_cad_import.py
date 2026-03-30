#!/usr/bin/env python
"""
GMSH CAD Import Workflow

Demonstrates:
- GMSH can directly import CAD files (STEP, IGES, etc.)
- No need for Cubit license
- Fully open source workflow
- Automatic surface element generation

Workflow:
    CAD (STEP) -> GMSH -> .msh -> NGSolve -> Radia

Advantages over Cubit:
- [OK] Free (no license required)
- [OK] Open source
- [OK] Cross-platform
- [OK] Python scriptable

Requirements:
- gmsh (pip install gmsh)
- ngsolve
- CAD file (.step, .iges, .brep, .stl)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import gmsh
from ngsolve import Mesh


def import_cad_and_mesh(step_file, mesh_size=0.01):
    """
    Import CAD file with GMSH and generate mesh.

    Args:
        step_file: Path to STEP file
        mesh_size: Maximum element size (meters)

    Returns:
        Path to generated .msh file
    """
    print("="*60)
    print("GMSH CAD Import and Meshing")
    print("="*60)

    if not os.path.exists(step_file):
        print(f"\n[INFO] Example STEP file not found: {step_file}")
        print("       This script demonstrates the workflow.")
        print("       To use, provide your own STEP file.")
        return None

    # Initialize GMSH
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    # Import CAD file
    print(f"\n[1] Importing CAD file: {step_file}")
    gmsh.model.add("imported_model")

    try:
        # GMSH can import: STEP, IGES, BREP, STL
        gmsh.merge(step_file)
        print("    [OK] CAD file imported successfully")
    except Exception as e:
        print(f"    ❌ Failed to import: {e}")
        gmsh.finalize()
        return None

    # Synchronize
    gmsh.model.geo.synchronize()

    # Set mesh size
    print(f"\n[2] Setting mesh size: {mesh_size} m")
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)

    # Define physical groups (important for NGSolve)
    print("\n[3] Defining physical groups...")

    # Get all volumes
    volumes = gmsh.model.getEntities(3)
    if volumes:
        volume_tags = [v[1] for v in volumes]
        gmsh.model.addPhysicalGroup(3, volume_tags, 1)
        gmsh.model.setPhysicalName(3, 1, "domain")
        print(f"    Volumes: {len(volume_tags)}")

    # Get all surfaces
    surfaces = gmsh.model.getEntities(2)
    if surfaces:
        surface_tags = [s[1] for s in surfaces]
        gmsh.model.addPhysicalGroup(2, surface_tags, 1)
        gmsh.model.setPhysicalName(2, 1, "boundary")
        print(f"    Surfaces: {len(surface_tags)}")

    # Generate mesh
    print("\n[4] Generating mesh...")
    gmsh.model.mesh.generate(3)

    # Get statistics
    nodes = gmsh.model.mesh.getNodes()
    elements = gmsh.model.mesh.getElements(3)  # Volume elements

    print(f"    Nodes: {len(nodes[0])}")
    if elements[1]:
        print(f"    Volume elements: {len(elements[1][0])}")

    # Check surface elements
    surf_elements = gmsh.model.mesh.getElements(2)  # Surface elements
    if surf_elements[1]:
        total_surf = sum(len(e) for e in surf_elements[1])
        print(f"    Surface elements: {total_surf}")
        print("    [OK] Surface elements present - Netgen GUI compatible")

    # Save mesh
    msh_file = step_file.replace('.step', '.msh').replace('.stp', '.msh')
    gmsh.write(msh_file)
    print(f"\n[5] Saved mesh: {msh_file}")

    # Optionally open in GMSH GUI
    # gmsh.fltk.run()

    gmsh.finalize()

    return msh_file


def workflow_example():
    """Demonstrate complete workflow."""
    print("="*60)
    print("Complete Open Source Workflow")
    print("="*60)
    print("\nWorkflow:")
    print("  CAD (STEP) -> GMSH -> .msh -> NGSolve -> Radia")
    print("\nAdvantages:")
    print("  [OK] 100% open source (no Cubit license)")
    print("  [OK] Cross-platform (Windows/Linux/Mac)")
    print("  [OK] Python scriptable")
    print("  [OK] Surface elements automatically generated")
    print("  [OK] Can import: STEP, IGES, BREP, STL")

    # Example: If you have a STEP file
    step_file = 'motor_rotor.step'  # Replace with your STEP file

    if not os.path.exists(step_file):
        print(f"\n[INFO] Example STEP file '{step_file}' not found")
        print("       To use this workflow:")
        print("       1. Export CAD model as STEP file")
        print("       2. Run: python demo_gmsh_cad_import.py")
        print("       3. GMSH will import and mesh automatically")
        return

    # Import and mesh
    msh_file = import_cad_and_mesh(step_file, mesh_size=0.005)

    if msh_file:
        print("\n" + "="*60)
        print("Next steps:")
        print("="*60)
        print(f"\n1. View in GMSH GUI:")
        print(f"   gmsh {msh_file}")

        print(f"\n2. Import to NGSolve:")
        print(f"   from ngsolve import Mesh")
        print(f"   mesh = Mesh('{msh_file}')")

        print(f"\n3. Convert to .vol (optional):")
        print(f"   mesh.ngmesh.Save('model.vol')")

        print(f"\n4. View in Netgen GUI:")
        print(f"   python utils/netgen_vol_viewer.py model.vol")

        print(f"\n5. Use with Radia:")
        print(f"   from netgen_mesh_import import netgen_mesh_to_radia")
        print(f"   radia_obj = netgen_mesh_to_radia(mesh, ...)")


def comparison_table():
    """Print comparison of mesh generation tools."""
    print("\n" + "="*60)
    print("Mesh Tool Comparison")
    print("="*60)
    print("\n| Tool | License | GUI | Python API | CAD Import | Platforms |")
    print("|------|---------|-----|------------|------------|-----------|")
    print("| **GMSH** | Free/OSS | [OK] | [OK] | [OK] STEP/IGES | Win/Lin/Mac |")
    print("| **Netgen** | Free/OSS | [OK] | [OK] | [OK] STEP/OCC | Win/Lin/Mac |")
    print("| **Cubit** | Commercial* | [OK] | [OK] | [OK] STEP/IGES | Win/Lin/Mac |")
    print("\n*Cubit: Free for academic use with restrictions")

    print("\n" + "="*60)
    print("Recommendation for Open Source Workflow")
    print("="*60)
    print("\nSimple geometry:")
    print("  -> Netgen (Python API: netgen.occ)")
    print("\nComplex geometry from CAD:")
    print("  -> GMSH (import STEP, auto-mesh)")
    print("\nHexahedral meshes:")
    print("  -> Cubit (required, Netgen cannot do 3D hex)")
    print("\nFully open source preference:")
    print("  -> GMSH + Netgen")


def main():
    workflow_example()
    comparison_table()

    print("\n" + "="*60)
    print("Summary")
    print("="*60)
    print("\nGMSH経由は完全に有効な選択肢です:")
    print("  [OK] オープンソース（無料）")
    print("  [OK] クロスプラットフォーム")
    print("  [OK] CAD直接読込（STEP, IGES）")
    print("  [OK] 表面要素自動生成")
    print("  [OK] NGSolve統合")
    print("\n特に、Cubitライセンスが不要な点が大きなメリットです。")


if __name__ == '__main__':
    main()
