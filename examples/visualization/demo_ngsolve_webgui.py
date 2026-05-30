#!/usr/bin/env python
"""
NGSolve webgui Visualization Demo for Radia-NGSolve

Demonstrates:
- Radia field as NGSolve CoefficientFunction
- GridFunction projection for visualization
- Interactive webgui in browser (Jupyter不要)
- Mesh + field simultaneous display

Use case: Interactive field exploration, browser-based visualization

Requirements:
- radia (v2.5.0+ includes RadiaField)
- NGSolve with webgui
- Browser (Chrome/Firefox recommended)

Note: webguiは通常のPythonスクリプト（.py）から使用可能です。
      Jupyterノートブック不要で、ブラウザが自動的に開きます。
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad
from ngsolve import *
from ngsolve.webgui import Draw
from ngsolve import TaskManager
from netgen.occ import Box, Pnt, OCCGeometry

# Check if rad.RadiaField is available
try:
    rad.RadiaField
    HAS_RADIA_NGSOLVE = True
except AttributeError:
    print("Warning: rad.RadiaField not available")
    print("This demo requires radia v2.5.0+")
    HAS_RADIA_NGSOLVE = False


def main():
    print("="*60)
    print("NGSolve webgui Visualization Demo")
    print("="*60)

    if not HAS_RADIA_NGSOLVE:
        print("\n[SKIP] rad.RadiaField not available")
        print("Requires radia v2.5.0+ to run this demo")
        return


    # Create Radia magnet
    print("\n[1] Creating Radia magnet...")
    magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])
    print("    Magnet: 40x40x20 mm, Br=1.2T (NdFeB)")

    # Create NGSolve mesh (observation region)
    print("\n[2] Creating NGSolve mesh...")
    box = Box(Pnt(-0.1, -0.1, 0.02), Pnt(0.1, 0.1, 0.15))
    geo = OCCGeometry(box)
    mesh = Mesh(geo.GenerateMesh(maxh=0.02))
    print(f"    Mesh: {mesh.nv} vertices, {mesh.ne} elements")

    # Radia field as CoefficientFunction
    print("\n[3] Creating Radia CoefficientFunction...")
    B_cf = rad.RadiaField(magnet, 'b', units='m')
    print("    CoefficientFunction created")

    # Project to GridFunction for visualization
    print("\n[4] Projecting to GridFunction...")
    fes = HDiv(mesh, order=2)
    B_gf = GridFunction(fes)
    with TaskManager():
        B_gf.Set(B_cf)
        print("    GridFunction projection complete")

        # Visualization 1: Field magnitude
        print("\n[5] Visualization 1: B field magnitude")
        print("    Opening browser...")
        B_mag = sqrt(InnerProduct(B_gf, B_gf))
        Draw(B_mag, mesh, 'B_magnitude')

        # Visualization 2: Vector field
        print("\n[6] Visualization 2: B field vectors")
        print("    Opening browser...")
        Draw(B_gf, mesh, 'B_field', vectors={'grid_size': 10})

        # Visualization 3: Mesh only
        print("\n[7] Visualization 3: Mesh structure")
        print("    Opening browser...")
        Draw(mesh)

        print("\n[OK] Demo complete")
        print("\nNote: webgui opens in browser")
        print("      Close browser tabs to exit")


if __name__ == '__main__':
    main()
