#!/usr/bin/env python
"""
NGSolve webgui with Accurate Geometry

Advantage over ParaView:
- OCC shapes displayed exactly (no mesh approximation)
- Browser-based (no ParaView installation needed)
- Direct from Python script (no manual GUI operation)

Workflow:
1. Create Radia magnet
2. Convert to OCC shape (exact geometry)
3. Display in webgui (browser opens automatically)

Requirements:
- radia (v2.5.0+ includes RadiaField)
- NGSolve with webgui
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad

try:
    rad.RadiaField
    HAS_RADIA_NGSOLVE = True
except AttributeError:
    HAS_RADIA_NGSOLVE = False

from ngsolve import *
from ngsolve.webgui import Draw
from netgen.occ import Box, Cylinder, Sphere, Pnt, OCCGeometry


def radia_recmag_to_occ(center, dimensions):
    """
    Convert Radia ObjRecMag to OCC Box (exact geometry).

    Args:
        center: [x, y, z] in meters
        dimensions: [dx, dy, dz] in meters

    Returns:
        OCC Box shape
    """
    x0, y0, z0 = center
    dx, dy, dz = dimensions

    # OCC Box from corner points
    p1 = Pnt(x0 - dx/2, y0 - dy/2, z0 - dz/2)
    p2 = Pnt(x0 + dx/2, y0 + dy/2, z0 + dz/2)

    return Box(p1, p2)


def main():
    print("="*60)
    print("NGSolve webgui with Accurate Geometry")
    print("="*60)


    # Create Radia magnet
    print("\n[1] Creating Radia magnet...")
    center = [0, 0, 0]
    dimensions = [0.04, 0.04, 0.02]  # 40x40x20 mm
    magnetization = [0, 0, 954930]   # 1.2T NdFeB

    magnet = rad.ObjRecMag(center, dimensions, magnetization)
    print(f"    Magnet: {dimensions[0]*1000}x{dimensions[1]*1000}x{dimensions[2]*1000} mm, Br=1.2T")

    # Convert to OCC shape (exact geometry)
    print("\n[2] Converting to OCC shape...")
    occ_magnet = radia_recmag_to_occ(center, dimensions)
    print("    OCC Box created (exact geometry, no approximation)")

    # Display geometry in webgui
    print("\n[3] Displaying geometry in webgui...")
    print("    Opening browser...")
    Draw(occ_magnet, name='Magnet_Geometry')

    if not HAS_RADIA_NGSOLVE:
        print("\n[INFO] rad.RadiaField not available")
        print("       Geometry displayed, but field visualization requires radia v2.5.0+")
        return

    # Create observation mesh
    print("\n[4] Creating observation mesh...")
    obs_box = Box(Pnt(-0.1, -0.1, 0.02), Pnt(0.1, 0.1, 0.15))
    geo = OCCGeometry(obs_box)
    mesh = Mesh(geo.GenerateMesh(maxh=0.02))
    print(f"    Mesh: {mesh.nv} vertices, {mesh.ne} elements")

    # Radia field as CoefficientFunction
    print("\n[5] Creating Radia field...")
    B_cf = rad.RadiaField(magnet, 'b', units='m')

    # Project to GridFunction
    print("\n[6] Projecting to GridFunction...")
    fes = HDiv(mesh, order=2)
    B_gf = GridFunction(fes)
    B_gf.Set(B_cf)

    # Display field
    print("\n[7] Displaying field in webgui...")
    print("    Opening browser...")
    B_mag = sqrt(InnerProduct(B_gf, B_gf))
    Draw(B_mag, mesh, 'B_magnitude')
    Draw(B_gf, mesh, 'B_field', vectors={'grid_size': 10})

    print("\n[OK] Visualization complete")
    print("\nAdvantages:")
    print("  - Magnet geometry: exact OCC Box (no approximation)")
    print("  - Field: accurate GridFunction projection")
    print("  - Browser-based: no ParaView needed")
    print("  - Script-based: no manual GUI operation")


if __name__ == '__main__':
    main()
