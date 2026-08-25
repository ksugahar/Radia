"""Test fixture: 2D axisymmetric workpiece mesh in (r, z).

Half-section of the cylinder R=25mm H=25mm IH workpiece.

Boundaries::

    outer : r=R curve (heating face)
    top   : z=+H/2 curve
    bot   : z=-H/2 curve
    axis  : r=0 curve (symmetry; no BC needed)

Run::

    python validation_test/panels/fixtures/generate_heat_cylinder_axisym.py

to (re)write::

    validation_test/panels/fixtures/heat_workpiece_cylinder_R25_H25_axisym.vol
"""

from __future__ import annotations

import os
import sys


def main():
    """Build the cylinder fixture as a strictly-rectangular structured
    grid in (r, z).  The Henrotte axisymmetric basis (radia 4.32.0+)
    requires axis-aligned quads exactly; Netgen's quad_dominated mesher
    produces slightly skewed quads (corner at z=-0.0125 vs opposite at
    -0.00937624) which break the closed-form Q1/Q2 element matrices.
    Per CLAUDE.md "Axisymmetric FE: Henrotte Basis Only" policy.
    """
    import netgen.meshing as ngm_pkg
    from netgen.meshing import (
        Mesh as NgMesh, MeshPoint, EdgeDescriptor, FaceDescriptor, Element1D,
        Element2D, Pnt,
    )

    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here,
                       "heat_workpiece_cylinder_R25_H25_axisym.vol")

    R, H = 0.025, 0.025
    NR, NZ = 9, 9   # quad count per side; 9 -> 8 cells -> 64 quads total

    ngm = NgMesh(dim=2)
    ngm.SetGeometry(None)

    # Vertex grid (NR+1) x (NZ+1)
    pids = {}
    for j in range(NZ + 1):
        z = -H / 2 + H * j / NZ
        for i in range(NR + 1):
            r = R * i / NR
            pids[(i, j)] = ngm.Add(MeshPoint(Pnt(r, z, 0)))

    # Face descriptor for material 'workpiece' (index 1)
    fd = ngm.Add(FaceDescriptor(bc=1, domin=1, surfnr=1))
    ngm.SetMaterial(1, "workpiece")

    # Quad elements
    for j in range(NZ):
        for i in range(NR):
            v0 = pids[(i,     j)]
            v1 = pids[(i + 1, j)]
            v2 = pids[(i + 1, j + 1)]
            v3 = pids[(i,     j + 1)]
            ngm.Add(Element2D(fd, [v0, v1, v2, v3]))

    # Boundary segments: one FaceDescriptor per named edge.
    # NGSolve looks at the bcname/index on Element1D to populate
    # mesh.GetBoundaries().
    edge_specs = [
        ("bot",   [(i,   0,    i + 1, 0    ) for i in range(NR)]),
        ("outer", [(NR,  j,    NR,    j + 1) for j in range(NZ)]),
        ("top",   [(i,   NZ,   i + 1, NZ   ) for i in range(NR)]),
        ("axis",  [(0,   j,    0,     j + 1) for j in range(NZ)]),
    ]
    for bc_idx, (name, edges) in enumerate(edge_specs, start=1):
        ngm.SetBCName(bc_idx - 1, name)
        edge_descriptor = EdgeDescriptor()
        edge_descriptor.edgenr = bc_idx
        edge_descriptor.surfnr = (bc_idx, -1)
        edge_descriptor.domin = 1
        edge_descriptor.domout = 0
        edge_descriptor.name = name
        ngm.Add(edge_descriptor)
        for (i0, j0, i1, j1) in edges:
            seg = Element1D([pids[(i0, j0)], pids[(i1, j1)]], index=bc_idx)
            ngm.Add(seg)

    ngm.Save(out)
    print(f"wrote {out}  (structured rect grid: {NR}x{NZ} quads = {NR*NZ})")


if __name__ == "__main__":
    sys.exit(main() or 0)
