"""Verify NGSolve's edge[k] convention for ET_TRIG.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

import netgen.meshing as ng_meshing
from netgen.meshing import (
    Mesh as NgMesh, Element1D, Element2D, FaceDescriptor, MeshPoint, Pnt,
)
from ngsolve import Mesh

m = NgMesh()
m.dim = 2
m.SetMaterial(1, "conductor")
m.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
m.SetBCName(0, "outer")
p0 = m.Add(MeshPoint(Pnt(1, 0, 0)))
p1 = m.Add(MeshPoint(Pnt(5, 0, 0)))
p2 = m.Add(MeshPoint(Pnt(5, 2, 0)))
m.Add(Element2D(1, [p0, p1, p2]))
m.Add(Element1D([p0, p1], index=1))
m.Add(Element1D([p1, p2], index=1))
m.Add(Element1D([p2, p0], index=1))
mesh = Mesh(m)

for el in mesh.Elements():
    verts = list(el.vertices)
    edges = list(el.edges)
    print(f"Triangle vertex indices: {[v.nr for v in verts]}")
    print(f"Triangle edges (in element-local order):")
    for k, e in enumerate(edges):
        # NGSolve global edge: look up by index
        ge = mesh[e]
        ev = list(ge.vertices)
        print(f"  edges[{k}] = edge index {e.nr}, connects vertices {[v.nr for v in ev]}")
