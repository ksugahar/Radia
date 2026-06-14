"""Determine NGSolve ET_TRIG reference triangle vertex convention.

For a known 1-triangle mesh, evaluate the element transformation at the
3 reference vertex positions (0,0), (1,0), (0,1) and see which physical
mesh vertex each one corresponds to.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

import netgen.meshing as ng_meshing
from netgen.meshing import (
    Mesh as NgMesh, Element1D, Element2D, FaceDescriptor, MeshPoint, Pnt,
)
from ngsolve import Mesh, IntegrationRule
from ngsolve.fem import ET

# Build a 1-triangle mesh with very distinguishable vertex coords.
m = NgMesh()
m.dim = 2
m.SetMaterial(1, "conductor")
m.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
m.SetBCName(0, "outer")
# Mesh vertex 0 at (10, 0); 1 at (20, 0); 2 at (15, 30).
p0 = m.Add(MeshPoint(Pnt(10.0, 0.0,  0)))
p1 = m.Add(MeshPoint(Pnt(20.0, 0.0,  0)))
p2 = m.Add(MeshPoint(Pnt(15.0, 30.0, 0)))
m.Add(Element2D(1, [p0, p1, p2]))
m.Add(Element1D([p0, p1], index=1))
m.Add(Element1D([p1, p2], index=1))
m.Add(Element1D([p2, p0], index=1))

mesh = Mesh(m)
print(f"Mesh vertex 0 (supplied via p0): {tuple(mesh.vertices[0].point)}")
print(f"Mesh vertex 1 (supplied via p1): {tuple(mesh.vertices[1].point)}")
print(f"Mesh vertex 2 (supplied via p2): {tuple(mesh.vertices[2].point)}")

# Element transformation for the only triangle.
for el in mesh.Elements():
    trafo = mesh.GetTrafo(el)  # may exist? Otherwise need different API
    print(f"\nElement vertices reported: {[v.nr for v in el.vertices]}")
    print(f"Querying transformation at reference points:")
    # Use IntegrationRule via integration_rule construction or directly
    # call eltrans on (xi, eta, w)-tuples? Simplest: use IntegrationRule.
    irule = IntegrationRule(ET.TRIG, 4)
    for ip_obj in irule:
        pt = trafo(ip_obj).point
        try:
            xi, eta = float(ip_obj.point[0]), float(ip_obj.point[1])
        except Exception:
            xi, eta = -1, -1
        print(f"  ref ({xi:.3f}, {eta:.3f})  w={float(ip_obj.weight):.4f}  -> phys ({pt[0]:.3f}, {pt[1]:.3f})")
