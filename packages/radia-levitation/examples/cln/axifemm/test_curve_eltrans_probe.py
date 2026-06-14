"""Probe whether mesh.Curve(p) actually moves the eltrans evaluation at
mid-edge reference points for an OCC half-circle mesh.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

import math
from netgen.occ import OCCGeometry, WorkPlane
from ngsolve import Mesh, IntegrationRule
from ngsolve.fem import ET
from ngsolve import TaskManager

R = 0.01
wp = WorkPlane().MoveTo(0, -R).Direction(1, 0).Arc(R, 180).LineTo(0, -R)
face = wp.Face()
face.name = "cond"
geo = OCCGeometry(face, dim=2)
ngmesh = geo.GenerateMesh(maxh=0.005)
mesh = Mesh(ngmesh)
print(f"mesh.ne={mesh.ne}, mesh.nv={mesh.nv}")

# Find an element whose vertices touch the curved arc (r large).
target = None
for el in mesh.Elements():
    pts = [mesh[v].point for v in el.vertices]
    if any(p[0] > 0.7 * R for p in pts):
        target = el
        break
verts_pts = [mesh[v].point for v in target.vertices]
print(f"Target triangle vertices: {[(round(p[0], 6), round(p[1], 6)) for p in verts_pts]}")

# Reference midpoints
ir = IntegrationRule(ET.TRIG, 4)

print(f"\n{'Curve':>7s}  midpoint xi-eta  -> physical (r, z)")
print("-" * 60)
for order in [0, 1, 2, 3, 4, 5]:
    with TaskManager():
        mesh.Curve(order)
        trafo = mesh.GetTrafo(target)
        samples = {}
        for ip in ir:
            xi, eta = float(ip.point[0]), float(ip.point[1])
            if abs(xi - 0.5) < 0.05 and abs(eta - 0.5) < 0.05:
                samples["(0.5, 0.5)"] = trafo(ip).point
            elif abs(xi) < 0.05 and abs(eta - 0.5) < 0.05:
                samples["(0.0, 0.5)"] = trafo(ip).point
            elif abs(xi - 0.5) < 0.05 and abs(eta) < 0.05:
                samples["(0.5, 0.0)"] = trafo(ip).point
        for k, v in sorted(samples.items()):
            print(f"  Curve({order})  {k}    -> ({v[0]:.6f}, {v[1]:.6f})")
