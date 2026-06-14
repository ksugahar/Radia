"""Probe whether mesh.Curve(p) actually moves the trafo's mid-edge position.

Builds a small SplineGeometry circle mesh, picks one boundary element,
and reports trafo evaluation at ref(0.5, 0.5) for Curve order 0..5.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

from netgen.geom2d import SplineGeometry
from ngsolve import Mesh, IntegrationRule
from ngsolve.fem import ET

g = SplineGeometry()
R = 1.0
p_bot = g.AppendPoint(0.0, -R)
p_apex = g.AppendPoint(R, 0.0)
p_top = g.AppendPoint(0.0, R)
g.Append(["spline3", p_bot, p_apex, p_top], bc="arc",
          leftdomain=1, rightdomain=0)
g.Append(["line", p_top, p_bot], bc="axis",
          leftdomain=1, rightdomain=0)
g.SetMaterial(1, "conductor")

m = Mesh(g.GenerateMesh(maxh=0.5))

# Find a triangle whose edge touches the "arc" boundary
target_el = None
for el in m.Elements():
    if any(v.point[0]**2 + v.point[1]**2 > 0.7**2 for v in el.vertices):
        target_el = el
        break

print(f"Triangle vertices: {[(v.point[0], v.point[1]) for v in target_el.vertices]}")

ir = IntegrationRule(ET.TRIG, 6)

for order in [0, 1, 2, 3, 4, 5]:
    m.Curve(order)
    trafo = m.GetTrafo(target_el)
    # Sample at (0.5, 0.5), (0.5, 0), (0, 0.5)
    samples = {}
    for ip in ir:
        xi, eta = float(ip.point[0]), float(ip.point[1])
        if abs(xi - 0.5) < 0.05 and abs(eta - 0.5) < 0.05:
            samples["m_V0_V1"] = trafo(ip).point
        elif abs(xi - 0.5) < 0.05 and abs(eta) < 0.05:
            samples["m_V2_V0"] = trafo(ip).point
        elif abs(xi) < 0.05 and abs(eta - 0.5) < 0.05:
            samples["m_V1_V2"] = trafo(ip).point
    print(f"Curve({order}): samples = {samples}")
