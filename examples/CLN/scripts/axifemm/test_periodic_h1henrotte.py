"""Smoke test: does Periodic(H1Henrotte(...)) work?

Builds a trivial 2-domain mesh with Identify on one boundary edge,
constructs Periodic(H1Henrotte) and verifies ndof, FreeDofs, and that
a simple BilinearForm assembles.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"S:/Radia/01_GitHub/src/radia")

print("step 1: import", flush=True)
from math import pi
from netgen.occ import (
    WorkPlane, MoveTo, Vertex, Glue, OCCGeometry, Pnt, IdentificationType,
)
from ngsolve import (
    Mesh, BilinearForm, CoefficientFunction, Periodic, TaskManager, ngsglobals,
)
from radia.radia_axifemm import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)
print("step 2: build geometry", flush=True)

R = 0.05
Z_OFFSET = 5 * R

wp_in = WorkPlane()
inner_full = wp_in.Circle(R).Face()
inner_full.name = "air_inner"
wp_out = WorkPlane().MoveTo(0, Z_OFFSET)
outer_full = wp_out.Circle(R).Face()
outer_full.name = "kelvin"

margin = 0.01
cut_in = MoveTo(-R - margin, -R - margin).Rectangle(
    R + margin, 2 * R + 2 * margin).Face()
cut_out = MoveTo(-R - margin, Z_OFFSET - R - margin).Rectangle(
    R + margin, 2 * R + 2 * margin).Face()

inner = inner_full - cut_in
outer = outer_full - cut_out

# Name edges
from math import sqrt as msqrt
kelvin_int_edges = []
for edge in inner.edges:
    cx = edge.center.x
    try:
        v0, v1 = edge.vertices
        d0 = msqrt(v0.p.x ** 2 + v0.p.y ** 2)
        d1 = msqrt(v1.p.x ** 2 + v1.p.y ** 2)
        is_arc = abs(d0 - R) < 1e-3 and abs(d1 - R) < 1e-3 and cx > 1e-3
    except Exception:
        is_arc = False
    if cx < 1e-4:
        edge.name = "axis"
    elif is_arc:
        edge.name = "kelvin_int"
        kelvin_int_edges.append(edge)
    else:
        edge.name = "default"

kelvin_ext_edges = []
for edge in outer.edges:
    cx = edge.center.x
    try:
        v0, v1 = edge.vertices
        d0 = msqrt(v0.p.x ** 2 + (v0.p.y - Z_OFFSET) ** 2)
        d1 = msqrt(v1.p.x ** 2 + (v1.p.y - Z_OFFSET) ** 2)
        is_arc = abs(d0 - R) < 1e-3 and abs(d1 - R) < 1e-3 and cx > 1e-3
    except Exception:
        is_arc = False
    if cx < 1e-4:
        edge.name = "axis_ext"
    elif is_arc:
        edge.name = "kelvin_ext"
        kelvin_ext_edges.append(edge)
    else:
        edge.name = "default"

print(f"  kelvin_int edges: {len(kelvin_int_edges)}", flush=True)
print(f"  kelvin_ext edges: {len(kelvin_ext_edges)}", flush=True)

# Pair by y-sign
matched = 0
for int_e in kelvin_int_edges:
    iy = int_e.center.y
    for ext_e in kelvin_ext_edges:
        ey = ext_e.center.y - Z_OFFSET
        if (iy > 0 and ey > 0) or (iy < 0 and ey < 0):
            int_e.Identify(ext_e, "kelvin", IdentificationType.PERIODIC)
            matched += 1
            break
print(f"  matched periodic pairs: {matched}", flush=True)

gnd = Vertex(Pnt(0, Z_OFFSET, 0))
gnd.name = "GND"

shape = Glue([inner, outer, gnd])
geo = OCCGeometry(shape, dim=2)
print("step 3: mesh", flush=True)
ngsglobals.msg_level = 0
mesh = Mesh(geo.GenerateMesh(maxh=0.02))
print(f"  mesh.ne={mesh.ne}, materials={mesh.GetMaterials()}", flush=True)
print(f"  boundaries={mesh.GetBoundaries()}", flush=True)

# Inspect periodic node pairs at mesh level
from ngsolve.comp import NodeId
ngmesh_native = mesh.ngmesh
nperiodic = mesh.nperiodic if hasattr(mesh, "nperiodic") else None
print(f"  mesh.nperiodic (if any): {nperiodic}", flush=True)
from ngsolve.fem import NODE_TYPE
try:
    pairs_v = mesh.GetPeriodicNodePairs(NODE_TYPE.VERTEX)
    print(f"  vertex periodic pairs: {len(pairs_v)}", flush=True)
    if pairs_v:
        print(f"    first 3: {pairs_v[:3]}", flush=True)
except Exception as e:
    print(f"  GetPeriodicNodePairs(VERTEX) FAIL: {e}", flush=True)
try:
    pairs_e = mesh.GetPeriodicNodePairs(NODE_TYPE.EDGE)
    print(f"  edge periodic pairs: {len(pairs_e)}", flush=True)
    if pairs_e:
        print(f"    first 3: {pairs_e[:3]}", flush=True)
except Exception as e:
    print(f"  GetPeriodicNodePairs(EDGE) FAIL: {e}", flush=True)

print("step 4: fes = H1Henrotte(order=2)", flush=True)
fes_pre = H1Henrotte(
    mesh, order=2, dirichlet="axis|axis_ext", dirichlet_bbnd="GND")
print(f"  fes_pre.ndof = {fes_pre.ndof}", flush=True)
# Test GetDofNrs(NodeId) directly
from ngsolve.comp import NodeId as CompNodeId
from ngsolve.fem import NODE_TYPE as NT
try:
    # node 0 should map to DOF 0 (vertex), node 3 to DOF 3
    d0 = fes_pre.GetDofNrs(CompNodeId(NT.VERTEX, 0))
    d3 = fes_pre.GetDofNrs(CompNodeId(NT.VERTEX, 3))
    print(f"  GetDofNrs(VERTEX, 0) = {list(d0)}", flush=True)
    print(f"  GetDofNrs(VERTEX, 3) = {list(d3)}", flush=True)
    e0 = fes_pre.GetDofNrs(CompNodeId(NT.EDGE, 0))
    e7 = fes_pre.GetDofNrs(CompNodeId(NT.EDGE, 7))
    print(f"  GetDofNrs(EDGE, 0) = {list(e0)}", flush=True)
    print(f"  GetDofNrs(EDGE, 7) = {list(e7)}", flush=True)
except Exception as e:
    print(f"  Direct GetDofNrs(NodeId) FAIL: {type(e).__name__}: {e}", flush=True)

print("step 5: fes = Periodic(fes_pre)", flush=True)
try:
    fes = Periodic(fes_pre)
    print(f"  fes.ndof = {fes.ndof}", flush=True)
    n_before = sum(1 for d in fes_pre.FreeDofs() if d)
    n_after = sum(1 for d in fes.FreeDofs() if d)
    print(f"  free DOFs: {n_before} -> {n_after} (coupled: {n_before - n_after})",
          flush=True)
except Exception as e:
    print(f"  Periodic FAILED: {type(e).__name__}: {e}", flush=True)
    sys.exit(1)

print("step 6: assemble stiffness", flush=True)
mu0 = 4 * pi * 1e-7
mu_cf = CoefficientFunction(mu0)
sigma_cf = mesh.MaterialCF({}, default=0.0)

try:
    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager():
        a.Assemble()
    print(f"  K nzelems = {a.mat.nze}", flush=True)
except Exception as e:
    print(f"  Stiffness FAILED: {type(e).__name__}: {e}", flush=True)
    sys.exit(1)

print("OK", flush=True)
