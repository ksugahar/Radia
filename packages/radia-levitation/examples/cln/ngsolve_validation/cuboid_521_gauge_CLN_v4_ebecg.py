"""Historical field-recurrence experiment retired.

Only the independent geometry, mesh and linear-algebra helpers remain:
build_cuboid, classify_vertices, build_spanning_tree_interior.
There is no reduction experiment entry point or validation result in this file.

Historical method attribution: Nagamine.

Reference: cuboid_521_mathematica_reference.json (Nagamine high-precision).
"""
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, H1, BilinearForm, LinearForm, GridFunction, Preconditioner,
    CoefficientFunction, curl, grad, dx, x, y, z, BitArray, BND,
    Integrate, TaskManager, ngsglobals,
)
from ngsolve.krylovspace import CGSolver
from collections import deque
from math import pi
import json
import time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
V_box = ax * ay * az
N_STAGES = 12
H_COND = 0.10e-3   # start at 0.10 to verify EBE+CG matches v3, then push
ORDER = 3
CG_TOL = 1e-12
CG_MAXITER = 100000


def build_cuboid():
    box = Box(Pnt(0, 0, 0), Pnt(ax, ay, az))
    box.mat("conductor").bc("conductor_surface")
    box.maxh = H_COND
    return OCCGeometry(box)


def classify_vertices(mesh):
    boundary_v = set()
    for el in mesh.Elements(BND):
        for v in el.vertices:
            boundary_v.add(v.nr)
    return set(range(mesh.nv)) - boundary_v, boundary_v


def build_spanning_tree_interior(mesh):
    interior_v, _ = classify_vertices(mesh)
    nv = mesh.nv
    visited = [False] * nv
    for v in range(nv):
        if v not in interior_v:
            visited[v] = True
    tree_edges = []
    adj = [[] for _ in range(nv)]
    for ed in mesh.edges:
        v0, v1 = ed.vertices[0].nr, ed.vertices[1].nr
        if v0 in interior_v and v1 in interior_v:
            adj[v0].append((v1, ed.nr))
            adj[v1].append((v0, ed.nr))
    for start in interior_v:
        if visited[start]:
            continue
        visited[start] = True
        queue = deque([start])
        while queue:
            v = queue.popleft()
            for vn, edn in adj[v]:
                if not visited[vn]:
                    visited[vn] = True
                    tree_edges.append(edn)
                    queue.append(vn)
    return tree_edges







# CLN experiment entry point retired; independent helpers above are retained.
