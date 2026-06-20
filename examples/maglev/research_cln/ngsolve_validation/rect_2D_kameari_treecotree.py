"""Historical field-recurrence experiment retired.

Only the independent geometry, mesh and linear-algebra helpers remain:
build_slab, build_spanning_tree.
There is no reduction experiment entry point or validation result in this file.

Historical method attribution: Tanimoto.
"""
import json
from collections import deque
from math import pi
from pathlib import Path

from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, BitArray,
    Integrate, TaskManager, ngsglobals,
)

mu0 = 4 * pi * 1e-7
sigma = 5.8e7
ax, ay = 5e-3, 2e-3
DZ = 0.1e-3        # thin slab thickness in z (must be small for quasi-2D)
H_MESH = 0.10e-3
ORDER = 2
N_STAGES = 12


def build_slab():
    box = Box(Pnt(0, 0, 0), Pnt(ax, ay, DZ))
    box.mat("conductor")
    # Name only x and y faces so we Dirichlet only those
    for face in box.faces:
        n = face.center
        # Determine which face by checking normal-aligned coordinate extreme
        if abs(n.x - 0) < 1e-9:
            face.bc("xleft")
        elif abs(n.x - ax) < 1e-9:
            face.bc("xright")
        elif abs(n.y - 0) < 1e-9:
            face.bc("yleft")
        elif abs(n.y - ay) < 1e-9:
            face.bc("yright")
        elif abs(n.z - 0) < 1e-9:
            face.bc("zbot")
        elif abs(n.z - DZ) < 1e-9:
            face.bc("ztop")
    box.maxh = H_MESH
    return OCCGeometry(box)


def build_spanning_tree(mesh):
    nv = mesh.nv
    visited = [False] * nv
    tree_edges = []
    adj = [[] for _ in range(nv)]
    for ed in mesh.edges:
        v0, v1 = ed.vertices[0].nr, ed.vertices[1].nr
        adj[v0].append((v1, ed.nr))
        adj[v1].append((v0, ed.nr))
    visited[0] = True
    queue = deque([0])
    while queue:
        v = queue.popleft()
        for vn, edn in adj[v]:
            if not visited[vn]:
                visited[vn] = True
                tree_edges.append(edn)
                queue.append(vn)
    return tree_edges







# Historical method attributions retained: Tanimoto.

# CLN experiment entry point retired; independent helpers above are retained.
