"""Historical field-recurrence experiment retired.

Only the independent geometry, mesh and linear-algebra helpers remain:
build_spanning_tree, build_geometry.
There is no reduction experiment entry point or validation result in this file.

Historical method attribution: Kameari.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from netgen.occ import Box, Pnt, OCCGeometry, Glue
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, ngsglobals, InnerProduct,
)
from collections import deque
import numpy as np
from math import pi
import json, time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7

ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

AIR_SCALE = 5
H_COND = 0.30e-3
H_AIR = 1.5e-3
ORDER = 2
N_STAGES = 12

R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))


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


def build_geometry():
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.mat("conductor").bc("conductor_surface")
    cuboid.maxh = H_COND
    Bx, By, Bz = AIR_SCALE * ax / 2, AIR_SCALE * ay / 2, AIR_SCALE * az / 2
    air = Box(Pnt(-Bx, -By, -Bz), Pnt(Bx, By, Bz))
    air.mat("air").bc("outer_box")
    air.maxh = H_AIR
    return OCCGeometry(Glue([air - cuboid, cuboid]))







# Historical method attributions retained: Tanimoto, Kameari.

# CLN experiment entry point retired; independent helpers above are retained.
