"""Historical field-recurrence experiment retired.

Only the independent geometry, mesh and linear-algebra helpers remain:
build_spanning_tree, build_geo.
There is no reduction experiment entry point or validation result in this file.

Historical method attribution: Kameari.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"W:/00_CAE/Radia/01_GitHub/src/radia")
print("Starting...", flush=True)

from netgen.occ import Box, Sphere, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, Periodic, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, ngsglobals,
)
from kelvin_geometry import add_kelvin_exterior_domain
from kelvin_material import make_kelvin_nu_cf, NU_0
from collections import deque
from math import pi
import json, time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az
N_STAGES = 8
ORDER = 2
BONUS_INT = 8

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


def build_geo(R_K, h_air, h_cond=0.5e-3):
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.name = "conductor"
    cuboid.maxh = h_cond

    sphere_inner = Sphere(Pnt(0, 0, 0), R_K)
    for f in sphere_inner.faces:
        f.name = "kelvin_int"
    sphere_inner.maxh = h_air

    inner_air = sphere_inner - cuboid
    inner_air.name = "air"

    OFFSET = (2.5 * R_K, 0, 0)
    geo, info = add_kelvin_exterior_domain(
        [inner_air, cuboid], offset=OFFSET, R_K=R_K, inner_maxh=h_air)
    return OCCGeometry(geo), OFFSET







# CLN experiment entry point retired; independent helpers above are retained.
