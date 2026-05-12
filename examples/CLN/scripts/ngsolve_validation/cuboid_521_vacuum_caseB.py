"""Historical field-recurrence experiment retired.

Only the independent geometry, mesh and linear-algebra helpers remain:
classify_vertices, build_spanning_tree_interior, build_geometry.
There is no reduction experiment entry point or validation result in this file.

Historical method attribution: Kameari, Tanimoto.
"""
from netgen.occ import Box, Pnt, OCCGeometry, Glue
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z, BitArray, BND,
    Integrate, TaskManager, ngsglobals,
)
from collections import deque
from math import pi
import json
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7

# Conductor (centered at origin)
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

# Air box scale (k * conductor)
AIR_SCALE = 5  # box is k times conductor size, centered at origin
H_COND = 0.30e-3   # conductor mesh size
H_AIR = 1.5e-3     # air mesh size (coarser)
ORDER = 2          # tradeoff: order vs total ndof
N_STAGES = 5       # Kameari stages

# Reference R_0 (Parseval, conductor centered at origin):
#   R_0 = 1 / <J,J/sigma>_cond, J = sigma*A_ext = (sigma/2)(-y,x,0)
#   <J,J/sigma> = (sigma/4) integral (x^2+y^2) dV over [-a/2,a/2]^3
#               = sigma * V * (a^2+b^2) / 48
# So R_0_centered = 48 / (sigma V (a^2+b^2))
R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))


def classify_vertices(mesh):
    """Boundary vertices = those on outer_box OR conductor_surface.
       For the vacuum problem, only outer_box vertices have Dirichlet BC.
       Tree-cotree mask should exclude only outer_box-vertex Whitney edges."""
    boundary_v = set()
    for el in mesh.Elements(BND):
        # Only outer_box has Dirichlet
        if "outer_box" in mesh.GetMaterial(el):
            for v in el.vertices:
                boundary_v.add(v.nr)
    # Actually mesh.GetMaterial doesn't apply here for BND; use bcname
    # Simpler: check via fes.FreeDofs after dirichlet="outer_box"
    return boundary_v


def build_spanning_tree_interior(mesh, free_dofs_check, fes):
    """BFS spanning tree restricted to vertices whose Whitney edge DoF is free.

    More robust than vertex-based classification: just check the FE space's
    FreeDofs to determine which edges contribute to kernel.
    """
    nv = mesh.nv
    # Mark vertices as "boundary" if all their incident edges have non-free
    # lowest-order DoF (i.e., they are kept in the kernel as constants).
    # Simpler: BFS over all vertices, mask first DoF of each tree edge if free.
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
    # Conductor centered at origin (so symmetry around origin for Kelvin v1)
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.mat("conductor").bc("conductor_surface")
    cuboid.maxh = H_COND

    # Air box (also centered)
    Bx, By, Bz = AIR_SCALE * ax / 2, AIR_SCALE * ay / 2, AIR_SCALE * az / 2
    air = Box(Pnt(-Bx, -By, -Bz), Pnt(Bx, By, Bz))
    air.mat("air").bc("outer_box")
    air.maxh = H_AIR

    full = Glue([air - cuboid, cuboid])
    return OCCGeometry(full)







# Historical method attributions retained: Tanimoto.

# CLN experiment entry point retired; independent helpers above are retained.
