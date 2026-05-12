"""Check if NGSolve reorders triangle vertices vs Python-supplied order."""
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

# Try several triangles with deliberately non-canonical orderings.
cases = [
    ("CCW [0,1,2]",   [(1e-3, 0), (5e-3, 0), (3e-3, 2e-3)], [0, 1, 2]),
    ("CW  [0,2,1]",   [(1e-3, 0), (3e-3, 2e-3), (5e-3, 0)], [0, 1, 2]),
    ("any [3,4,5]",   [(2e-3, 1e-3), (6e-3, 1e-3), (4e-3, 4e-3)], [0, 1, 2]),
]
pids_all = []
for name, coords, _ in cases:
    pid_set = []
    for r, z in coords:
        pid_set.append(m.Add(MeshPoint(Pnt(r, z, 0))))
    pids_all.append(pid_set)
    m.Add(Element2D(1, pid_set))

# Add boundary edges around first triangle so the mesh is valid.
m.Add(Element1D([pids_all[0][0], pids_all[0][1]], index=1))
m.Add(Element1D([pids_all[0][1], pids_all[0][2]], index=1))
m.Add(Element1D([pids_all[0][2], pids_all[0][0]], index=1))

mesh = Mesh(m)
print(f"Mesh has {mesh.ne} triangles")
for k, el in enumerate(mesh.Elements()):
    name, _, pid_set_input = cases[k]
    verts_returned = list(el.vertices)
    # pids supplied: 1-based pnums from netgen.Add().
    pids_supplied = pids_all[k]
    print(f"\nTri {k}: {name}")
    print(f"  Supplied pids: {pids_supplied}")
    print(f"  ngel.vertices: {verts_returned}")
    # Compare ordering
    supplied_indices = [int(p.nr) - 1 for p in pids_supplied]  # 0-based
    returned_indices = [int(v.nr) for v in verts_returned]
    print(f"  Supplied 0-based: {supplied_indices}")
    print(f"  Returned 0-based: {returned_indices}")
    if supplied_indices == returned_indices:
        print(f"  -> ORDER PRESERVED")
    else:
        print(f"  -> ORDER REORDERED!")
