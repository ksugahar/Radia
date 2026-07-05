"""Small .vol-backed mesh helpers for axifem tests.

The production route for axifem is a Netgen ``.vol`` file loaded by NGSolve.
These helpers intentionally save the generated Netgen mesh to ``.vol`` and
read it back with ``ngsolve.Mesh(path)`` before returning it to tests.
"""

from __future__ import annotations

import os
from pathlib import Path
import uuid

import numpy as np
from netgen.meshing import Element1D, Element2D, FaceDescriptor, Mesh as NgMesh
from netgen.meshing import MeshPoint, Pnt
from ngsolve import Mesh


def _vol_tmp_dir() -> Path:
    if os.name == "nt":
        base = Path(os.environ.get("RADIA_AXIFEM_TEST_TMP", r"C:\temp\radia_axifem_tests"))
    else:
        base = Path(os.environ.get("RADIA_AXIFEM_TEST_TMP", "/tmp/radia_axifem_tests"))
    base.mkdir(parents=True, exist_ok=True)
    return base


def reload_via_vol(mesh: Mesh, stem: str) -> Mesh:
    """Save an NGSolve mesh as Netgen .vol and reload it through Mesh(path)."""

    path = _vol_tmp_dir() / f"{stem}_{uuid.uuid4().hex}.vol"
    mesh.ngmesh.Save(str(path))
    return Mesh(str(path))


def structured_rect_vol_mesh(
    ra: float,
    rb: float,
    za: float,
    zb: float,
    *,
    nx: int = 1,
    ny: int = 1,
    quads: bool = True,
    mapping=None,
    stem: str = "axifem_rect",
) -> Mesh:
    """Build a structured rectangle/parallelogram .vol and reload it."""

    if mapping is None:
        mapping = lambda xi, eta: (ra + (rb - ra) * xi, za + (zb - za) * eta)

    ngmesh = NgMesh()
    ngmesh.dim = 2
    ngmesh.SetMaterial(1, "domain")

    ngmesh.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    ngmesh.Add(FaceDescriptor(surfnr=2, domin=0, bc=2))
    ngmesh.Add(FaceDescriptor(surfnr=3, domin=0, bc=3))
    ngmesh.Add(FaceDescriptor(surfnr=4, domin=0, bc=4))
    ngmesh.SetBCName(0, "bottom")
    ngmesh.SetBCName(1, "right")
    ngmesh.SetBCName(2, "top")
    ngmesh.SetBCName(3, "left")

    pids = np.empty((ny + 1, nx + 1), dtype=object)
    for j in range(ny + 1):
        eta = j / ny
        for i in range(nx + 1):
            xi = i / nx
            r, z = mapping(xi, eta)
            pids[j, i] = ngmesh.Add(MeshPoint(Pnt(r, z, 0.0)))

    for j in range(ny):
        for i in range(nx):
            if quads:
                ngmesh.Add(Element2D(1, [
                    pids[j, i],
                    pids[j, i + 1],
                    pids[j + 1, i + 1],
                    pids[j + 1, i],
                ]))
            else:
                ngmesh.Add(Element2D(1, [
                    pids[j, i],
                    pids[j, i + 1],
                    pids[j + 1, i + 1],
                ]))
                ngmesh.Add(Element2D(1, [
                    pids[j, i],
                    pids[j + 1, i + 1],
                    pids[j + 1, i],
                ]))

    for i in range(nx):
        ngmesh.Add(Element1D([pids[0, i], pids[0, i + 1]], index=1))
        ngmesh.Add(Element1D([pids[ny, i], pids[ny, i + 1]], index=3))
    for j in range(ny):
        ngmesh.Add(Element1D([pids[j, nx], pids[j + 1, nx]], index=2))
        ngmesh.Add(Element1D([pids[j, 0], pids[j + 1, 0]], index=4))

    path = _vol_tmp_dir() / f"{stem}_{uuid.uuid4().hex}.vol"
    ngmesh.Save(str(path))
    return Mesh(str(path))
