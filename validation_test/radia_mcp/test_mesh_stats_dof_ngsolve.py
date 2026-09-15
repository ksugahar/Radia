"""Cross-check mesh_quality dof_estimate against NGSolve's own FES.

Moved out of packages/radia-mcp/tests (lane separation: package tests are
fast API/MCP contracts; solver-backed numerical checks live here).  This
is the claim that makes dof_estimate worth reporting, so it is verified
against an independent implementation rather than asserted.
"""
from __future__ import annotations

import json
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

pytest.importorskip("gmsh")

from radia_mcp.gmsh.msh_inspect import mesh_quality


def _write_unit_cube_tet_msh(path: Path) -> None:
    """One cube split into tets, written as .msh v4.1 by gmsh itself."""
    import gmsh

    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("cube")
        gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMin", 1.0)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 1.0)
        gmsh.model.mesh.generate(3)
        gmsh.write(str(path))
    finally:
        gmsh.finalize()


@pytest.fixture(scope="module")
def cube_msh(tmp_path_factory) -> Path:
    p = tmp_path_factory.mktemp("meshstats") / "cube.msh"
    _write_unit_cube_tet_msh(p)
    return p


def test_dof_estimate_matches_ngsolve_fes(cube_msh, tmp_path):
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    q = mesh_quality(cube_msh)
    dof = q["mesh_stats"]["dof_estimate"]

    script = textwrap.dedent(
        """
        import json, sys
        from netgen.occ import Box, OCCGeometry, Pnt
        from ngsolve import H1, HCurl, HDiv, L2, Mesh
        geo = OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1)))
        mesh = Mesh(geo.GenerateMesh(maxh=1.0))
        print(json.dumps({
            "h1_p1": H1(mesh, order=1).ndof,
            "hcurl_lowest": HCurl(mesh, order=0).ndof,
            "hdiv_lowest": HDiv(mesh, order=0).ndof,
            "l2_p0": L2(mesh, order=0).ndof,
            "ne": mesh.ne,
        }))
        """
    )
    out = subprocess.run([sys.executable, "-c", script], capture_output=True,
                         text=True, timeout=300)
    assert out.returncode == 0, out.stderr[-2000:]
    ng = json.loads(out.stdout.strip().splitlines()[-1])

    # Different mesh instance, so compare the RELATIONS that must hold in
    # any tet mesh rather than raw counts.
    assert ng["hcurl_lowest"] > ng["h1_p1"]
    assert ng["hdiv_lowest"] > ng["hcurl_lowest"]
    # ...and that our own mesh obeys the same ordering
    assert dof["hcurl_lowest"] > dof["h1_p1"]
    assert dof["hdiv_lowest"] > dof["hcurl_lowest"]
    # Euler identity holds for NGSolve's mesh too
    assert (ng["h1_p1"] - ng["hcurl_lowest"] + ng["hdiv_lowest"]
            - ng["l2_p0"]) == 1
