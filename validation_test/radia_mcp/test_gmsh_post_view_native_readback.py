"""VOL post views through Radia's exporter load cleanly in native Gmsh 5.

Guards: a checked two-region Netgen VOL converted by `convert_post_view`
(materials and boundaries) is published only after native Gmsh reads the
published bytes back with the selected groups and element count, no warning
and no unconnected node; the node count equals an independent NGSolve count
(all vertices for materials, the distinct vertices of the BND elements for
boundaries).
Assumes: the qualified Gmsh 5 runtime (RADIA_GMSH_PROVENANCE), Netgen OCC,
NGSolve and a Radia exporter that compacts boundary nodes (a759e448d or
later). An older exporter FAILS the boundary case: that is the regression
this lane exists to catch, not a skip.
"""
import hashlib
import json
import os
import tempfile
from pathlib import Path

import pytest

pytestmark = pytest.mark.skipif(not os.environ.get("RADIA_GMSH_PROVENANCE"),
                                reason="needs the qualified Gmsh 5 runtime (RADIA_GMSH_PROVENANCE)")


@pytest.fixture(scope="module")
def checked_vol(tmp_path_factory):
    occ = pytest.importorskip("netgen.occ")
    ng = pytest.importorskip("ngsolve")
    from radia_mcp.mesh_viewer.post_view import _load_exporter
    _load_exporter()  # fails loudly without the production MSH 4.2 exporter
    base = tmp_path_factory.mktemp("post_view")
    source_dir = base / "in"; source_dir.mkdir()
    scratch = base / "out"; scratch.mkdir()
    iron = occ.Box((0, 0, 0), (0.01, 0.01, 0.01)); iron.mat("iron"); iron.faces.name = "iron_wall"
    air = occ.Box((0.01, 0, 0), (0.02, 0.01, 0.01)); air.mat("air"); air.faces.name = "air_wall"
    source = source_dir / "two_region.vol"
    occ.OCCGeometry(occ.Glue([iron, air])).GenerateMesh(maxh=0.005).Save(str(source))
    mesh = ng.Mesh(str(source))
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    # The handoff is caller-asserted (documented contract): counts come from NGSolve itself.
    report = {"schema": "cubit-mesh-export.vol-check.v1", "passed": True,
              "mesh": {"dimension": mesh.dim, "curve_order": 1, "n_elements": mesh.ne, "n_points": mesh.nv},
              "labels": {"schema": "radia.vol-label-contract.v1", "passed": True,
                         "materials": sorted(set(mesh.GetMaterials())),
                         "boundaries": sorted(set(mesh.GetBoundaries()))}}
    handoff = base / "handoff.json"
    handoff.write_text(json.dumps({"schema": "radia.checked-vol-handoff.v1", "source_sha256": digest,
                                   "report": report}), encoding="utf-8")
    return source, scratch, handoff, digest, mesh, ng


@pytest.mark.parametrize("scope", ["materials", "boundaries"])
def test_checked_vol_post_view_loads_cleanly_in_native_gmsh(checked_vol, scope, monkeypatch):
    from radia_mcp.gmsh import readback as module
    from radia_mcp.mesh_viewer import post_view
    source, scratch, handoff, digest, mesh, ng = checked_vol
    monkeypatch.setenv("RADIA_MESH_VIEWER_SCRATCH", str(scratch))
    boundary = scope == "boundaries"
    dim = mesh.dim - int(boundary)
    elements = list(mesh.Elements(ng.BND if boundary else ng.VOL))
    names = mesh.GetBoundaries() if boundary else mesh.GetMaterials()
    # Independent node count from NGSolve: the distinct vertices the selected elements use.
    vertices = {v.nr for element in elements for v in element.vertices}
    calls = []
    original = module.require_native_readback
    monkeypatch.setattr(module, "require_native_readback",
                        lambda path, **kwargs: calls.append(kwargs) or original(path, **kwargs))
    operation = Path(tempfile.mkdtemp(prefix=".radia-post-receipt-", dir=scratch))
    result = post_view.convert_post_view(source, str(scratch / f"view_{scope}.msh"), str(handoff), digest,
                                         scope, staging_dir=operation / "staging")
    # The exporter's companion stays in staging; the readback read the MSH bytes only.
    assert (operation / "staging" / "view.msh.opt").exists()
    readback = result["native_readback"]
    assert calls == [{"expected_groups": {(dim, n) for n in names},
                      "expected_elements": {dim: len(elements)}}], calls
    assert readback["elements"] == {str(dim): len(elements)}, readback
    assert {tuple(g) for g in readback["physical_groups"]} == {(dim, n) for n in names}
    assert (readback["nodes"], readback["unreferenced_nodes"]) == (len(vertices), 0), (readback, len(vertices))
    assert readback["sha256"] == result["sha256"] == hashlib.sha256(Path(result["output"]).read_bytes()).hexdigest()
