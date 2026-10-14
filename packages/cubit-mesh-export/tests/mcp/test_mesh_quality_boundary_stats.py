"""mesh_quality boundary statistics on high-order and mixed-element meshes.

Gmsh 5 changed the face API (faces of any size with ``faceSizes``); the
statistics must stay right on Gmsh 4 and 5.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("gmsh")
pytest.importorskip("numpy")

from cubit_mesh_export.mcp._support.mesh_quality import mesh_quality

def _mesh_and_truth(path: Path, kind: str) -> dict:
    """Write a .msh with Gmsh and count its boundary independently.

    The truth comes from the geometry: surfaces that bound exactly one
    volume are the boundary, so their 2D elements are the boundary faces
    and their nodes (high-order ones included) the boundary nodes.
    """
    import gmsh

    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(kind)
        if kind in ("tet_p1", "tet_p2"):
            gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
            gmsh.model.occ.synchronize()
            gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
            gmsh.model.mesh.generate(3)
            if kind == "tet_p2":
                gmsh.model.mesh.setOrder(2)
        elif kind == "hex_prism":
            # A quad patch and a triangle patch sharing an edge, extruded
            # in one layer: hexes next to prisms, sharing quad faces.
            a = gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
            b = gmsh.model.occ.addRectangle(1, 0, 0, 1, 1)
            gmsh.model.occ.fragment([(2, a)], [(2, b)])
            gmsh.model.occ.synchronize()
            quad, tri = sorted(s for _, s in gmsh.model.getEntities(2))
            gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
            gmsh.model.mesh.setRecombine(2, quad)
            gmsh.model.occ.extrude([(2, quad), (2, tri)], 0, 0, 1,
                                   numElements=[2], recombine=True)
            gmsh.model.occ.synchronize()
            gmsh.model.mesh.generate(3)
        else:
            raise ValueError(kind)
        gmsh.write(str(path))

        owners: dict[int, int] = {}
        for _, vol in gmsh.model.getEntities(3):
            for _, surf in gmsh.model.getBoundary([(3, vol)], combined=False,
                                                  oriented=False):
                owners[abs(surf)] = owners.get(abs(surf), 0) + 1
        outer = [s for s, n in owners.items() if n == 1]
        faces = 0
        nodes: set[int] = set()
        for s in outer:
            _, tags, _ = gmsh.model.mesh.getElements(2, s)
            faces += sum(len(t) for t in tags)
            node_tags, _, _ = gmsh.model.mesh.getNodes(2, s, includeBoundary=True)
            nodes.update(int(t) for t in node_tags)
        types = {int(t) for t in gmsh.model.mesh.getElementTypes(3)}
        return {"faces": faces, "nodes": len(nodes), "types": types}
    finally:
        gmsh.finalize()


@pytest.mark.parametrize("kind", ["tet_p1", "tet_p2", "hex_prism"])
def test_boundary_matches_geometry(tmp_path, kind):
    """Boundary faces and nodes against the geometry's own boundary.

    tet_p2 has six nodes per face and hex_prism has hex/prism faces in
    common; both were miscounted when faces were read as three or four
    nodes per corner count and matched only within one element type.
    """
    path = tmp_path / f"{kind}.msh"
    truth = _mesh_and_truth(path, kind)
    if kind == "hex_prism":
        assert {5, 6} <= truth["types"], truth  # hexahedra and prisms
    stats = mesh_quality(path)["mesh_stats"]
    assert stats["n_boundary_faces"] == truth["faces"], (stats, truth)
    assert stats["n_boundary_nodes"] == truth["nodes"], (stats, truth)
