"""Mesh-origin orientation regression for GmshPostExport volume elements.

Vertex-listing orientation is MESH-ORIGIN dependent (measured
2026-08-06): netgen-generated meshes list tet corners NEGATIVELY
oriented, while Cubit-plugin `.vol` files list them POSITIVELY -- same
reference-corner correspondence, mirrored geometry.  The static corner
permutation introduced by the 2026-08-05 inversion fix was therefore
correct for netgen-origin meshes and inverted EVERY element of
cubit-origin meshes (Tet4 496/496 and Hex8 64/64 negative on the
reference sphere/brick).

GmshPostExport now selects the corner permutation PER ELEMENT from the
measured GMSH corner-frame determinant and fails loud when neither
orientation table yields a positive frame.  These tests lock both
orientations for every 3D element family using hand-built netgen meshes
(no Cubit needed): each written .msh must re-open in gmsh with ZERO
negative Jacobian points.
"""

import os
import sys

import pytest

ng = pytest.importorskip("ngsolve")
gmsh = pytest.importorskip("gmsh")

from netgen.csg import Pnt
from netgen.meshing import Element2D, Element3D, FaceDescriptor, MeshPoint
from netgen.meshing import Mesh as NgMesh
from ngsolve import Mesh

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from radia.gmsh_post_export import GmshPostExport


def _count_negative_jacobians(path):
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(path)
        negative = 0
        total = 0
        etypes, _etags, _ntags = gmsh.model.mesh.getElements(3)
        for et in etypes:
            pts, _w = gmsh.model.mesh.getIntegrationPoints(int(et), "Gauss2")
            _jac, det, _ = gmsh.model.mesh.getJacobians(int(et), pts)
            negative += sum(1 for d in det if d <= 0)
            total += len(det)
        return negative, total
    finally:
        gmsh.finalize()


def _mesh_from_elements(points, elements):
    m = NgMesh(3)
    m.Add(FaceDescriptor(bc=1, domin=1))
    pnums = [m.Add(MeshPoint(Pnt(*p))) for p in points]
    for conn in elements:
        m.Add(Element3D(1, [pnums[i] for i in conn]))
    m.SetMaterial(1, "mat")
    return Mesh(m)


# Base (negatively-listed, netgen-native) single elements; the mirrored
# variant swaps two vertices to flip the listing orientation.
_CASES = {
    "tet": {
        "points": [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
        "neg": [0, 2, 1, 3],
        "pos": [0, 1, 2, 3],
    },
    "hex": {
        "points": [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                   (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)],
        "neg": [0, 3, 2, 1, 4, 7, 6, 5],
        "pos": [0, 1, 2, 3, 4, 5, 6, 7],
    },
    "prism": {
        "points": [(0, 0, 0), (1, 0, 0), (0, 1, 0),
                   (0, 0, 1), (1, 0, 1), (0, 1, 1)],
        "neg": [0, 2, 1, 3, 5, 4],
        "pos": [0, 1, 2, 3, 4, 5],
    },
    "pyramid": {
        "points": [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                   (0.5, 0.5, 1)],
        "neg": [0, 3, 2, 1, 4],
        "pos": [0, 1, 2, 3, 4],
    },
}

_TMP = os.environ.get("TEMP", "C:\\temp")


@pytest.mark.parametrize("family", sorted(_CASES))
@pytest.mark.parametrize("orientation", ["neg", "pos"])
def test_single_element_both_orientations(family, orientation):
    case = _CASES[family]
    mesh = _mesh_from_elements(case["points"], [case[orientation]])
    out = os.path.join(_TMP, f"orient_{family}_{orientation}.msh")
    GmshPostExport(mesh).write(out)
    negative, total = _count_negative_jacobians(out)
    assert total > 0
    assert negative == 0, (family, orientation, negative, total)


def test_mixed_orientations_in_one_mesh():
    """Two tets of OPPOSITE listing orientation in the same mesh: the
    per-element selection must handle both within one export."""
    points = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, -1)]
    elements = [
        [0, 2, 1, 3],   # negatively listed (netgen-native)
        [0, 1, 2, 4],   # positively listed w.r.t. its own frame
    ]
    mesh = _mesh_from_elements(points, elements)
    out = os.path.join(_TMP, "orient_mixed_tets.msh")
    GmshPostExport(mesh).write(out)
    negative, total = _count_negative_jacobians(out)
    assert total > 0
    assert negative == 0, (negative, total)
