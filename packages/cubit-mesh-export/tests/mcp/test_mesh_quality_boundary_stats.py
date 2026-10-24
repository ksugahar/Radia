"""Guards: Gmsh 5 post-processing counts high-order/shared boundary faces.
Assumes: stored explicit MSH 4.2 connectivity; no Gmsh meshing or solver input.
"""
from pathlib import Path
import pytest
from cubit_mesh_export.mcp._support.mesh_quality import mesh_quality

FIXTURES = Path(__file__).parent / "fixtures" / "msh42"

@pytest.mark.parametrize("kind,faces,nodes", [("tet_p1",4,4),("tet_p2",4,10),("hex_prism",9,10)])
def test_boundary_matches_geometry(kind, faces, nodes):
    # Four tetrahedron faces; hex+prism share one quad (6+5-2 faces).
    stats = mesh_quality(FIXTURES / f"{kind}.msh")["mesh_stats"]
    assert stats["n_boundary_faces"] == faces, stats
    assert stats["n_boundary_nodes"] == nodes, stats
