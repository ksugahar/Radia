# -*- coding: utf-8 -*-
"""Lightweight checks for the promoted build123d pipeline API."""

from __future__ import annotations

import os
import sys


_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)


def test_build123d_pipeline_public_api_imports_without_geometry_runtime():
    from radia_mcp.build123d import pipeline

    assert pipeline.__all__ == ["run_pipeline", "run_pipeline_multi", "save_record"]
    for name in pipeline.__all__:
        assert callable(getattr(pipeline, name))

def test_pipeline_gmsh41_reader_maps_volume_physical_tags(tmp_path):
    from radia_mcp.build123d import pipeline

    msh = tmp_path / "one_tet_v41.msh"
    msh.write_text(
        """$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 0 1
10 0 0 0 1 1 1 1 7 0
$EndEntities
$Nodes
1 4 1 4
3 10 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
1 1 1 1
3 10 4 1
1 1 2 3 4
$EndElements
""",
        encoding="utf-8",
    )

    nodes, elements = pipeline._read_gmsh_nodes_elements(msh)

    assert nodes == [
        (1, 0.0, 0.0, 0.0),
        (2, 1.0, 0.0, 0.0),
        (3, 0.0, 1.0, 0.0),
        (4, 0.0, 0.0, 1.0),
    ]
    assert elements == [{
        "id": 1,
        "type": 4,
        "tags": [7, 10],
        "nodes": [1, 2, 3, 4],
    }]

    post = pipeline._stage_post_multi(msh, ["core"], tmp_path, "one_tet")
    text = (tmp_path / "one_tet_post.msh").read_text(encoding="utf-8")

    assert post["n_regions"] == 1
    assert post["regions"][0]["physical_tag"] == 7
    assert '"core"' in text
    assert "$ElementData" in text
