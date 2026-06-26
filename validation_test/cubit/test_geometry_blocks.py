"""Geometry-block export tests for the Coreform Cubit 2025.12 plugin.

The public contract is Cubit's APREPRO command:

    cubit.cmd('export netgen "model.vol" order N overwrite')

The retired Python helper ``extract_mesh_data`` is intentionally not used here.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

from tests.cubit.cubit_202512_helpers import (
    export_netgen,
    init_cubit,
    load_ngsolve_mesh,
)


def _mesh_brick(cubit) -> None:
    cubit.cmd("reset")
    cubit.cmd("create brick x 1 y 1 z 1")
    cubit.cmd("volume 1 scheme tetmesh")
    cubit.cmd("volume 1 size 0.5")
    cubit.cmd("mesh volume 1")


def test_geometry_volume_and_surface_blocks_export_to_netgen():
    """Volume/surface geometry blocks survive the 2025.12 Netgen export."""
    cubit = init_cubit()
    _mesh_brick(cubit)

    expected_tets = cubit.get_tet_count()
    cubit.cmd("block 1 add volume 1")
    cubit.cmd('block 1 name "solid"')
    cubit.cmd("block 2 add surface all")
    cubit.cmd('block 2 name "boundary"')

    mesh = load_ngsolve_mesh(export_netgen(cubit, "geometry_blocks", order=1))

    assert mesh.ne == expected_tets
    assert "solid" in mesh.GetMaterials()
    assert "boundary" in mesh.GetBoundaries()


def test_mesh_element_blocks_do_not_cross_contaminate():
    """Cubit still keeps volume and boundary element blocks separated."""
    cubit = init_cubit()
    _mesh_brick(cubit)

    cubit.cmd("block 1 add tet all")
    cubit.cmd('block 1 name "solid"')
    cubit.cmd("block 2 add tri all")
    cubit.cmd('block 2 name "boundary"')

    assert len(cubit.get_block_tets(1)) > 0
    assert len(cubit.get_block_tris(1)) == 0
    assert len(cubit.get_block_tets(2)) == 0
    assert len(cubit.get_block_tris(2)) > 0

    mesh = load_ngsolve_mesh(export_netgen(cubit, "element_blocks", order=1))
    assert mesh.ne == cubit.get_tet_count()


def test_gmsh_geometry_volume_block_exports_elements():
    """Geometry-based volume blocks are also accepted by Cubit's Gmsh export."""
    cubit = init_cubit()
    _mesh_brick(cubit)

    expected_tets = cubit.get_tet_count()
    cubit.cmd("block 1 add volume 1")
    cubit.cmd('block 1 name "solid"')

    out = Path(tempfile.mkdtemp(prefix="radia_cubit_gmsh_")) / "geometry.msh"
    cmd_path = str(out).replace("\\", "/")
    cubit.cmd(f'export gmsh "{cmd_path}" overwrite')

    content = out.read_text(encoding="utf-8", errors="replace")
    tet_count = 0
    in_elements = False
    for line in content.splitlines():
        if line.strip() == "$Elements":
            in_elements = True
            continue
        if line.strip() == "$EndElements":
            break
        if not in_elements or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 2 and parts[1] in {"4", "11"}:
            tet_count += 1

    assert tet_count == expected_tets
