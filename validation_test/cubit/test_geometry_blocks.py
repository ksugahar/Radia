"""Geometry-block export tests for the Coreform Cubit 2025.12 plugin.

The public contract is Cubit's APREPRO command:

    cubit.cmd('export vol "model.vol" order N overwrite')
"""

from __future__ import annotations

from validation_test.cubit.cubit_202512_helpers import (
    export_vol,
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

    mesh = load_ngsolve_mesh(export_vol(cubit, "geometry_blocks", order=1))

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

    mesh = load_ngsolve_mesh(export_vol(cubit, "element_blocks", order=1))
    assert mesh.ne == cubit.get_tet_count()

