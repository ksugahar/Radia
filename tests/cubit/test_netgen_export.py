"""Netgen export tests for the Coreform Cubit 2025.12 command path."""

from ngsolve import H1

from tests.cubit.cubit_202512_helpers import (
    export_netgen,
    init_cubit,
    load_ngsolve_mesh,
)


def test_basic_tet_export_loads_in_ngsolve():
    cubit = init_cubit()
    cubit.cmd("reset")
    cubit.cmd("create brick x 1 y 1 z 1")
    cubit.cmd("volume 1 scheme tetmesh")
    cubit.cmd("volume 1 size 0.5")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add tet all")
    cubit.cmd('block 1 name "solid"')
    cubit.cmd("block 2 add tri all in surface all")
    cubit.cmd('block 2 name "boundary"')

    mesh = load_ngsolve_mesh(export_netgen(cubit, "brick", order=1))
    assert mesh.ne == len(cubit.get_block_tets(1))
    assert tuple(mesh.GetMaterials()) == ("solid",)
    assert "boundary" in tuple(mesh.GetBoundaries())
    assert H1(mesh, order=1).ndof > 0


def test_multiple_material_blocks_survive_export():
    cubit = init_cubit()
    cubit.cmd("reset")
    cubit.cmd("create brick x 1 y 1 z 1")
    cubit.cmd("volume 1 move -0.75 0 0")
    cubit.cmd("create brick x 1 y 1 z 1")
    cubit.cmd("volume 2 move 0.75 0 0")
    cubit.cmd("volume all scheme tetmesh")
    cubit.cmd("volume all size 0.5")
    cubit.cmd("mesh volume all")
    cubit.cmd("block 1 add tet all in volume 1")
    cubit.cmd('block 1 name "region1"')
    cubit.cmd("block 2 add tet all in volume 2")
    cubit.cmd('block 2 name "region2"')
    cubit.cmd("block 3 add tri all in surface all")
    cubit.cmd('block 3 name "boundary"')

    mesh = load_ngsolve_mesh(export_netgen(cubit, "two_regions", order=1))
    assert set(mesh.GetMaterials()) == {"region1", "region2"}
