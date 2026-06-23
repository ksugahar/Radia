"""Mixed element Netgen export through Cubit 2025.12."""

from tests.cubit.cubit_202512_helpers import (
    export_netgen,
    init_cubit,
    load_ngsolve_mesh,
)


def test_mixed_tet_hex_blocks_export_and_load():
    cubit = init_cubit()
    cubit.cmd("reset")
    cubit.cmd("brick x 1 y 1 z 1")
    cubit.cmd("volume 1 move -0.75 0 0")
    cubit.cmd("brick x 1 y 1 z 1")
    cubit.cmd("volume 2 move 0.75 0 0")

    cubit.cmd("volume 1 scheme map")
    cubit.cmd("volume 1 size 0.5")
    cubit.cmd("mesh volume 1")
    cubit.cmd("volume 2 scheme tetmesh")
    cubit.cmd("volume 2 size 0.5")
    cubit.cmd("mesh volume 2")

    cubit.cmd("block 1 add hex all in volume 1")
    cubit.cmd('block 1 name "hex_region"')
    cubit.cmd("block 2 add tet all in volume 2")
    cubit.cmd('block 2 name "tet_region"')
    cubit.cmd("block 3 add face all in surface all")
    cubit.cmd('block 3 name "boundary"')

    expected = len(cubit.get_block_hexes(1)) + len(cubit.get_block_tets(2))
    mesh = load_ngsolve_mesh(export_netgen(cubit, "mixed", order=2))
    assert mesh.ne == expected
    assert set(mesh.GetMaterials()) == {"hex_region", "tet_region"}
