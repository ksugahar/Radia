"""Current first-order Netgen export contract for Cubit 2025.12."""

from tests.cubit.cubit_202512_helpers import (
    export_netgen,
    init_cubit,
    load_ngsolve_mesh,
)


def test_order1_export_uses_linear_topology_even_after_tetra10_assignment():
    cubit = init_cubit()
    cubit.cmd("reset")
    cubit.cmd("create sphere radius 1")
    cubit.cmd("volume 1 scheme tetmesh")
    cubit.cmd("volume 1 size 0.5")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add tet all in volume 1")
    cubit.cmd('block 1 name "sphere"')
    cubit.cmd("block 1 element type tetra10")
    cubit.cmd("block 2 add tri all in surface all")
    cubit.cmd('block 2 name "boundary"')
    cubit.cmd("block 2 element type tri6")

    tet_id = cubit.get_block_tets(1)[0]
    assert len(cubit.get_connectivity("tet", tet_id)) == 4
    assert len(cubit.get_expanded_connectivity("tet", tet_id)) == 10

    mesh = load_ngsolve_mesh(export_netgen(cubit, "sphere_linear", order=1))
    assert mesh.ne == len(cubit.get_block_tets(1))


def test_higher_order_export_loads_as_ngsolve_mesh():
    cubit = init_cubit()
    cubit.cmd("reset")
    cubit.cmd("create sphere radius 1")
    cubit.cmd("volume 1 scheme tetmesh")
    cubit.cmd("volume 1 size 0.5")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add tet all in volume 1")
    cubit.cmd('block 1 name "sphere"')
    cubit.cmd("block 2 add tri all in surface all")
    cubit.cmd('block 2 name "boundary"')

    for order in (1, 2, 3):
        mesh = load_ngsolve_mesh(export_netgen(cubit, "sphere", order=order))
        assert mesh.ne == len(cubit.get_block_tets(1))
        assert mesh.nv > 0
