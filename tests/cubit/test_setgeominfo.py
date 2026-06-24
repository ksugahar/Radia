"""Current Cubit 2025.12 export path for former SetGeomInfo coverage."""

from tests.cubit.cubit_202512_helpers import (
    export_netgen,
    init_cubit,
    load_ngsolve_mesh,
)


def test_hex_mesh_curves_through_export_netgen():
    cubit = init_cubit()
    cubit.cmd("reset")
    cubit.cmd("brick x 1 y 1 z 1")
    cubit.cmd("volume 1 scheme map")
    cubit.cmd("volume 1 size 0.25")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add hex all")
    cubit.cmd('block 1 name "domain"')
    cubit.cmd("block 2 add quad all in surface all")
    cubit.cmd('block 2 name "boundary"')

    mesh = load_ngsolve_mesh(export_netgen(cubit, "hex", order=2))
    assert mesh.ne == cubit.get_hex_count()
    assert mesh.nv > 0
