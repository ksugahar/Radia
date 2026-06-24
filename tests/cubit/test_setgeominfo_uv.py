"""Curved cylinder accuracy through Cubit 2025.12 export netgen."""

import math

from ngsolve import BND, CF, Integrate, TaskManager

from tests.cubit.cubit_202512_helpers import (
    export_netgen,
    init_cubit,
    load_ngsolve_mesh,
)


def test_export_netgen_uses_geometry_info_for_curved_cylinder():
    cubit = init_cubit()
    radius = 0.5
    height = 2.0
    expected_volume = math.pi * radius * radius * height
    expected_area = 2 * math.pi * radius * height + 2 * math.pi * radius * radius

    cubit.cmd("reset")
    cubit.cmd(f"create cylinder height {height} radius {radius}")
    cubit.cmd("volume all scheme tetmesh")
    cubit.cmd("volume all size 0.2")
    cubit.cmd("mesh volume all")
    cubit.cmd("block 1 add tet all")
    cubit.cmd('block 1 name "domain"')
    cubit.cmd("block 2 add tri all")
    cubit.cmd('block 2 name "boundary"')

    with TaskManager():
        mesh1 = load_ngsolve_mesh(export_netgen(cubit, "cylinder_uv", order=1))
        vol1 = Integrate(CF(1), mesh1)
        mesh3 = load_ngsolve_mesh(export_netgen(cubit, "cylinder_uv", order=3))
        vol3 = Integrate(CF(1), mesh3)
        area3 = Integrate(CF(1), mesh3, VOL_or_BND=BND)

    assert abs(vol3 - expected_volume) < abs(vol1 - expected_volume)
    assert abs(area3 - expected_area) / expected_area < 0.03
