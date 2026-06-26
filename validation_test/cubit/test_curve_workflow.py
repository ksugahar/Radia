"""Cubit 2025.12 -> export netgen -> NGSolve high-order workflow."""

import math

from ngsolve import BND, CF, Integrate, TaskManager

from tests.cubit.cubit_202512_helpers import (
    export_netgen,
    init_cubit,
    load_ngsolve_mesh,
)


def test_cylinder_export_netgen_orders_have_reasonable_geometry():
    cubit = init_cubit()
    radius = 0.5
    height = 2.0
    expected_area = 2 * math.pi * radius * height + 2 * math.pi * radius * radius
    expected_volume = math.pi * radius * radius * height

    cubit.cmd("reset")
    cubit.cmd(f"create cylinder height {height} radius {radius}")
    cubit.cmd("volume all scheme tetmesh")
    cubit.cmd("volume all size 0.2")
    cubit.cmd("mesh volume all")
    cubit.cmd("block 1 add tet all")
    cubit.cmd('block 1 name "domain"')
    cubit.cmd("block 2 add tri all")
    cubit.cmd('block 2 name "boundary"')

    errors = {}
    for order in (1, 2, 3):
        mesh = load_ngsolve_mesh(export_netgen(cubit, "cylinder", order=order))
        with TaskManager():
            area = Integrate(CF(1), mesh, VOL_or_BND=BND)
            volume = Integrate(CF(1), mesh)
        errors[order] = abs(volume - expected_volume) / expected_volume
        assert area > 0

    assert errors[2] < errors[1]
    assert errors[3] <= errors[2] * 1.5
    assert errors[3] < 0.02
