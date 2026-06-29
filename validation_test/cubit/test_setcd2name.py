"""BBND edge-label propagation through export netgen."""

import os

from tests.cubit.cubit_202512_helpers import (
    export_netgen,
    init_cubit,
    load_ngsolve_mesh,
)


MODEL = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        "../../validation_test/induction_heating/cubit_panels_legacy/inductance_torus.cub5",
    )
)


def test_curve_sideset_name_reaches_ngsolve_bboundaries():
    cubit = init_cubit()
    cubit.cmd("reset")
    cubit.cmd(f'open "{MODEL.replace(os.sep, "/")}"')

    curve_ids = list(cubit.parse_cubit_list("curve", "all"))
    assert curve_ids
    cubit.cmd(f"sideset 100 add curve {curve_ids[0]}")
    cubit.cmd('sideset 100 name "test_edge"')

    mesh = load_ngsolve_mesh(export_netgen(cubit, "cd2", order=2))
    assert "test_edge" in tuple(mesh.GetBBoundaries())
