"""Fast native guard for the retired mesh-less nonlinear relaxation path."""

import pytest


@pytest.mark.basic
def test_meshless_relaxation_rejected_in_favor_of_hdiv_vim(radia_clean):
    """Mesh-less soft iron must fail loudly instead of approximating FEM."""
    rad = radia_clean
    material = rad.MatSatIsoFrm(
        [[1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759]]
    )
    elements = []
    for i in range(10):
        for j in range(10):
            magnet = rad.magnet_box(
                [-0.050 + i * 0.010, -0.050 + j * 0.010, 0],
                [0.008, 0.008, 0.010],
                [0, 0, 0],
            )
            rad.MatApl(magnet, material)
            elements.append(magnet)

    with pytest.raises(RuntimeError, match="[Mm]esh-less soft iron"):
        rad.Solve(rad.ObjCnt(elements), 0.0001, 1000, 0)
