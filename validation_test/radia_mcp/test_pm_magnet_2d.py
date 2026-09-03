"""2D planar magnetostatic A_z with a PERMANENT-MAGNET (magnetization) source -- the PM extension
of solve_magnetostatic_az (currents only). Gated against the trusted in-repo closed form
radia.analytical_formulas.rect_magnet_2d (A_z, B of a 2D uniformly magnetised bar).
"""
import math
import os
import sys

import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
_CORE = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "src"))
for p in (_SRC, _CORE):
    if p not in sys.path:
        sys.path.insert(0, p)

from radia_mcp.radia_ngsolve.scalar_fem2d import solve_magnetostatic_az_magnet
rmag = pytest.importorskip("radia.analytical_formulas.rect_magnet_2d")   # trusted in-repo gate

MU0 = 4e-7 * math.pi


def _solve(a_h=0.1, b_h=0.1, M=(0.0, None), L=3.0, order=3):
    My = M[1] if M[1] is not None else 1.0 / MU0
    geo = SplineGeometry()
    geo.AddRectangle((-L, -L), (L, L), bc="far", leftdomain=1, rightdomain=0)
    geo.AddRectangle((-a_h, -b_h), (a_h, b_h), leftdomain=2, rightdomain=1)
    geo.SetMaterial(1, "air"); geo.SetMaterial(2, "mag")
    geo.SetDomainMaxH(2, 0.012)
    mesh = ng.Mesh(geo.GenerateMesh(maxh=0.4))
    gfu = solve_magnetostatic_az_magnet(mesh, 1.0 / MU0, {"mag": (M[0], My)}, {"far": 0.0}, order=order)
    B = ng.CF((ng.grad(gfu)[1], -ng.grad(gfu)[0]))
    return mesh, B, (a_h, b_h, (M[0], My))


def test_magnet_center_internal_field():
    import numpy as np
    mesh, B, (a_h, b_h, M) = _solve()
    bx, by = B(mesh(0.0, 0.0))
    ban = rmag.rect_magnet_2d_B(np.array([0.0, 0.0]), a_h, b_h, M)
    assert by == pytest.approx(float(ban[1]), rel=3e-3)        # square 2D magnet center By = mu0 M/2
    assert abs(bx) < 1e-2 * abs(by)                            # symmetry: Bx ~ 0


def test_magnet_exterior_field_matches_closed_form():
    import numpy as np
    mesh, B, (a_h, b_h, M) = _solve()
    for px, py in [(0.0, 0.25), (0.25, 0.0), (0.0, -0.3)]:
        bx, by = B(mesh(px, py))
        ban = rmag.rect_magnet_2d_B(np.array([px, py]), a_h, b_h, M)
        den = max(abs(ban[0]), abs(ban[1]))
        assert max(abs(bx - ban[0]), abs(by - ban[1])) / den < 2e-2   # finite box ~1%


def test_field_scales_with_magnetization():
    # B is linear in M: doubling M doubles the center field
    mesh1, B1, _ = _solve(M=(0.0, 1.0 / MU0))
    mesh2, B2, _ = _solve(M=(0.0, 2.0 / MU0))
    by1 = B1(mesh1(0.0, 0.0))[1]
    by2 = B2(mesh2(0.0, 0.0))[1]
    assert by2 / by1 == pytest.approx(2.0, rel=1e-3)
