"""Golden test: 3D HCurl transverse (m=1) eddy polarizability of a conducting
ellipsoid, via CompactAMS + COCR on a graded fine-air-shell mesh.

Locks docs/maglev/demos/ellipsoid/ellipsoid_alpha_tensor_3d.py:
  - the 3D alpha_xx and alpha_zz reproduce the analytic sphere 4 pi a^3 G(x)
    (so the transverse m=1 machinery is quantitatively correct), and
  - alpha_xx == alpha_zz (isotropy) to < 2%.

The key to getting the magnitude right is a fine AIR shell around the body
(Re[alpha] is set by the reaction dipole field just outside the conductor);
order-1 suffices.  Needs ngsolve + radia (CompactAMS); skipped otherwise.
One 3D solve per direction (~30-60 s each).
"""
import math
import os
import sys

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
pytest.importorskip("radia.sparsesolv_ngsolve")

_HERE = os.path.dirname(os.path.abspath(__file__))
# CLN was absorbed into radia.maglev ->
# docs/maglev/demos/{sphere,ellipsoid}. Add both so sphere + ellipsoid modules resolve.
_LEV = os.path.join(_HERE, "..", "docs", "maglev", "demos")
sys.path.insert(0, os.path.join(_LEV, "sphere"))
sys.path.insert(0, os.path.join(_LEV, "ellipsoid"))

import ellipsoid_alpha_tensor_3d as M3  # noqa: E402


@pytest.fixture(scope="module")
def sphere_alpha():
    # graded mesh: fine air shell around the 5 mm sphere (the magnitude lever)
    mesh = M3.build_mesh((M3.A_SPH,) * 3, box=0.030, maxh=0.012,
                         maxh_cond=8e-4, maxh_shell=1.0e-3)
    w = 2 * math.pi * 300.0
    azz = M3.alpha_tensor_component(mesh, w, "z")
    axx = M3.alpha_tensor_component(mesh, w, "x")
    an = M3.sphere_analytic(w)
    return azz, axx, an


def test_transverse_matches_analytic_sphere(sphere_alpha):
    azz, axx, an = sphere_alpha
    # the TRANSVERSE (x) response is the new capability the axisym solve can't reach
    assert abs(axx - an) / abs(an) < 0.06
    assert abs(azz - an) / abs(an) < 0.06


def test_sphere_is_isotropic(sphere_alpha):
    azz, axx, _ = sphere_alpha
    assert abs(azz - axx) / abs(azz) < 0.02


def test_real_part_is_negative_lift(sphere_alpha):
    azz, axx, _ = sphere_alpha
    # Re[alpha] < 0 (diamagnetic exclusion) in both directions
    assert azz.real < 0.0 and axx.real < 0.0
