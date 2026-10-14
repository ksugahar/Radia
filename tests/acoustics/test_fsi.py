"""Fast regression for radia.acoustics.fsi (DtN FSI, NGSolve VectorH1 interior).

Coarse-mesh sanity: the coupled solve runs, the spherical DtN gate passes, and the
P2 interior beats P1 on the same mesh (the whole point of the NGSolve order-p
interior).  The full convergence-to-Faran study lives in validation_test/acoustics.
"""
import numpy as np
import pytest
from ngsolve import Mesh, TaskManager

from radia.acoustics import elastic_sphere_scattering, fsi, rigid_sphere_scattering

K, R = 2.0, 1.0
OBS = np.array([[0.0, 0.0, -3.0], [3.0, 0.0, 0.0], [0.0, 0.0, 3.0], [0.0, 0.0, 2.0]])


def _rel(a, b):
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


def test_p2_more_accurate_than_p1_on_same_mesh():
    with TaskManager():
        mesh = fsi.sphere_mesh(R, maxh=0.4)
        far = elastic_sphere_scattering(K, R, OBS, longitudinal_speed=2.0,
                                        shear_speed=1.0, density_ratio=1.5)["scattered"]
        s1 = fsi.fsi_dtn_solve(mesh, K, order=1, obs=OBS)
        s2 = fsi.fsi_dtn_solve(mesh, K, order=2, obs=OBS)

    assert s2["residual"] < 1e-8
    assert s1["dtn"]["deviation"] < 3e-2                     # spherical truncation gate
    assert len(s2["c"]) == (s2["dtn"]["degree"] + 1) ** 2    # (N+1)^2 harmonic coeffs
    assert s2["ndof_u"] > s1["ndof_u"]                       # P2 has more interior DOFs
    # P2 (O(h^2)) is strictly more accurate than P1 (O(h)) on the same mesh
    assert _rel(s2["scattered"], far) < _rel(s1["scattered"], far)


def test_stiff_limit_approaches_rigid_sphere():
    with TaskManager():
        mesh = fsi.sphere_mesh(R, maxh=0.32)
        s = fsi.fsi_dtn_solve(mesh, K, cL=50.0, cT=25.0, rho_s=50.0, order=2, obs=OBS)
        rig = rigid_sphere_scattering(K, R, OBS)["scattered"]
    # a very stiff/heavy elastic sphere radiates like a sound-hard sphere
    assert _rel(s["scattered"], rig) < 0.15


def test_non_spherical_truncation_raises():
    from netgen.occ import Box, OCCGeometry
    with TaskManager():
        box = Box((-1, -1, -1), (1, 1, 1))
        box.faces.name = "gamma"
        mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.7))
        with pytest.raises(ValueError):                      # DtN is fail-loud off-sphere
            fsi.fsi_dtn_solve(mesh, K, order=1, obs=OBS)
