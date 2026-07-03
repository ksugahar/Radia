"""Golden lock for per-region (multi-material) 2D planar MMMM (radia.mmmm2d dict materials).

A motor cross-section body has several soft-iron regions (grades / rotor+stator); mu_r and bh_table
may be passed as {region_name: value} dicts.  Mesh = an inner disk ("inner") + an annulus ("outer").
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.geom2d import SplineGeometry

import radia.mmmm2d as m2


def _two_region_mesh(maxh=0.12):
    geo = SplineGeometry()
    geo.AddCircle((0.0, 0.0), r=2.0, leftdomain=1, rightdomain=0, bc="outer")
    geo.AddCircle((0.0, 0.0), r=1.0, leftdomain=2, rightdomain=1)
    geo.SetMaterial(1, "outer")     # annulus 1 < r < 2
    geo.SetMaterial(2, "inner")     # disk r < 1
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def _region_absM(mesh, M):
    """Mean |M| over each named region."""
    mats = m2._element_materials(mesh)
    out = {}
    for name in set(mats):
        ids = [i for i, m in enumerate(mats) if m == name]
        out[name] = float(np.mean(np.linalg.norm(M[ids], axis=1)))
    return out


def test_uniform_dict_matches_scalar():
    """A per-region dict with the SAME mu_r everywhere == the scalar-mu_r solve."""
    with ng.TaskManager():
        mesh = _two_region_mesh()
        rs = m2.solve_planar_demag(mesh, mu_r=5.0, H_ext=(1000.0, 0.0))
        rd = m2.solve_planar_demag(mesh, mu_r={"inner": 5.0, "outer": 5.0}, H_ext=(1000.0, 0.0))
    assert rd["per_region"] and not rs["per_region"]
    assert np.allclose(rd["M"], rs["M"], rtol=1e-10, atol=1e-8)


def test_uniform_dict_torque_sweep_matches_scalar():
    """The factor-once torque sweep must accept per-region mu_r dicts, too."""
    angles = np.deg2rad([15.0, 40.0, 80.0])
    with ng.TaskManager():
        mesh = _two_region_mesh()
        rs = m2.torque_angle_sweep(mesh, 1000.0, angles, Rc=3.0, mu_r=5.0, n=360)
        rd = m2.torque_angle_sweep(
            mesh, 1000.0, angles, Rc=3.0,
            mu_r={"inner": 5.0, "outer": 5.0}, n=360,
        )
    assert rs["factored_once"] and rd["factored_once"]
    assert rd["per_region"] and not rs["per_region"]
    assert np.allclose(rd["M_avg"], rs["M_avg"], rtol=1e-10, atol=1e-8)
    assert np.allclose(rd["torque"], rs["torque"], rtol=1e-10, atol=1e-12)


def test_heterogeneous_regions():
    """A high-permeability inner disk magnetizes much more strongly than a near-air outer ring."""
    with ng.TaskManager():
        mesh = _two_region_mesh()
        r = m2.solve_planar_demag(mesh, mu_r={"inner": 200.0, "outer": 1.5}, H_ext=(1000.0, 0.0))
    absM = _region_absM(mesh, r["M"])
    assert absM["inner"] > 5.0 * absM["outer"], absM
    assert r["demag_factors"] is None       # not a single-body concept


def test_missing_region_raises():
    """A mesh region absent from the dict raises with the available regions (no silent fallback)."""
    with ng.TaskManager():
        mesh = _two_region_mesh()
        with pytest.raises(ValueError, match="inner|outer|region"):
            m2.solve_planar_demag(mesh, mu_r={"inner": 100.0}, H_ext=(1000.0, 0.0))


def test_nonlinear_per_region():
    """Per-region bh_table dict: two soft-iron grades converge; the softer grade magnetizes more."""
    soft = np.array([[0, 0.0], [100, 1.2], [1000, 1.6], [1e5, 2.0], [1e6, 2.2]])
    hard = np.array([[0, 0.0], [500, 0.6], [5000, 1.0], [1e5, 1.5], [1e6, 1.8]])
    with ng.TaskManager():
        mesh = _two_region_mesh()
        r = m2.solve_planar_demag(mesh, bh_table={"inner": soft, "outer": hard},
                                  H_ext=(3000.0, 0.0), nl_tol=1e-4)
    assert r["nonlinear"] and r["per_region"] and r["iters"] >= 1
    absM = _region_absM(mesh, r["M"])
    assert absM["inner"] > absM["outer"], absM
