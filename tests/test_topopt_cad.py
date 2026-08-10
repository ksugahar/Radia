"""radia.topopt_cad: the design-side half of the shape-regeneration bridge.

Pure-Python contracts (no Cubit): nodal averaging identities, the Exodus
level-set writer round-trip, and the grid iso-surface route on a known
body.  The Cubit-coupled halves (`create tri iso`, `cubit_stl_to_vol`)
are exercised by validation_test/isochronous_topopt/test_shape_regen_lane.py.
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")

from netgen.occ import OCCGeometry, Pnt, Sphere
from ngsolve import Integrate, Mesh, TaskManager

from radia.topopt_cad import (
    iso_stl_from_grid,
    nodal_from_element_density,
    write_levelset_exodus,
    write_vfrac_exodus,
)


@pytest.fixture(scope="module")
def ball():
    with TaskManager():
        return Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0))
                    .GenerateMesh(maxh=0.35))


def test_nodal_average_is_exact_for_constant_density(ball):
    nodal = nodal_from_element_density(ball, np.full(ball.ne, 0.37))
    assert nodal.shape == (ball.nv,)
    np.testing.assert_allclose(nodal, 0.37, rtol=0, atol=1e-14)


def test_nodal_average_is_a_convex_combination(ball):
    rng = np.random.default_rng(7)
    rho = rng.uniform(0.0, 1.0, ball.ne)
    nodal = nodal_from_element_density(ball, rho)
    assert nodal.min() >= rho.min() - 1e-14
    assert nodal.max() <= rho.max() + 1e-14


def test_nodal_average_rejects_bad_input(ball):
    with pytest.raises(ValueError, match="entries"):
        nodal_from_element_density(ball, np.ones(3))
    bad = np.ones(ball.ne)
    bad[0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        nodal_from_element_density(ball, bad)


def test_levelset_exodus_roundtrip(ball, tmp_path):
    pytest.importorskip("netCDF4")
    from netCDF4 import Dataset

    nodal = nodal_from_element_density(ball, np.full(ball.ne, 0.8))
    out = tmp_path / "nested" / "lsd.exo"
    info = write_levelset_exodus(ball, nodal, out, level=0.5)
    assert info["n_tets"] == ball.ne and info["n_nodes"] == ball.nv

    ds = Dataset(str(out), "r")
    try:
        conn = np.asarray(ds.variables["connect1"][:])
        assert conn.shape == (ball.ne, 4)
        assert conn.min() >= 1 and conn.max() <= ball.nv
        raw = ds.variables["name_nod_var"][0, :]
        if hasattr(raw, "filled"):
            raw = raw.filled(b"\x00")
        name = np.asarray(raw).tobytes().replace(b"\x00", b"").decode().strip()
        assert name == "LSD"
        vals = np.asarray(ds.variables["vals_nod_var1"][0, :])
        # level shift: stored field = nodal - level
        np.testing.assert_allclose(vals, 0.8 - 0.5, atol=1e-14)
        # coordinates must be the mesh vertices, same order
        x = np.asarray(ds.variables["coordx"][:])
        assert x.shape == (ball.nv,)
    finally:
        ds.close()


def test_levelset_exodus_rejects_nonfinite_and_invalid_name(ball, tmp_path):
    nodal = np.ones(ball.nv)
    bad = nodal.copy()
    bad[0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        write_levelset_exodus(ball, bad, tmp_path / "bad.exo")
    with pytest.raises(ValueError, match="level"):
        write_levelset_exodus(
            ball, nodal, tmp_path / "bad.exo", level=float("inf"))
    with pytest.raises(ValueError, match="varname"):
        write_levelset_exodus(ball, nodal, tmp_path / "bad.exo",
                              varname="x" * 33)


def test_grid_iso_stl_recovers_a_known_sphere(tmp_path):
    pytest.importorskip("skimage")
    pytest.importorskip("trimesh")
    pytest.importorskip("scipy")
    # The coarse module fixture cannot represent an r=0.6 body in its nodal
    # field.  A finer mesh is required; this is a property of P0->P1
    # averaging, not a bug.
    import ngsolve as ngs
    with TaskManager():
        fine = Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0))
                    .GenerateMesh(maxh=0.15))
    vols = np.asarray(Integrate(ngs.CoefficientFunction(1.0), fine, ngs.VOL,
                                element_wise=True), float)
    cx = np.asarray(Integrate(ngs.x, fine, ngs.VOL, element_wise=True),
                    float) / vols
    cy = np.asarray(Integrate(ngs.y, fine, ngs.VOL, element_wise=True),
                    float) / vols
    cz = np.asarray(Integrate(ngs.z, fine, ngs.VOL, element_wise=True),
                    float) / vols
    r = np.sqrt(cx ** 2 + cy ** 2 + cz ** 2)
    rho = (r <= 0.6).astype(float)
    nodal = nodal_from_element_density(fine, rho)

    out = tmp_path / "sphere.stl"
    info = iso_stl_from_grid(fine, nodal, out, level=0.5, resolution=64)
    assert info["watertight"] is True
    v_exact = 4.0 / 3.0 * np.pi * 0.6 ** 3
    # Allow the expected one-layer inward bias of averaging a P0 step to P1.
    assert abs(info["volume"] - v_exact) / v_exact < 0.25, info["volume"]


def test_grid_iso_stl_rejects_flat_field(ball, tmp_path):
    pytest.importorskip("skimage")
    pytest.importorskip("trimesh")
    pytest.importorskip("scipy")
    nodal = nodal_from_element_density(ball, np.full(ball.ne, 0.9))
    # everything iron at level 0.99 -> nothing to extract must RAISE,
    # not return an empty surface... level above the range:
    with pytest.raises(ValueError, match="outside the sampled"):
        iso_stl_from_grid(ball, nodal, tmp_path / "x.stl", level=0.99)


def test_smoothing_iteration_cap(ball, tmp_path):
    pytest.importorskip("skimage")
    pytest.importorskip("trimesh")
    pytest.importorskip("scipy")
    nodal = nodal_from_element_density(ball, np.ones(ball.ne))
    with pytest.raises(ValueError, match="smooth_iterations"):
        iso_stl_from_grid(ball, nodal, tmp_path / "x.stl",
                          smooth_iterations=10)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"resolution": 15}, "resolution"),
        ({"resolution": 129}, "resolution"),
        ({"resolution": 64.0}, "integer"),
        ({"cutoff_factor": 0.0}, "cutoff_factor"),
        ({"level": float("nan")}, "level"),
        ({"target_faces": 50}, "target_faces"),
    ],
)
def test_grid_iso_stl_rejects_unsafe_arguments(ball, tmp_path, kwargs,
                                                message):
    nodal = np.ones(ball.nv)
    with pytest.raises(ValueError, match=message):
        iso_stl_from_grid(ball, nodal, tmp_path / "x.stl", **kwargs)


def test_grid_iso_stl_rejects_nonfinite_nodal_field(ball, tmp_path):
    nodal = np.ones(ball.nv)
    nodal[0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        iso_stl_from_grid(ball, nodal, tmp_path / "x.stl")


def test_vfrac_exodus_contract_and_sphere_volume(tmp_path):
    pytest.importorskip("netCDF4")
    pytest.importorskip("scipy")
    from netCDF4 import Dataset, chartostring
    import ngsolve as ngs

    # Same body definition as test_grid_iso_stl_recovers_a_known_sphere:
    # the coarse module fixture cannot represent an r=0.6 body in its
    # nodal field (mostly surface vertices), so use the finer mesh --
    # the two regeneration routes then share one A/B-comparable body.
    with TaskManager():
        fine = Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0))
                    .GenerateMesh(maxh=0.15))
    vols = np.asarray(Integrate(ngs.CoefficientFunction(1.0), fine, ngs.VOL,
                                element_wise=True), float)
    cx = np.asarray(Integrate(ngs.x, fine, ngs.VOL, element_wise=True),
                    float) / vols
    cy = np.asarray(Integrate(ngs.y, fine, ngs.VOL, element_wise=True),
                    float) / vols
    cz = np.asarray(Integrate(ngs.z, fine, ngs.VOL, element_wise=True),
                    float) / vols
    rho = (np.sqrt(cx ** 2 + cy ** 2 + cz ** 2) <= 0.6).astype(float)
    nodal = nodal_from_element_density(fine, rho)

    info = write_vfrac_exodus(fine, nodal, tmp_path / "vf", level=0.5,
                              cells=32, supersample=3)
    out = tmp_path / "vf.e.1.0"
    assert info["path"] == str(out) and out.is_file()
    v_exact = 4.0 / 3.0 * np.pi * 0.6 ** 3
    # One-layer inward bias of averaging the P0 step to P1, same band as
    # the STL sibling: the two routes regenerate the SAME blurred body.
    assert abs(info["v_vfrac"] / v_exact - 1.0) < 0.25

    ds = Dataset(str(out))
    try:
        names = [str(s).strip() for s in
                 chartostring(ds.variables["name_elem_var"][:])]
        assert names == ["VOID", "MAT_1"]
        void = np.asarray(ds.variables["vals_elem_var1eb1"][0])
        mat = np.asarray(ds.variables["vals_elem_var2eb1"][0])
        np.testing.assert_allclose(void + mat, 1.0, rtol=0, atol=1e-15)
        assert mat.min() >= 0.0 and mat.max() <= 1.0
        gnames = [str(s).strip() for s in
                  chartostring(ds.variables["name_glo_var"][:])]
        gvals = dict(zip(gnames, np.asarray(ds.variables["vals_glo_var"][0])))
        nx, ny, nz = (int(gvals["gxint"]), int(gvals["gyint"]),
                      int(gvals["gzint"]))
        assert [nx, ny, nz] == info["nel"]
        h = [(gvals["xmax"] - gvals["xmin"]) / nx,
             (gvals["ymax"] - gvals["ymin"]) / ny,
             (gvals["zmax"] - gvals["zmin"]) / nz]
        # cubic cells and the file-recomputed integral matches the report
        np.testing.assert_allclose(h, info["cell_size"], rtol=1e-12)
        np.testing.assert_allclose(mat.sum() * np.prod(h), info["v_vfrac"],
                                   rtol=1e-12)
    finally:
        ds.close()


def test_vfrac_exodus_rejects_bad_input(ball, tmp_path):
    pytest.importorskip("netCDF4")
    pytest.importorskip("scipy")
    with pytest.raises(ValueError, match="entries"):
        write_vfrac_exodus(ball, np.ones(3), tmp_path / "bad")
    with pytest.raises(ValueError, match="zero volume"):
        write_vfrac_exodus(ball, np.zeros(ball.nv), tmp_path / "empty",
                           level=0.5)
