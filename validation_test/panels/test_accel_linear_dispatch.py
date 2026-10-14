"""Exercise accelerator solver selection on unchanged H1 and full HCurl spaces."""
import numpy as np
import pytest

from radia.panels.calc_accel_magnet import _accel_inverse
from radia.panels.calc_common import apply_fe_inverse


@pytest.mark.parametrize("form,order,backend", [
    ("omega", 2, "sparsecholesky"), ("omega", 2, "bddc"),
    ("omega", 2, "iccg"), ("a", 1, "iccg"),
    ("a", 1, "bddc_ams"), ("a", 2, "bddc_ams"), ("a", 3, "bddc_ams"),
])
def test_accel_backend_matches_direct(form, order, backend):
    from netgen.csg import unit_cube
    from ngsolve import Mesh, H1, HCurl, BilinearForm, LinearForm, GridFunction, TaskManager, CF, dx, grad, curl
    mesh = Mesh(unit_cube.GenerateMesh(maxh=.45))
    fes = (H1 if form == "omega" else HCurl)(mesh, order=order, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    if form == "omega":
        a += (grad(u) * grad(v) + u * v) * dx
    else:
        a += (curl(u) * curl(v) + 1e-6 * u * v) * dx
    f = LinearForm(fes)
    f += (v if form == "omega" else CF((1, 2, 3)) * v) * dx
    with TaskManager():
        a.Assemble()
        f.Assemble()
        reference = GridFunction(fes)
        reference.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    expected = reference.vec.FV().NumPy().copy()
    inverse = _accel_inverse(a, fes, backend)
    result = GridFunction(fes)
    with TaskManager():
        residual = apply_fe_inverse(a.mat, inverse, f.vec, result.vec, fes.FreeDofs())
    assert residual < 1e-8
    assert np.linalg.norm(result.vec.FV().NumPy() - expected) / np.linalg.norm(expected) < 1e-5


@pytest.mark.parametrize("order", [1, 2, 3])
def test_accel_application_energy_matches_direct(tmp_path, order):
    from netgen.occ import Box, Pnt, OCCGeometry
    from radia.em_material import EMMaterial
    from radia.panels.calc_accel_magnet import solve_accel

    body = Box(Pnt(-.01, -.01, -.01), Pnt(.01, .01, .01))
    body.mat("yoke")
    body.faces.name = "outer"
    mesh = OCCGeometry(body).GenerateMesh(maxh=.012)
    vol = tmp_path / "yoke.vol"
    mesh.Save(str(vol))
    coil = tmp_path / "coil.py"
    coil.write_text(
        "from radia.coil_builder import CoilBuilder\n"
        "def build_coil():\n"
        "    return (CoilBuilder(current=100).set_start([.03,-.05,0])\n"
        "            .set_cross_section(.002,.002).add_straight(.1))\n",
        encoding="utf-8",
    )
    kwargs = dict(coil_script=str(coil), vol_file=str(vol), formulation="a",
                  fes_order=order, mat=EMMaterial("linear", 0, 100))
    direct = solve_accel(**kwargs, solver="sparsecholesky")
    iterative = solve_accel(**kwargs, solver="auto")
    assert iterative["linear_solver"] == "bddc_ams"
    assert direct["ndof"] == iterative["ndof"]
    for result in (direct, iterative):
        assert result["converged"]
        assert result["linear_relative_residual"] < 1e-8
        assert result["W_mag"] > 0
    assert iterative["W_mag"] == pytest.approx(direct["W_mag"], rel=1e-6)
    np.testing.assert_allclose(iterative["B_origin"], direct["B_origin"], rtol=1e-6, atol=1e-12)


def test_accel_cadless_flat_mesh_keeps_stored_order(tmp_path):
    """A CAD-less order-1 .vol without a Kelvin shell runs at FE order 2.

    FE order and geometry order are independent; the panel must neither
    re-curve the loaded mesh nor reject it when no curved boundary needs it.
    """
    from netgen.occ import Box, Pnt, OCCGeometry
    from radia.em_material import EMMaterial
    from radia.panels.calc_accel_magnet import solve_accel

    body = Box(Pnt(-.01, -.01, -.01), Pnt(.01, .01, .01))
    body.mat("yoke")
    body.faces.name = "outer"
    with_cad = tmp_path / "yoke.vol"
    OCCGeometry(body).GenerateMesh(maxh=.012).Save(str(with_cad))
    lines = with_cad.read_text().splitlines()
    cadless = tmp_path / "yoke_cadless.vol"
    cadless.write_text("\n".join(lines[:lines.index("endmesh") + 1]) + "\n")
    coil = tmp_path / "coil.py"
    coil.write_text(
        "from radia.coil_builder import CoilBuilder\n"
        "def build_coil():\n"
        "    return (CoilBuilder(current=100).set_start([.03,-.05,0])\n"
        "            .set_cross_section(.002,.002).add_straight(.1))\n",
        encoding="utf-8",
    )
    common = dict(coil_script=str(coil), formulation="a", fes_order=2,
                  mat=EMMaterial("linear", 0, 100), solver="sparsecholesky")
    reference = solve_accel(vol_file=str(with_cad), **common)
    result = solve_accel(vol_file=str(cadless), **common)
    assert "error" not in result, result.get("error")
    assert result["converged"]
    assert result["W_mag"] == pytest.approx(reference["W_mag"], rel=1e-9)
