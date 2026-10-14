"""Public semantic validation for the origin-anchored nonlinear BH route."""

import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve.meshes import MakeStructured3DMesh

from radia import vim


def _mesh():
    return MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1)


@pytest.mark.parametrize(
    "table, message",
    [
        ([[-1.0, 0.0], [0.0, 1.0]], "H=0"),
        ([[0.0, 0.1], [1.0, 1.0]], "B=0"),
        ([[0.0, 0.0], [1.0, 1.0], [1.0, 1.1]], "strictly increasing"),
        ([[0.0, 0.0], [1.0, 1.0], [2.0, 0.9]], "non-decreasing"),
    ],
)
def test_soft_iron_bh_table_rejects_noncanonical_curves(table, message):
    with pytest.raises(ValueError, match=message):
        vim.Solve(
            _mesh(),
            bh_table=table,
            H_ext=ng.CoefficientFunction((0.0, 0.0, 0.0)),
        )


def test_image_plane_crossing_mesh_is_rejected():
    """A reduced mesh must lie on one side of every image plane."""
    mesh = MakeStructured3DMesh(
        hexes=True, nx=2, ny=1, nz=1, mapping=lambda x, y, z: (2.0 * x - 1.0, y, z))
    with ng.TaskManager():
        with pytest.raises(ValueError, match="crosses the image plane x=0"):
            vim.Solve(mesh, mu_r=100.0, H_ext=ng.CoefficientFunction((0.0, 0.0, 1.0)),
                      image="-x")


def test_mesh_soft_iron_rejects_a_filter_it_cannot_solve():
    """rad.Solve solves the whole registered mesh; a strict subset must fail early."""
    mesh = MakeStructured3DMesh(hexes=True, nx=2, ny=1, nz=1)
    with pytest.raises(NotImplementedError, match="material_filter selected"):
        vim.MeshSoftIron(mesh, mu_r=100.0, material_filter=["nothing"])


def test_curve_mesh_keeps_a_higher_order_geometry():
    """A mesh already curved above the requested order is not re-curved."""
    from netgen.occ import OCCGeometry, Sphere, Pnt
    from radia.vim._vim import _curve_mesh

    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=0.8))
    mesh.Curve(3)
    _curve_mesh(mesh, 2)
    assert mesh.GetCurveOrder() == 3


def test_falling_magnetization_tail_is_reported_not_silently_capped():
    """dB/dH < mu0 makes M fall; the capped law is reported, not hidden."""
    falling = [[0.0, 0.0], [500.0, 1.3], [5000.0, 1.85], [5e5, 2.4]]
    physical = [[0.0, 0.0], [500.0, 1.3], [5000.0, 1.85], [5e5, 2.5]]
    applied = ng.CoefficientFunction((0.0, 0.0, 1.0e4))
    with ng.TaskManager():
        with pytest.warns(RuntimeWarning, match="falls between H=5000 and H=500000"):
            capped = vim.Solve(_mesh(), bh_table=falling, H_ext=applied)
        clean = vim.Solve(_mesh(), bh_table=physical, H_ext=applied)
    assert capped["bh_table_falling_magnetization"] == {None: [5000.0, 500000.0]}
    assert clean["bh_table_falling_magnetization"] is None
