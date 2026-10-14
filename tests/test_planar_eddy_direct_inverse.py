"""Direct-solver contract for the planar reduced-potential eddy operator."""
from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")
MU0 = 4e-7 * np.pi


def _eddy_inverse_name():
    source_path = Path(__file__).parents[1] / "src" / "radia" / "planar_eddy.py"
    tree = ast.parse(source_path.read_text(encoding="utf-8"))
    function = next(
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "eddy_operator")
    inverse_calls = [
        node for node in ast.walk(function)
        if (isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "Inverse")]
    assert len(inverse_calls) == 1
    keyword = next(
        item for item in inverse_calls[0].keywords if item.arg == "inverse")
    return ast.literal_eval(keyword.value)


def test_eddy_operator_uses_sparsecholesky():
    assert _eddy_inverse_name() == "sparsecholesky"


def test_sparsecholesky_solves_the_complex_planar_eddy_system():
    from netgen.geom2d import unit_square

    mesh = ng.Mesh(unit_square.GenerateMesh(maxh=0.25))
    mesh.Curve(2)
    sigma = 3.7e7
    freq = 50.0
    omega = 2.0 * np.pi * freq
    fes = ng.H1(
        mesh, order=2, complex=True,
        dirichlet="left|right|top|bottom")
    u, v = fes.TnT()
    matrix = ng.BilinearForm(fes, symmetric=True)
    matrix += ng.grad(u) * ng.grad(v) * ng.dx
    matrix += 1j * omega * MU0 * sigma * u * v * ng.dx("default")
    rhs = ng.LinearForm(fes)
    rhs += -1j * omega * MU0 * sigma * ng.x * v * ng.dx("default")
    matrix.Assemble()
    rhs.Assemble()
    inverse = matrix.mat.Inverse(fes.FreeDofs(), inverse=_eddy_inverse_name())

    assert type(inverse).__name__.startswith("SparseCholesky")

    solution = rhs.vec.CreateVector()
    solution.data = inverse * rhs.vec
    residual = rhs.vec.CreateVector()
    residual.data = matrix.mat * solution - rhs.vec
    free = np.asarray(list(fes.FreeDofs()), dtype=bool)
    residual_values = residual.FV().NumPy()[free]
    rhs_values = rhs.vec.FV().NumPy()[free]
    relative_residual = np.linalg.norm(residual_values) / np.linalg.norm(rhs_values)
    assert relative_residual < 1e-12


def test_eddy_operator_factor_is_checked_and_rejects_a_singular_system():
    """The reused factor is verified once; a Neumann-only magnetostatic system
    (freq 0, no Dirichlet boundary) is singular and must raise, not solve."""
    from netgen.geom2d import unit_square
    from radia.planar_eddy import eddy_operator

    mesh = ng.Mesh(unit_square.GenerateMesh(maxh=0.25))
    fes, inverse, omega = eddy_operator(mesh, 3.7e7, 50.0, order=2,
                                        conductor="default", dirichlet="left|right|top|bottom")
    assert fes.ndof > 0 and omega == pytest.approx(2 * np.pi * 50.0)
    with pytest.raises(RuntimeError, match="planar eddy SparseCholesky factor"):
        eddy_operator(ng.Mesh(unit_square.GenerateMesh(maxh=0.25)), 3.7e7, 0.0, order=2,
                      conductor="default", dirichlet="")
