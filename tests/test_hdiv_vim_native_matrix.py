"""Fast contract for the native NGSolve HDiv demagnetization matrix."""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

from netgen.csg import unit_cube  # noqa: E402
from radia.vim import DemagOperator  # noqa: E402


def test_demag_operator_is_native_ngsolve_matrix_and_matches_configured_apply():
    mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=1.0))
    fes = ng.HDiv(mesh, order=1)
    with ng.TaskManager():
        operator = DemagOperator(fes, eps=1e-10)
        magnetization = ng.GridFunction(fes)
        magnetization.Set(ng.CF((0.2, -0.1, 0.3)))
        result = magnetization.vec.CreateVector()
        result.data = operator.mat * magnetization.vec

        coefficients = np.ascontiguousarray(
            magnetization.vec.FV().NumPy(), dtype=np.float64)
        reference = operator._G.apply_configured_demag(coefficients, True)

        trial, test = fes.TnT()
        mass = ng.BilinearForm(fes)
        mass += trial * test * ng.dx
        mass.Assemble()
        mass_reference = magnetization.vec.CreateVector()
        mass_reference.data = mass.mat * magnetization.vec
        mass_native = operator._G.apply_configured_geometry_mass(coefficients)

        added = magnetization.vec.CreateVector()
        added[:] = 1.25
        operator.mat.MultAdd(0.4, magnetization.vec, added)

        transposed = magnetization.vec.CreateVector()
        transposed[:] = -0.2
        operator.mat.MultTransAdd(-0.3, magnetization.vec, transposed)

    assert isinstance(operator.mat, ng.BaseMatrix)
    assert type(operator.mat).__module__ == "radia._radia_pybind"
    assert operator.mat.height == fes.ndof == operator.mat.width
    assert np.array_equal(result.FV().NumPy(), reference)
    assert np.allclose(mass_native, mass_reference.FV().NumPy(), rtol=2e-15, atol=1e-15)
    assert np.array_equal(added.FV().NumPy(), 1.25 + 0.4 * reference)
    assert np.array_equal(transposed.FV().NumPy(), -0.2 - 0.3 * reference)
