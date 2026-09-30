"""All three Kelvin vector direct routes use the shared gate and agree with CG.

This is a small algebraic manufactured-source comparison, not coil accuracy.
The periodic geometry is shared; HCurl is mass-gauged, not point-grounded.
"""
import math

import numpy as np
import pytest
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import cg
import ngsolve as ng
from netgen.occ import Sphere, Pnt, Glue, OCCGeometry, IdentificationType

from radia import kelvin_solver as solver
from radia import _residual_gate as gate


@pytest.fixture(scope='module')
def mesh():
    inner = Sphere(Pnt(0, 0, 0), 1)
    inner.mat('iron')
    inner.faces.name = 'kelvin_int'
    outer = Sphere(Pnt(3, 0, 0), 1)
    outer.mat('kelvin')
    outer.faces.name = 'kelvin_ext'
    inner.faces[0].Identify(outer.faces[0], 'kelvin', IdentificationType.PERIODIC)
    with ng.TaskManager():
        return ng.Mesh(OCCGeometry(Glue([inner, outer])).GenerateMesh(maxh=.8))


def _run(mesh, route):
    options = dict(R_K=1., offset=(3., 0., 0.), nu_0=1., gauge_eps=1e-3)
    if route == 'full':
        return solver.solve_full_A_kelvin(mesh, ng.CF((-ng.y, ng.x, 0)),
                                          source_material='iron', **options)
    source = solver.project_kelvin_A_source(mesh, ng.CF((-ng.y, ng.x, 0)), order=1)
    if route == 'reduced':
        return solver.solve_reduced_A_kelvin(mesh, source, **options)
    return solver.solve_magnetostatic_reduced_A_kelvin(
        mesh, source, mu_r_by_material={'iron': 5.}, **options)


@pytest.mark.parametrize('route', ['full', 'reduced', 'magnetostatic_reduced'])
def test_kelvin_direct_matches_independent_cg_and_calls_shared_gate(mesh, route, monkeypatch):
    captured = []
    check = gate.check_true_residual

    def observe(matrix, residual, solution, rhs, free, what, **kwargs):
        value = check(matrix, residual, solution, rhs, free, what, **kwargs)
        values, columns, pointers = matrix.CSR()
        matrix_copy = csr_matrix((np.asarray(values).copy(), np.asarray(columns).copy(),
                                  np.asarray(pointers).copy()), shape=(len(rhs), len(rhs)))
        captured.append((matrix_copy, rhs.copy(), free.copy(), value))
        return value

    monkeypatch.setattr(gate, 'check_true_residual', observe)
    with ng.TaskManager():
        result = _run(mesh, route)
    assert len(captured) == 1
    matrix, rhs, free, residual = captured[0]
    assert residual <= gate.RELATIVE_LIMIT
    assert np.linalg.norm(rhs[free]) > 0
    reduced = matrix[free][:, free]
    reference, info = cg(reduced, rhs[free], M=diags(1./reduced.diagonal()),
                         rtol=1e-11, atol=0., maxiter=4000)
    assert info == 0
    ref = ng.GridFunction(result['fes'])
    ref.vec.FV().NumPy()[free] = reference
    direct = result['gfu' if route == 'full' else 'gfu_r']
    with ng.TaskManager():
        b = ng.curl(direct)
        delta = b-ng.curl(ref)
        field_error = ng.Integrate(ng.InnerProduct(delta, delta), mesh)
        field_scale = ng.Integrate(ng.InnerProduct(b, b), mesh)
    assert field_scale > 0
    assert math.sqrt(abs(field_error)/field_scale) < 1e-8
    actual = direct.vec.FV().NumPy()[free]
    energy = float(actual @ (reduced @ actual))
    reference_energy = float(reference @ (reduced @ reference))
    assert abs(energy/reference_energy-1) < 1e-8


@pytest.mark.parametrize('route', ['full', 'reduced', 'magnetostatic_reduced'])
def test_kelvin_direct_rejects_residual_over_shared_limit(mesh, route, monkeypatch):
    check = gate.check_true_residual

    def inject(matrix, residual, solution, rhs, free, what, **kwargs):
        bad = np.zeros_like(residual)
        bad[np.flatnonzero(free)[0]] = 2e-6*np.linalg.norm(rhs[free])
        return check(matrix, bad, solution, rhs, free, what, **kwargs)

    monkeypatch.setattr(gate, 'check_true_residual', inject)
    with ng.TaskManager(), pytest.raises(RuntimeError, match='true relative residual'):
        _run(mesh, route)
