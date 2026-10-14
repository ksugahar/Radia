"""Polynomial reproduction, patch failure, and mixed-source preservation."""

import itertools

import numpy as np
import pytest

from radia.patch_recovery import recover_vertex_patches, recover_mixed_omega_p1


def cube_tets(n=3):
    points = np.array(list(itertools.product(range(n + 1), repeat=3)), float) / n
    indices = np.arange((n + 1) ** 3).reshape((n + 1,) * 3)
    cells = []
    for origin in itertools.product(range(n), repeat=3):
        for permutation in itertools.permutations(range(3)):
            node = np.array(origin)
            cell = [indices[tuple(node)]]
            for axis in permutation:
                node = node.copy()
                node[axis] += 1
                cell.append(indices[tuple(node)])
            cells.append(cell)
    return points, np.array(cells)


def test_linear_patch_reproduction_with_translation_and_scaling():
    points, cells = cube_tets()
    for length, origin in ((1., 0.), (1.e-5, 1.), (1.e4, -2.e3)):
        coords = origin + length * points
        centers = coords[cells].mean(axis=1)
        matrix = np.array([[1., 2., 0.], [0., -1., .5], [3., 0., 1.]])
        samples = ((centers - origin) / length) @ matrix.T + [1., 2., 3.]
        patch = recover_vertex_patches(coords, cells, samples)
        expected = points[patch['vertex_ids']] @ matrix.T + [1., 2., 3.]
        np.testing.assert_allclose(patch['values'], expected, rtol=1e-9, atol=1e-9)
        assert patch['expanded_patches'] > 0  # corners need a larger patch


def test_separate_material_patches_keep_an_interface_jump():
    points, cells = cube_tets(n=4)
    centers = points[cells].mean(axis=1)
    recovered = []
    for mask, shift in ((centers[:, 0] < .5, 0.), (centers[:, 0] > .5, 10.)):
        patch = recover_vertex_patches(points, cells[mask], centers[mask] + shift)
        ids = patch['vertex_ids']
        interface = points[ids, 0] == .5
        recovered.append(patch['values'][interface])
    np.testing.assert_allclose(recovered[1] - recovered[0], 10., atol=1e-12)


def test_rank_deficient_patch_does_not_silently_average():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    with pytest.raises(ValueError, match='rank-deficient'):
        recover_vertex_patches(points, np.array([[0, 1, 2, 3]]), np.ones((1, 3)))


def test_mixed_recovery_preserves_analytic_source_and_raw_fields():
    import ngsolve as ng
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    left = Box(Pnt(0, 0, 0), Pnt(.5, 1, 1)); left.mat('reduced')
    right = Box(Pnt(.5, 0, 0), Pnt(1, 1, 1)); right.mat('total')
    mesh = ng.Mesh(OCCGeometry(Glue([left, right])).GenerateMesh(maxh=.22))
    space = ng.H1(mesh, order=1)
    source = ng.CF((ng.y**2, 0, 0))
    gradient_left = ng.CF((1, 2, 3))
    gradient_right = ng.CF((4, 5, 6))
    raw = mesh.MaterialCF({'reduced': source-gradient_left, 'total': ng.CF((0, ng.z**2, 0))-gradient_right})
    mu = mesh.MaterialCF({'reduced': 1., 'total': 5.})
    left_phi = ng.GridFunction(space)
    right_phi = ng.GridFunction(space)
    with ng.TaskManager():
        left_phi.Set(ng.x + 2*ng.y + 3*ng.z)
        right_phi.Set(4*ng.x + 5*ng.y + 6*ng.z)
    result = {'fes': space, 'H_cf': raw, 'B_cf': mu*raw, 'mu_cf': mu,
              'phi_reduced': left_phi, 'phi_total': right_phi}
    recovered = recover_mixed_omega_p1(mesh, result, reduced_materials=('reduced',),
                                      total_materials=('total',))
    assert recovered['raw_H_cf'] is raw
    assert result['H_cf'] is raw
    for coords in ((.2, .123, .27), (.8, .31, .24)):
        point = mesh(*coords)
        np.testing.assert_allclose(recovered['H_cf'](point), raw(point), atol=1e-10)
        np.testing.assert_allclose(recovered['B_cf'](point), result['B_cf'](point), atol=1e-10)
    assert not recovered['diagnostics']['accuracy_accepted']
    with pytest.raises(ValueError, match='P1'):
        recover_mixed_omega_p1(mesh, dict(result, fes=ng.H1(mesh, order=2)),
                              reduced_materials=('reduced',), total_materials=('total',))

@pytest.mark.parametrize('rings', [0, -1, 1.5, True, float('nan')])
def test_patch_ring_count_is_a_positive_integer(rings):
    points, cells = cube_tets()
    with pytest.raises(ValueError, match='positive integer'):
        recover_vertex_patches(points, cells, np.ones((len(cells), 3)), max_patch_rings=rings)


def test_degenerate_tetrahedra_are_rejected():
    with pytest.raises(ValueError, match='degenerate'):
        recover_vertex_patches(np.zeros((4, 3)), np.array([[0, 1, 2, 3]]), np.ones((1, 3)))


def test_condensed_solver_recovery_keeps_solution_and_source():
    import ngsolve as ng
    from netgen.occ import Box, Glue, OCCGeometry, Pnt, X
    from radia.kelvin_solver import solve_magnetostatic_matching_trace_total_reduced_omega

    left = Box(Pnt(-1, -1, -1), Pnt(0, 1, 1)); left.mat('reduced')
    right = Box(Pnt(0, -1, -1), Pnt(1, 1, 1)); right.mat('total')
    for solid in (left, right):
        solid.faces.name = 'natural'
    left.faces.Min(X).name = 'fixed'; right.faces.Max(X).name = 'fixed'
    left.faces.Max(X).name = 'interface'; right.faces.Min(X).name = 'interface'
    mesh = ng.Mesh(OCCGeometry(Glue([left, right])).GenerateMesh(maxh=.4))
    trace = ng.GridFunction(ng.H1(mesh, order=1, definedon=mesh.Boundaries('interface')))
    with ng.TaskManager():
        result = solve_magnetostatic_matching_trace_total_reduced_omega(
            mesh, ng.CF((1, 0, 0)), trace, reduced_materials=('reduced',),
            total_materials=('total',), interface_boundary='interface',
            dirichlet_boundary='fixed', order=1, inverse='sparsecholesky',
            mu_r_by_material={'reduced': 1., 'total': 2.})
        before = result['solution'].vec.FV().NumPy().copy()
        recovered = recover_mixed_omega_p1(mesh, result, reduced_materials=('reduced',),
                                          total_materials=('total',))
    np.testing.assert_array_equal(result['solution'].vec.FV().NumPy(), before)
    for point in [mesh(-.3, .13, .21), mesh(.4, -.2, .1)]:
        np.testing.assert_allclose(recovered['H_cf'](point), result['H_cf'](point), atol=1e-9)
    for override, message in [({'nonlinear_stats': {}}, 'linear'),
                              ({'kelvin_materials': ('total',)}, 'linear')]:
        with pytest.raises(ValueError, match=message):
            recover_mixed_omega_p1(mesh, dict(result, **override),
                                  reduced_materials=('reduced',), total_materials=('total',))
