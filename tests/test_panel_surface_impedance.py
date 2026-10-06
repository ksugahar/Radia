"""Independent flat-triangle forms, source placement and conservative heat."""
import numpy as np
import pytest
import ngsolve as ng
import netgen.meshing as nm

from radia.surface_impedance import (
    PanelSurfaceImpedance, assemble_panel_impedance_stiffness,
    panel_tangential_field_rms,
)
from radia.bem_sibc_solver import ScalarBIESIBCSolver, _project_sibc_surface_heat, MU_0
from radia.workpiece_surface import _complete_sibc_reaction


def tetra_surface():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    triangles = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
    mesh = nm.Mesh(dim=3)
    face = mesh.Add(nm.FaceDescriptor(surfnr=1, domin=1, bc=1))
    vertices = [mesh.Add(nm.MeshPoint(nm.Pnt(*p))) for p in points]
    for tri in triangles:
        mesh.Add(nm.Element2D(face, [vertices[i] for i in tri]))
    return ng.Mesh(mesh), points, triangles


def independent_form(points, triangles, z):
    matrix = np.zeros((len(points), len(points)), dtype=complex)
    areas, gradients = [], []
    for tri, value in zip(triangles, z):
        p = points[tri]
        normal = np.cross(p[1]-p[0], p[2]-p[0])
        twice_area = np.linalg.norm(normal)
        normal /= twice_area
        gradients.append(np.array([np.cross(normal, p[(k+2)%3]-p[(k+1)%3])
                                   / twice_area for k in range(3)]))
        areas.append(twice_area/2)
        matrix[np.ix_(tri, tri)] += value * twice_area/2 * (gradients[-1] @ gradients[-1].T)
    return matrix, np.array(areas), np.array(gradients)


def test_weighted_form_and_local_heat_match_independent_panel_integrals():
    mesh, points, triangles = tetra_surface()
    z = np.array([1+2j, 3+.5j, .1+4j, 2+1j])
    phi = np.array([1+2j, -2+.3j, .5-1j, 3+0j])
    with ng.TaskManager():
        fes = ng.H1(mesh, order=1)
        panel = PanelSurfaceImpedance(z)
        actual = assemble_panel_impedance_stiffness(fes, panel)
        expected, areas, gradients = independent_form(points, triangles, z)
        np.testing.assert_allclose(actual, expected, atol=2e-14)
        np.testing.assert_allclose(actual, actual.T, atol=1e-14)
        np.testing.assert_allclose(actual @ np.ones(4), 0, atol=2e-14)
        h = np.einsum('tik,ti->tk', gradients, phi[triangles])
        q = .5*z.real*np.sum(abs(h)**2, axis=1)
        heat, power = _project_sibc_surface_heat(fes, phi, panel)
        assert power == pytest.approx(areas @ q, rel=2e-14)
        assert np.min(heat.vec.FV().NumPy()) >= 0
        np.testing.assert_allclose(panel_tangential_field_rms(fes, phi),
                                   np.linalg.norm(h, axis=1), rtol=2e-14)
        source = np.array([.1+2j, 1., -3., 2j])
        reaction = _complete_sibc_reaction(1j, source, phi, None, panel, 7., 2.,
                                          weighted_stiffness=actual)
        assert reaction == pytest.approx(1j + source @ expected @ phi / (1j*7*4))


@pytest.mark.parametrize('bad', [[], [[1]], [1, -1], [1, np.nan], [1, np.inf]])
def test_bad_panel_impedance_rejected(bad):
    with pytest.raises(ValueError, match='finite passive'):
        PanelSurfaceImpedance(bad)


def test_panel_layout_is_not_a_nodal_array():
    mesh, _, _ = tetra_surface()
    with ng.TaskManager():
        with pytest.raises(ValueError, match='one value per BND'):
            assemble_panel_impedance_stiffness(ng.H1(mesh, order=1), PanelSurfaceImpedance([1, 2]))


def test_dense_and_hacapk_source_weighting_and_constant_reduction():
    mesh, _, _ = tetra_surface()
    with ng.TaskManager():
        solver = ScalarBIESIBCSolver(mesh, order=1, use_intree_bem=True,
                                    use_intree_hacapk=True, intree_geom_order=1,
                                    intree_singular_n_q=6, intree_regular_quad_degree=7)
        source = np.array([1+1j, 2-1j, -1., .3j])
        z = PanelSurfaceImpedance([.001+.002j, .003+.004j, .002+.001j, .005+.001j])
        omega = 1000.
        kz = solver.impedance_stiffness(z)
        a = .5*solver.M-solver.DL + solver.SL @ solver.M_inv @ kz/(1j*omega*MU_0)
        aug = np.block([[a, solver._c_gauge[:, None]],
                        [solver._c_gauge[None, :], np.zeros((1, 1))]])
        expected = np.linalg.solve(aug, np.r_[solver.M @ source, 0])[:solver.ndof]
        dense = solver.solve(source, z, omega)
        hacapk = solver.solve_hacapk(source, z, omega, tol=1e-11)
        np.testing.assert_allclose(dense['phi_vec'], expected, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(hacapk['phi_vec'], dense['phi_vec'], rtol=1e-9, atol=1e-9)
        assert hacapk['P_density'] == pytest.approx(dense['P_density'], rel=1e-9)
        scalar = solver.solve(source, .002+.003j, omega)
        constant = solver.solve(source, PanelSurfaceImpedance(np.full(4, .002+.003j)), omega)
        np.testing.assert_allclose(scalar['phi_vec'], constant['phi_vec'], rtol=1e-12, atol=1e-12)
        assert scalar['P_density'] == pytest.approx(constant['P_density'], rel=1e-12)
