"""Discrete source adjoint and native current-average checks."""
import importlib.util
from pathlib import Path

import ngsolve as ng
import numpy as np
import pytest

from radia.bem_complete_reaction import (electric_incident_vertex_load,
    surface_source_reaction_load, surface_current_average_maps)
from radia.bem_sibc_solver import (SurfacePoissonPhiInc,
    H_from_surface_J_complex, A_from_surface_J)
from radia.bem_coupled_solver import extract_element_J


@pytest.fixture(scope="module")
def surfaces():
    path = Path(__file__).resolve().parents[1]/"validation_test/induction_heating/run_loop_work_ring.py"
    spec = importlib.util.spec_from_file_location("reaction_ring", path)
    ring = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ring)
    mesh, points, tris = ring.ring_mesh(16)
    source, _, _ = ring.ring_mesh(16)
    return mesh, points, tris, source


def test_average_maps_match_native_real_and_complex(surfaces):
    _, _, _, mesh = surfaces
    rng = np.random.default_rng(214)
    with ng.TaskManager():
        fes = ng.HDivSurface(mesh, order=1)
        maps = surface_current_average_maps(fes)
        coefficients = rng.normal(size=fes.ndof)+1j*rng.normal(size=fes.ndof)
        real, imag = ng.GridFunction(fes), ng.GridFunction(fes)
        real.vec.FV().NumPy()[:] = coefficients.real
        imag.vec.FV().NumPy()[:] = coefficients.imag
        _, _, direct_real = extract_element_J(mesh, real)
        _, _, direct_imag = extract_element_J(mesh, imag)
    direct = direct_real+1j*direct_imag
    mapped = np.stack([m @ coefficients for m in maps], axis=1)
    np.testing.assert_allclose(mapped, direct, rtol=1e-12, atol=1e-11)
    np.testing.assert_allclose(np.stack([m @ coefficients.real for m in maps], axis=1),
                               direct_real, rtol=1e-12, atol=1e-11)


def test_complete_reaction_transpose_matches_explicit_source_work(surfaces):
    _, points, tri, _ = surfaces
    rng = np.random.default_rng(315)
    poisson = SurfacePoissonPhiInc(points, tri)
    centers = points[tri].mean(axis=1)
    area = poisson._areas
    # Arbitrary complex total tangent field tests all three Cartesian maps;
    # it need not be a single-valued scalar gradient.
    total_h = rng.normal(size=centers.shape)+1j*rng.normal(size=centers.shape)
    total_h -= np.einsum('td,td->t', total_h, poisson._n_hat)[:, None]*poisson._n_hat
    body_j = np.cross(poisson._n_hat, total_h)
    z = np.where(np.arange(len(tri)) % 2, 2e-4+3e-4j, 8e-5+1e-4j)
    source_centers = centers[::7]+[.09, 0, .03]
    source_area = area[::7]
    source_j = rng.normal(size=source_centers.shape)+1j*rng.normal(size=source_centers.shape)
    omega = 2*np.pi*50000
    hload = electric_incident_vertex_load(poisson, z, total_h)
    incident_h = H_from_surface_J_complex(points, source_centers, source_area, source_j)
    incident_phi, _ = poisson(incident_h)
    incident_ht = -np.einsum('tkd,tk->td', poisson._g, incident_phi[tri])
    incident_a = A_from_surface_J(centers, source_centers, source_area, source_j)
    explicit = np.sum(area*z*np.einsum('td,td->t', incident_ht, total_h))
    explicit += 1j*omega*np.sum(area*np.einsum('td,td->t', incident_a, body_j))
    load = surface_source_reaction_load(source_centers, source_area, points,
        centers, area, body_j, hload, omega, pair_budget=17)
    transpose = np.sum(source_j*load)
    assert transpose == pytest.approx(explicit, rel=1e-12, abs=1e-20)
    whole = surface_source_reaction_load(source_centers, source_area, points,
        centers, area, body_j, hload, omega)
    np.testing.assert_allclose(load, whole, rtol=1e-13, atol=1e-21)


def test_electric_transpose_uses_bilinear_phasor_work(surfaces):
    _, points, tri, _ = surfaces
    poisson = SurfacePoissonPhiInc(points, tri)
    rng = np.random.default_rng(416)
    incident = rng.normal(size=points.shape)+1j*rng.normal(size=points.shape)
    phi, _ = poisson(incident)
    total = rng.normal(size=(len(tri), 3))+1j*rng.normal(size=(len(tri), 3))
    z = 1e-4+2e-4j
    load = electric_incident_vertex_load(poisson, z, total)
    projected = -np.einsum('tkd,tk->td', poisson._g, phi[tri])
    direct = np.sum(poisson._areas*z*np.einsum('td,td->t', projected, total))
    assert np.sum(incident*load) == pytest.approx(direct, rel=1e-12, abs=1e-18)


def test_strong_output_uses_port_reaction_and_coil_loss_change(monkeypatch):
    from types import SimpleNamespace
    from radia.panels import calc_inductance as ci
    monkeypatch.setattr(ci, "_assemble_vacuum_output", lambda *_:
                        {"L_coil_nH": 99., "R_coil_mOhm": 99.})
    strong = dict(P_total=1., H_t_rms=2., L_air=2e-9, L_total=3e-9,
        Delta_L=1e-9, R_air=.2, R_total=2.6, Delta_R=2.4,
        body_reaction_power_W=1., body_power_balance_relative_error=0.,
        body_residual=1e-15, port_power_W=1.3, coil_loss_W=.3,
        coil_loss_air_W=.1, coil_loss_change_W=.2, converged=True,
        coupling_residual=1e-9, iterations=5, n_phi_wp=8, n_J_coil=2,
        wp_mesh_nv=8, wp_mesh_n_tris=12, Z_s=1e-4+1e-4j,
        delta_wp_mm=1., t_wp_mesh_s=0., t_coupled_solve_s=0.)
    out = ci._assemble_strong_output(SimpleNamespace(current=1., coil_solver="peec"),
                                     {"source_type": "filament"}, strong)
    assert out["L_coil_nH"] == pytest.approx(2.)
    assert out["L_total_nH"] == pytest.approx(3.)
    assert out["R_total_mOhm"] == pytest.approx(2600.)
    assert out["delta_R_mOhm"] == pytest.approx(2400.)
    assert out["delta_R_mOhm"] != pytest.approx(2*out["P_wp_W"]*1000)
    assert out["port_power_W"] == pytest.approx(out["body_reaction_power_W"]+out["coil_loss_W"])


@pytest.mark.parametrize("coil_model", ["ideal-surface", "resistive-filaments"])
def test_coupled_iteration_matches_independent_three_loop_circuit(coil_model):
    from radia.bem_complete_reaction import (MU_0, _iterate_complete_current,
                                             surface_coil_reaction_rhs)
    from radia.peec_bundle import solve_loop_bundle
    omega = 2*np.pi*50000
    inductance = np.array([[.4, .1], [.1, .5]])*1e-6
    mutual = np.array([.15, .12])*1e-6
    resistance = (np.zeros((2, 2)) if coil_model == "ideal-surface"
                  else np.diag([.001, .003]))
    body_z = .002+1j*omega*.3e-6
    reflected = omega**2*np.outer(mutual, mutual)/body_z
    coil_z = resistance+1j*omega*inductance
    # Independent physical coil/body/port system, before eliminating body current.
    block = np.zeros((4, 4), dtype=complex)
    block[:2, :2] = coil_z
    block[:2, 2] = 1j*omega*mutual
    block[2, :2] = 1j*omega*mutual
    block[2, 2] = body_z
    block[:2, 3] = -1
    block[3, :2] = 1
    exact = np.linalg.solve(block, [0., 0., 0., 1.])

    if coil_model == "ideal-surface":
        saddle = np.zeros((3, 3), dtype=complex)
        saddle[:2, :2] = inductance/MU_0
        saddle[:2, 2] = saddle[2, :2] = 1

        def mapped(emf):
            rhs = np.r_[surface_coil_reaction_rhs(emf, omega), 1.]
            return np.linalg.solve(saddle, rhs)[:2]
    else:
        def mapped(emf):
            return solve_loop_bundle(resistance, inductance, omega/(2*np.pi), emf=emf)[0]

    def response(currents):
        return dict(body=dict(residual=0.), emf=reflected @ currents)

    currents, state, residual, _ = _iterate_complete_current(mapped(np.zeros(2)),
        response, mapped, max_iter=60, tol=1e-12, relax=.7)
    np.testing.assert_allclose(currents, exact[:2], rtol=1e-10, atol=1e-12)
    assert residual <= 1e-12
    body_current = -1j*omega*(mutual @ currents)/body_z
    body_heat = .5*body_z.real*abs(body_current)**2
    body_reaction = .5*np.vdot(currents, state['emf']).real
    assert body_reaction == pytest.approx(body_heat, rel=1e-12)
    port = np.vdot(currents, coil_z @ currents+state['emf'])
    assert port == pytest.approx(exact[3], rel=1e-10)
    coil_loss = .5*np.vdot(currents, resistance @ currents).real
    assert .5*port.real == pytest.approx(body_heat+coil_loss, rel=1e-12)
    np.testing.assert_array_equal(reflected, reflected.T)
    with pytest.raises(RuntimeError, match="did not converge"):
        _iterate_complete_current(mapped(np.zeros(2)), response, mapped,
                                  max_iter=2, tol=1e-15, relax=.1)


def test_uniform_field_linear_vector_potential_exact_closed_surface_work(surfaces):
    """Uniform H / linear A makes both magnetic midpoint rules exact."""
    from radia.bem_complete_reaction import MU_0
    _, points, tri, _ = surfaces
    poisson = SurfacePoissonPhiInc(points, tri)
    rng = np.random.default_rng(519)
    phi = rng.normal(size=len(points))+1j*rng.normal(size=len(points))
    field = np.array([.4+.1j, -.3+.2j, .7-.5j])
    centers = points[tri].mean(axis=1)
    a = .5*MU_0*np.cross(field, centers)
    ht = -np.einsum('tkd,tk->td', poisson._g, phi[tri])
    current = np.cross(poisson._n_hat, ht)
    ja = np.sum(poisson._areas*np.einsum('td,td->t', current, a))
    phib = np.sum(poisson._areas*phi[tri].mean(axis=1)
                  *np.einsum('td,d->t', poisson._n_hat, MU_0*field))
    assert ja == pytest.approx(phib, rel=1e-12, abs=1e-20)


@pytest.mark.parametrize('loop_dof', [False, True])
def test_frozen_coil_complete_body_matches_identical_weak_work(surfaces, loop_dof):
    """First coupled body state equals the corresponding fixed-coil solve."""
    from radia.bem_sibc_solver import ScalarBIESIBCSolver
    from radia.bem_complete_reaction import solve_complete_body
    from radia.bem_loop_extension import solve_loop_extended, A_from_filaments
    from radia.biot_savart import h_segments_batch
    mesh, points, tri, _ = surfaces
    angle = np.arange(48)*2*np.pi/48
    pp = np.c_[.06*np.cos(angle), .06*np.sin(angle), np.full(48, .005)]
    path = list(zip(pp, np.roll(pp, -1, axis=0)))
    current = .7+.2j
    omega = 2*np.pi*50000
    impedance = (1+1j)*np.sqrt(omega*4e-7*np.pi/(2*5.8e7))
    poisson = SurfacePoissonPhiInc(points, tri)
    incident_h = current*h_segments_batch(path, points)
    incident_phi, _ = poisson(incident_h, max_grad_residual=.1)
    def a_inc(x):
        return A_from_filaments(x, [path], [current])
    with ng.TaskManager():
        solver = ScalarBIESIBCSolver(mesh, order=1, assemble_dense=True,
            use_intree_bem=True, intree_geom_order=1,
            intree_singular_n_q=6, intree_regular_quad_degree=7)
        body = solve_complete_body(solver, poisson, incident_phi, impedance,
                                   omega, a_inc, loop_dof=loop_dof)
        if loop_dof:
            weak = solve_loop_extended(solver, incident_phi, impedance, omega, a_inc)
            np.testing.assert_allclose(body['field'], weak['H_t_tri'], rtol=1e-12, atol=1e-12)
            assert body['loop_metadata']['wp_loop_alpha'] == pytest.approx(weak['alpha'], rel=1e-12)
            weak_work = 1j*omega*weak['reaction_integral']
        else:
            weak = solver.solve(incident_phi, impedance, omega)
            weak_ht = -np.einsum('tkd,tk->td', poisson._g, weak['phi_vec'][tri])
            np.testing.assert_allclose(body['field'], weak_ht, rtol=1e-12, atol=1e-12)
            inc_ht = -np.einsum('tkd,tk->td', poisson._g, incident_phi[tri])
            weak_work = impedance*np.sum(poisson._areas*np.einsum('td,td->t', inc_ht, weak_ht))
            weak_work += 1j*omega*np.sum(poisson._areas*np.einsum(
                'td,td->t', np.cross(poisson._n_hat, weak_ht), a_inc(points[tri].mean(axis=1))))
        load = electric_incident_vertex_load(poisson, impedance, body['field'])
        strong_work = np.sum(incident_h*load)+1j*omega*np.sum(
            poisson._areas*np.einsum('td,td->t', body['current'], a_inc(points[tri].mean(axis=1))))
        assert strong_work == pytest.approx(weak_work, rel=1e-12, abs=1e-18)


def test_peec_first_coil_iterate_matches_fixed_coil_loop_reaction(surfaces, monkeypatch):
    """Capture the actual first response callback, before any coil update."""
    from radia import bem_complete_reaction as reaction
    from radia.peec_coupled_bem_solver import CoupledPEECBEMSolver
    from radia.bem_loop_extension import solve_loop_extended, A_from_filaments
    from radia.biot_savart import h_segments_batch
    mesh, points, tri, _ = surfaces
    angle = np.arange(48)*2*np.pi/48
    p = np.c_[.06*np.cos(angle), .06*np.sin(angle), np.full(48, .005)]
    path = list(zip(p, np.roll(p, -1, axis=0)))
    omega = 2*np.pi*50000
    z = (1+1j)*np.sqrt(omega*4e-7*np.pi/(2*5.8e7))
    captured = {}
    class FirstStateCaptured(Exception):
        pass
    def capture(initial, response, map_current, **kwargs):
        captured['current'] = initial.copy()
        captured['state'] = response(initial)
        raise FirstStateCaptured
    monkeypatch.setattr(reaction, '_iterate_complete_current', capture)
    with ng.TaskManager():
        solver = CoupledPEECBEMSolver([path], np.array([[.002]]),
                                    np.array([[3e-7]]), mesh)
        with pytest.raises(FirstStateCaptured):
            solver.solve(z, omega, loop_dof=True)
        current = captured['current'][0]
        poisson = solver._phi_poisson
        phi, _ = poisson(current*h_segments_batch(path, points), max_grad_residual=.1)
        weak = solve_loop_extended(solver.wp_solver, phi, z, omega,
            lambda x: A_from_filaments(x, [path], [current]))
    assert current*captured['state']['emf'][0] == pytest.approx(
        1j*omega*weak['reaction_integral'], rel=1e-12, abs=1e-18)
    np.testing.assert_allclose(captured['state']['body']['field'], weak['H_t_tri'],
                               rtol=1e-12, atol=1e-12)
