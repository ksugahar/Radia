"""Face-weighted genus-1 work, heating and constant-panel reduction checks."""
import importlib.util
from pathlib import Path

import ngsolve as ng
import numpy as np
import pytest

from radia.bem_loop_extension import solve_loop_extended
from radia.bem_sibc_solver import ScalarBIESIBCSolver
from radia.surface_impedance import PanelSurfaceImpedance


@pytest.fixture(scope="module")
def ring():
    path = Path(__file__).resolve().parents[1] / "validation_test" / "induction_heating" / "run_loop_work_ring.py"
    spec = importlib.util.spec_from_file_location("ring_example", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    mesh, points, triangles = module.ring_mesh(16)
    with ng.TaskManager():
        bem = ScalarBIESIBCSolver(mesh, order=1, assemble_dense=True,
            use_intree_bem=True, intree_geom_order=1,
            intree_singular_n_q=6, intree_regular_quad_degree=7)
    areas = .5*np.linalg.norm(np.cross(points[triangles[:, 1]]-points[triangles[:, 0]],
                                       points[triangles[:, 2]]-points[triangles[:, 0]]), axis=1)
    mu0, omega = 4e-7*np.pi, 2*np.pi*5e4
    z = (1+1j)*np.sqrt(omega*mu0/(2*5.8e7))

    def a_inc(p):
        return .5*mu0*np.c_[-p[:, 1], p[:, 0], np.zeros(len(p))]

    angles = np.arange(96)*2*np.pi/96
    carrier = .03*np.c_[np.cos(angles), np.sin(angles), np.zeros(96)]
    return bem, points, triangles, areas, z, omega, a_inc, carrier


def solve(ring, impedance):
    bem, points, _, _, _, omega, a_inc, carrier = ring
    with ng.TaskManager():
        return solve_loop_extended(bem, -points[:, 2].astype(complex), impedance,
            omega, a_inc, section_anchor=(.03, 0), carrier_ring=carrier)


def test_constant_panels_reduce_to_scalar(ring):
    _, _, triangles, _, z, _, _, _ = ring
    scalar = solve(ring, z)
    panel = solve(ring, PanelSurfaceImpedance(np.full(len(triangles), z)))
    assert panel["loop_impedance_assembly"] == "scalar-equal-faces"
    np.testing.assert_allclose(panel["H_t_tri"], scalar["H_t_tri"], rtol=1e-10, atol=1e-12)
    for key in ["alpha", "P_total", "reaction_integral", "P_frozen"]:
        assert panel[key] == pytest.approx(scalar[key], rel=1e-10, abs=1e-20)


def test_nearly_constant_panels_still_use_weighted_assembly(ring):
    _, _, triangles, _, z, _, _, _ = ring
    values = np.full(len(triangles), z)
    values[0] *= 1+1e-12
    scalar = solve(ring, z)
    panel = solve(ring, PanelSurfaceImpedance(values))
    assert panel["loop_impedance_assembly"] == "weighted-panel"
    np.testing.assert_array_equal(panel["Z_s_per_panel"], values)
    # Fine-mesh derivative recovery amplified assembly rounding by ~1e6;
    # this perturbed-coefficient check has a conditioning allowance of 1e-9.
    assert np.linalg.norm(panel["H_t_tri"]-scalar["H_t_tri"])/np.linalg.norm(scalar["H_t_tri"]) < 1e-9
    for key in ["alpha", "P_total", "reaction_integral"]:
        assert panel[key] == pytest.approx(scalar[key], rel=1e-9, abs=1e-20)


@pytest.mark.parametrize("reactive_half", [False, True])
def test_face_jump_weights_field_heat_and_reaction(ring, reactive_half):
    bem, points, triangles, areas, z, omega, _, _ = ring
    # Inner/outer minor-angle halves; both interfaces lie on mesh edges.
    cell = np.arange(len(triangles)) // (2*16)
    values = np.where((cell >= 1) & (cell < 3), .25*z, 3*z)
    if reactive_half:
        values[(cell >= 1) & (cell < 3)] = .25j*z.imag
    field = PanelSurfaceImpedance(values)
    out = solve(ring, field)
    np.testing.assert_array_equal(out["Z_s_per_panel"], values)
    expected_heat = .5*values.real*np.sum(abs(out["H_t_tri"])**2, axis=1)
    np.testing.assert_allclose(out["q_tri"], expected_heat, rtol=1e-13)
    assert out["P_total"] == pytest.approx(areas @ expected_heat, rel=1e-12)
    with ng.TaskManager():
        plain = bem.solve(-points[:, 2].astype(complex), field, omega)
    assert out["P_frozen"] == pytest.approx(plain["P_density"]*plain["area"], rel=1e-10)
    averaged = solve(ring, np.average(values, weights=areas))
    difference = np.linalg.norm(out["H_t_tri"]-averaged["H_t_tri"])/np.linalg.norm(out["H_t_tri"])
    assert difference > 1e-4  # Fail if the panel field is silently averaged.
    assert out["linear_residual_rel"] < 1e-6
    assert out["faraday_residual_rel"] < 1e-6
    from radia.workpiece_surface import _check_sibc_reaction_power
    assert _check_sibc_reaction_power(out["P_total"], out["P_reaction"]) < .1
    # The output amplitude must retain alpha even on zero-real-impedance faces.
    from radia.panels.calc_inductance import _solved_panel_tangential_field
    amplitude = _solved_panel_tangential_field(bem, {"loop_H_t_tri":out["H_t_tri"]})
    np.testing.assert_allclose(amplitude, np.linalg.norm(out["H_t_tri"], axis=1), rtol=1e-14)
    if reactive_half:
        assert np.all(out["q_tri"][values.real == 0] == 0)
        assert np.all(amplitude[values.real == 0] > 0)


@pytest.mark.parametrize("wrong", ["untagged", "wrong_count"])
def test_ambiguous_or_misaligned_panel_values_fail(ring, wrong):
    _, _, triangles, _, z, _, _, _ = ring
    value = (np.full(len(triangles), z) if wrong == "untagged"
             else PanelSurfaceImpedance(np.full(len(triangles)-1, z)))
    with pytest.raises(ValueError, match="tagged face|one impedance"):
        solve(ring, value)


@pytest.mark.parametrize("backend",["intree-dense","hacapk"])
def test_prescribed_panels_allow_explicit_loop_before_source_build(monkeypatch,backend):
    from radia.panels import calc_inductance as ci
    args = ci.build_argparser().parse_args([
        "--coil-solver", "peec", "--coil-step", "unused.step",
        "--vol", "unused.vol", "--frequency", "50000", "--sigma", "5.8e7",
        "--panel-zs-file", "unused.json", "--wp-loop-dof", "on",
        "--wp-bem-backend",backend])

    def reached_source(_):
        raise RuntimeError("test reached source construction")

    monkeypatch.setattr(ci, "_solve_coil_peec", reached_source)
    with pytest.raises(RuntimeError, match="test reached source construction"):
        ci.run_inductance(args)


@pytest.mark.parametrize("law", ["smooth", "jump"])
def test_augmented_gradient_quadratic_heat(ring, monkeypatch, law):
    from scipy.sparse import coo_matrix, diags
    from radia import bem_loop_work
    bem, points, tri, areas, z, _, _, _ = ring
    centers = points[tri].mean(axis=1)
    values = (z*(1+.5*centers[:, 2]/.003) if law == "smooth"
              else np.where(centers[:, 0] >= 0, 3*z, .25*z))
    original = bem_loop_work._loop_work_row
    captured = {}

    def capture(*args, **kwargs):
        # Geometry/lift inputs, not the solver's returned electric matrix.
        captured["map"], captured["open_tri"] = args[4], args[3]
        captured["grad"], captured["theta"] = args[5], args[8]
        return original(*args, **kwargs)

    monkeypatch.setattr(bem_loop_work, "_loop_work_row", capture)
    out = solve(ring, PanelSurfaceImpedance(values))
    vertex_map, open_tri = captured["map"], captured["open_tri"]
    grad, theta = captured["grad"], captured["theta"]
    ndof = bem.ndof
    rows = np.repeat(np.arange(len(tri)), 3)
    columns = vertex_map[open_tri].ravel()
    quadratic = 0j
    x = np.r_[out["phi_u"], out["alpha"]]
    weights = diags(areas*values.real)
    for component in range(3):
        ordinary = coo_matrix((grad[:, :, component].ravel(), (rows, columns)),
                              shape=(len(tri), ndof)).tocsr()
        lifted = np.einsum('tk,tk->t', grad[:, :, component], theta[open_tri])
        from scipy.sparse import hstack
        augmented = hstack([ordinary, coo_matrix(lifted[:, None])]).tocsr()
        matrix = augmented.T @ weights @ augmented
        quadratic += .5*np.vdot(x, matrix @ x)
    assert abs(quadratic.imag) < 1e-15*max(abs(quadratic.real), 1e-30)
    assert quadratic.real == pytest.approx(out["P_total"], rel=1e-12)
    assert quadratic.real == pytest.approx(areas @ out["q_tri"], rel=1e-12)
    from radia.workpiece_surface import _check_sibc_reaction_power
    assert _check_sibc_reaction_power(out["P_total"], out["P_reaction"]) < .1


def test_forced_fmm_matches_dense_and_refreshes_prepared_system(ring, monkeypatch):
    bem, _, triangles, _, z, omega, _, _ = ring
    from radia.bem_loop_extension import solve_loop_extended
    _, points, _, _, _, _, a_inc, carrier = ring
    field = PanelSurfaceImpedance(np.where(np.arange(len(triangles)) % 2, z, 2*z))

    def prepared():
        with ng.TaskManager():
            return solve_loop_extended(bem, -points[:, 2].astype(complex), field,
                omega, a_inc, section_anchor=(.03, 0), carrier_ring=carrier,
                _reuse_prepared=True)

    monkeypatch.setattr(bem, "loop_work_backend", "dense")
    direct = prepared()
    system = bem._prepared_loop_system
    monkeypatch.setattr(bem, "loop_work_backend", "fmm")
    compressed = prepared()
    assert bem._prepared_loop_system is not system
    np.testing.assert_allclose(compressed["H_t_tri"], direct["H_t_tri"], rtol=1e-7, atol=1e-9)
    assert compressed["P_total"] == pytest.approx(direct["P_total"], rel=1e-7)
    diagnostic = compressed["loop_work_diagnostics"]
    assert diagnostic["single_layer_backend"] == "ngsolve-galerkin-fmm"
    assert diagnostic["single_layer_native_storage_bytes"] is None
    assert diagnostic["single_layer_dense_array_bytes"] == 0
    assert diagnostic["single_layer_fmm_parameters"]["fmm_minorder"] == 20
    assert diagnostic["ngsolve_version"]
    again = prepared()
    assert bem._prepared_loop_system is not system
    np.testing.assert_array_equal(again["H_t_tri"], compressed["H_t_tri"])
