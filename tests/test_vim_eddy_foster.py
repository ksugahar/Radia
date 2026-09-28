"""Foster modal form of the reduced HCurl eddy-current model.

Every check compares against an independent dense solve of
``(R + sL + Zs Ms) c = P u`` so the modal form is verified, not assumed.
"""

import json

import numpy as np
import pytest

from radia import vim


def _spd(rng, n, shift):
    a = rng.standard_normal((n, n))
    return a @ a.T + shift * np.eye(n)


def _system(seed=0, n=6, ports=2):
    rng = np.random.default_rng(seed)
    resistance = _spd(rng, n, 0.1)
    inductance = _spd(rng, n, 0.5) * 1.0e-3
    port = rng.standard_normal((n, ports))
    return resistance, inductance, port


def _model(resistance, inductance, port, surface_mass=None):
    n = resistance.shape[0]
    return vim.HCurlEddyFosterModel(
        resistance=resistance,
        inductance=inductance,
        surface_mass=np.zeros((n, n)) if surface_mass is None else surface_mass,
        port_rhs=port,
    )


@pytest.mark.parametrize("s", [1.0, 2j * np.pi * 50.0, 2j * np.pi * 5.0e3, 3.0 + 40.0j])
def test_foster_solve_equals_direct_reduced_solve(s):
    resistance, inductance, port = _system()
    model = _model(resistance, inductance, port)
    drive = np.array([1.5, -0.25])
    direct = np.linalg.solve(resistance + s * inductance, port @ drive)
    np.testing.assert_allclose(model.solve(s, drive), direct, rtol=1e-10, atol=0)
    current = np.array([2.0, 1.0])
    np.testing.assert_allclose(
        model.solve_vector_potential_drive(s, current),
        np.linalg.solve(resistance + s * inductance, -s * port @ current),
        rtol=1e-10, atol=0)
    np.testing.assert_allclose(
        model.port_admittance(s),
        port.T @ np.linalg.solve(resistance + s * inductance, port),
        rtol=1e-10, atol=0)


def test_foster_state_space_is_diagonal_passive_and_matches_transfer_function():
    resistance, inductance, port = _system(seed=3)
    model = _model(resistance, inductance, port)
    ss = model.derivative_input_state_space()
    a = ss["A"]
    np.testing.assert_array_equal(a, np.diag(np.diag(a)))
    assert np.all(np.diag(a) < 0.0)
    rates = np.linalg.eigvals(np.linalg.solve(inductance, resistance)).real
    np.testing.assert_allclose(np.sort(model.decay_rates), np.sort(rates), rtol=1e-10)
    for s in (0.5, 2j * np.pi * 200.0):
        transfer = ss["C"] @ np.linalg.solve(s * np.eye(a.shape[0]) - a, ss["B"]) + ss["D"]
        direct = port.T @ np.linalg.solve(resistance + s * inductance, port)
        np.testing.assert_allclose(transfer, direct, rtol=1e-10, atol=0)
    diagnostics = model.diagnostics()
    assert diagnostics["passive"] is True
    assert diagnostics["form"] == "foster_modal"
    assert diagnostics["min_decay_rate"] > 0.0


def test_foster_modal_force_operator_matches_reduced_force():
    resistance, inductance, port = _system(seed=5)
    model = _model(resistance, inductance, port)
    rng = np.random.default_rng(11)
    operator = rng.standard_normal((3, model.state_order, model.port_count))
    s = 2j * np.pi * 50.0
    current = np.array([20.0, -5.0])
    c = model.solve_vector_potential_drive(s, current)
    z = np.linalg.solve(model.modes, c)
    reduced = 0.5 * np.real(np.einsum("kab,a,b->k", operator, c, np.conj(current)))
    modal = 0.5 * np.real(np.einsum(
        "kjb,j,b->k", model.modal_force_operator(operator), z, np.conj(current)))
    np.testing.assert_allclose(modal, reduced, rtol=1e-10, atol=1e-14)


def test_sibc_block_uses_the_direct_solve_and_blocks_state_space():
    resistance, inductance, port = _system(seed=7)
    rng = np.random.default_rng(2)
    surface = _spd(rng, resistance.shape[0], 0.2) * 1.0e-2
    model = _model(resistance, inductance, port, surface_mass=surface)
    s = 2j * np.pi * 1.0e3
    zs = 0.3 + 0.3j
    drive = np.array([1.0, 0.0])
    direct = np.linalg.solve(resistance + s * inductance + zs * surface, port @ drive)
    np.testing.assert_allclose(model.solve(s, drive, surface_impedance=zs), direct,
                               rtol=1e-10, atol=0)
    with pytest.raises(ValueError, match="SIBC termination"):
        model.derivative_input_state_space()


def test_foster_form_rejects_non_passive_or_non_hermitian_matrices():
    resistance, inductance, port = _system(seed=9)
    with pytest.raises(ValueError, match="positive definite"):
        _model(resistance, -inductance, port)
    skew = inductance.copy()
    skew[0, 1] += 1.0
    with pytest.raises(ValueError, match="Hermitian"):
        _model(resistance, skew, port)


def test_foster_exchange_and_family_share_one_mode_set(tmp_path):
    resistance, inductance, port = _system(seed=13)
    rng = np.random.default_rng(4)
    models, operators = [], []
    for scale in (1.0, 0.8, 0.6):
        models.append(_model(resistance, inductance, scale * port))
        operators.append(rng.standard_normal((3, 6, 2)))
    single = vim.ExportHCurlEddyFosterJSON(models[0], tmp_path / "single.json",
                                           force_operator=operators[0])
    data = json.loads(single.read_text(encoding="utf-8"))
    assert data["schema"] == "radia.hcurl.eddy_foster.exchange.v1"
    np.testing.assert_allclose(data["decay_rates"]["values"], models[0].decay_rates)
    modal = np.asarray(data["arrays"]["modal_port_rhs"]["values"]).reshape(6, 2)
    np.testing.assert_allclose(modal, models[0].modal_port_rhs)

    family = vim.ExportHCurlEddyFosterFamilyJSON(
        [{"height_m": h, "model": m, "force_operator": k}
         for h, m, k in zip((0.012, 0.011, 0.010), models, operators)],
        tmp_path / "family.json")
    fam = json.loads(family.read_text(encoding="utf-8"))
    assert fam["schema"] == "radia.hcurl.eddy_foster.family.v1"
    assert [snap["height_m"] for snap in fam["snapshots"]] == [0.010, 0.011, 0.012]
    reference_modes = models[0].modes
    for snap, model in zip(fam["snapshots"], models[::-1]):
        modal = np.asarray(snap["arrays"]["modal_port_rhs"]["values"]).reshape(6, 2)
        np.testing.assert_allclose(modal, reference_modes.T @ model.port_rhs, rtol=1e-12)

    moved = _model(resistance * 1.01, inductance, port)
    with pytest.raises(ValueError, match="shared R/L"):
        vim.ExportHCurlEddyFosterFamilyJSON(
            [{"height_m": 0.0, "model": models[0]}, {"height_m": 1.0, "model": moved}],
            tmp_path / "bad.json")


def test_foster_model_from_vim_system_copies_the_reduced_blocks():
    resistance, inductance, port = _system(seed=17)
    n = resistance.shape[0]
    system = vim.HybridVIMSystem(
        resistance=resistance,
        inductance=inductance,
        surface_mass=np.zeros((n, n)),
        basis_names=tuple(f"m{i}" for i in range(n)),
        blocks={"volume": (0, n)},
    )
    model = vim.HCurlEddyFosterModelFromVIM(system, port)
    assert model.basis_names == system.basis_names
    assert model.blocks == {"volume": (0, n)}
    s = 2j * np.pi * 50.0
    np.testing.assert_allclose(
        model.solve(s, np.ones(2)),
        np.linalg.solve(system.impedance(s), port @ np.ones(2)), rtol=1e-10, atol=0)
