import numpy as np
import pytest

import radia.vim as vim


SIGMA = 1.0e7
THICKNESS_M = 0.5e-3
MU = 100.0 * vim.MU0


def _volumetric_gate() -> vim.EddySIBCApplicability:
    return vim.EddySIBCApplicability(
        frequency_hz=100.0,
        sigma=SIGMA,
        characteristic_thickness_m=THICKNESS_M,
        mu=MU,
    )


def _single_mode_system(
    gate: vim.EddySIBCApplicability,
) -> vim.CoupledHDivHybridVIMSystem:
    points = np.array([[0.0, 0.0, 0.0]])
    weights = np.array([1.0])
    modes = np.array([[[1.0, 0.0, 0.0]]])
    current = vim.VolumeCurrentBasis(points, weights, modes, names=("j0",))
    magnetization = vim.MagnetizationBasis(
        points,
        weights,
        modes,
        names=("m0",),
    )
    eddy = vim.HybridVIMSystem(
        resistance=np.array([[2.0]]),
        inductance=np.array([[0.4]]),
        surface_mass=np.zeros((1, 1)),
        basis_names=("j0",),
        blocks={"volume": (0, 1)},
    )
    return vim.CoupledHDivHybridVIMSystem(
        magnetization_basis=magnetization,
        eddy_system=eddy,
        eddy_bases=(current,),
        coupling=np.zeros((1, 1)),
        magnetic_operator=np.array([[3.0]]),
        magnetic_rhs=np.array([1.5]),
        eddy_rhs=np.array([0.75]),
        sibc_applicability=gate,
        conductivity=SIGMA,
    )


def test_frequency_route_reports_same_regime_and_crossing():
    gate = _volumetric_gate()

    same_route = gate.frequency_route(200.0)
    crossing = gate.frequency_route(10_000.0)

    assert same_route["assembled_frequency_hz"] == 100.0
    assert same_route["requested_frequency_hz"] == 200.0
    assert same_route["assembled_model"] == "volumetric"
    assert same_route["requested_model"] == "volumetric"
    assert same_route["route_compatible"] is True
    assert same_route["requested_thickness_to_skin_depth"] == pytest.approx(
        gate.at_frequency(200.0).thickness_to_skin_depth
    )
    assert same_route["requested_curvature_radius_to_skin_depth"] is None
    assert crossing["assembled_model"] == "volumetric"
    assert crossing["requested_model"] == "sibc"
    assert crossing["route_compatible"] is False


def test_frequency_solve_allows_same_route_but_rejects_stale_model():
    system = _single_mode_system(_volumetric_gate())

    solved = system.solve_frequency(200.0, surface_impedance=0.0)

    assert np.all(np.isfinite(solved.reduced_solution))
    with pytest.raises(
        ValueError,
        match="frequency sweep crosses the assembled eddy-model route",
    ):
        system.solve_frequency(10_000.0, surface_impedance=0.0)


def test_coupled_diagnostics_preserve_frequency_gate():
    gate = _volumetric_gate()
    system = _single_mode_system(gate)

    diagnostics = system.diagnostics()

    assert diagnostics["sibc_applicability"]["selected_model"] == "volumetric"
    assert diagnostics["sibc_applicability"]["frequency_hz"] == 100.0


def test_coupled_gate_conductivity_must_match_assembled_system():
    with pytest.raises(
        ValueError,
        match="conductivity must match sibc_applicability.sigma",
    ):
        _single_mode_system(
            vim.EddySIBCApplicability(
                frequency_hz=100.0,
                sigma=2.0 * SIGMA,
                characteristic_thickness_m=THICKNESS_M,
                mu=MU,
            )
        )


def test_coupling_factory_preserves_frequency_gate():
    gate = _volumetric_gate()
    points = np.array([[0.0, 0.0, 0.0]])
    weights = np.array([1.0])
    modes = np.array([[[1.0, 0.0, 0.0]]])
    current = vim.VolumeCurrentBasis(points, weights, modes, names=("j0",))
    magnetization = vim.MagnetizationBasis(
        points,
        weights,
        modes,
        names=("m0",),
    )
    eddy = vim.AssembleHybridVIM(
        current,
        sigma=SIGMA,
        kernel_epsilon=0.1,
    )

    coupled = vim.CoupleHybridVIMWithHDivMMM(
        magnetization,
        eddy,
        (current,),
        magnetic_operator=np.array([[3.0]]),
        magnetic_rhs=np.array([1.0]),
        eddy_rhs=np.array([1.0]),
        sibc_applicability=gate,
        conductivity=SIGMA,
        kernel_epsilon=0.1,
    )

    assert coupled.sibc_applicability is gate
    assert coupled.diagnostics()["sibc_applicability"]["selected_model"] == (
        "volumetric"
    )
