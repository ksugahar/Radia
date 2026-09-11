from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from radia.force import integrate_time_average_lorentz_force
from radia.maglev.ecb import lorentz


class _Model:
    def __init__(self, coefficients):
        self.coefficients = np.asarray(coefficients, dtype=complex)
        self.calls = []

    def solve_vector_potential_drive(self, s, drive):
        self.calls.append((s, drive))
        return self.coefficients


def _sampled_case():
    basis = SimpleNamespace(
        modes=np.asarray(
            [
                [[1.0 + 1.0j, 0.0, 0.0], [0.0, 2.0, 0.0]],
                [[0.0, 1.0, 0.0], [3.0j, 0.0, 0.0]],
            ],
            dtype=complex,
        ),
        weights=np.asarray([0.25, 0.75]),
        points=np.asarray([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]]),
    )
    magnetic_flux_density = np.asarray(
        [[0.0, 0.0, 2.0 - 1.0j], [0.0, 0.0, 0.5 + 0.25j]],
        dtype=complex,
    )
    return _Model([0.5 - 0.25j, -0.75j]), basis, magnetic_flux_density


def test_hcurl_vim_adapter_uses_the_common_peak_phasor_integrator():
    model, basis, magnetic_flux_density = _sampled_case()
    force = lorentz.compute_lorentz_force_via_hcurl_vim(
        model, basis, magnetic_flux_density, 4.0j, drive=2.5
    )
    current = np.einsum("a,aik->ik", model.coefficients, basis.modes)
    expected = integrate_time_average_lorentz_force(
        current,
        magnetic_flux_density,
        basis.weights,
        amplitude="peak",
    )
    assert force == pytest.approx(expected)
    assert model.calls == [(4.0j, 2.5)]


def test_hcurl_vim_result_records_action_reaction_and_3d_provenance():
    model, basis, magnetic_flux_density = _sampled_case()
    result = lorentz.compute_lorentz_force_result_via_hcurl_vim(
        model,
        basis,
        magnetic_flux_density,
        4.0j,
        frame="plate",
    )
    assert result["schema"] == "radia.maglev-force-pair/v1"
    assert result["validated_against_3d"] is True
    assert np.add(
        result["conductor"]["force_N"], result["source_pm"]["force_N"]
    ).tolist() == pytest.approx([0.0, 0.0, 0.0])
    assert result["conductor"]["method"] == (
        "hcurl_vim_time_average_lorentz_body_force"
    )
    assert result["conductor"]["torque_Nm"] is not None
    assert np.add(
        result["conductor"]["torque_Nm"], result["source_pm"]["torque_Nm"]
    ).tolist() == pytest.approx([0.0, 0.0, 0.0])
    assert result["action_reaction_torque_residual_Nm"] == [0.0, 0.0, 0.0]
    assert result["conductor"]["phasor_amplitude"] == "peak"
    assert result["conductor"]["frame"] == "plate"


def _verified_call(**kwargs):
    return lorentz.compute_lorentz_force_via_foster_verified(
        None,
        np.arange(4.0),
        None,
        None,
        1.0,
        1.0,
        1.0j,
        1.0,
        1.0,
        0.0,
        **kwargs,
    )


def test_verified_foster_keeps_a_converged_reduction(monkeypatch):
    monkeypatch.setattr(
        lorentz, "compute_lorentz_force_via_foster", lambda *a, **k: (0, 0, -1.005)
    )
    monkeypatch.setattr(
        lorentz, "compute_lorentz_force_via_direct_scalar", lambda *a, **k: (0, 0, -1.0)
    )
    with pytest.warns(RuntimeWarning, match="not a quantitatively validated"):
        result = _verified_call(rtol=0.01)
    assert result["foster_converged"] is True
    assert result["backend"] == "foster"
    assert result["force_N"] == [0.0, 0.0, -1.005]
    assert result["foster_mode_count"] == 4
    assert result["validated_against_3d"] is False


def test_verified_foster_falls_back_when_high_frequency_modes_are_insufficient(
    monkeypatch,
):
    monkeypatch.setattr(
        lorentz, "compute_lorentz_force_via_foster", lambda *a, **k: (0, 0, -0.7)
    )
    monkeypatch.setattr(
        lorentz, "compute_lorentz_force_via_direct_scalar", lambda *a, **k: (0, 0, -1.0)
    )
    with pytest.warns(RuntimeWarning):
        result = _verified_call(rtol=0.01)
    assert result["foster_converged"] is False
    assert result["backend"] == "direct_scalar_fallback"
    assert result["force_N"] == [0.0, 0.0, -1.0]
    assert result["relative_difference"] == pytest.approx(0.3)


def test_verified_foster_can_fail_loudly_instead_of_falling_back(monkeypatch):
    monkeypatch.setattr(
        lorentz, "compute_lorentz_force_via_foster", lambda *a, **k: (0, 0, -0.7)
    )
    monkeypatch.setattr(
        lorentz, "compute_lorentz_force_via_direct_scalar", lambda *a, **k: (0, 0, -1.0)
    )
    with pytest.warns(RuntimeWarning), pytest.raises(RuntimeError, match="did not converge"):
        _verified_call(rtol=0.01, fallback_to_direct=False)


@pytest.mark.parametrize("name,value", [("rtol", -1.0), ("atol_N", -1.0)])
def test_verified_foster_rejects_invalid_tolerances(name, value):
    with pytest.raises(ValueError, match=name):
        _verified_call(**{name: value})


def test_hcurl_vim_adapter_rejects_incompatible_mode_shape():
    model = _Model([1.0, 2.0])
    basis = SimpleNamespace(
        modes=np.zeros((1, 2, 3)),
        weights=np.ones(2),
        points=np.zeros((2, 3)),
    )
    with pytest.raises(ValueError, match="n_modes"):
        lorentz.compute_lorentz_force_via_hcurl_vim(
            model, basis, np.zeros((2, 3)), 1.0j
        )


def test_verified_foster_live_high_frequency_falls_back_to_direct_scalar():
    ngsolve = pytest.importorskip("ngsolve")
    from ngsolve.meshes import MakeStructured3DMesh

    from radia.maglev.mixed_galerkin.alpha import _dirichlet_eigenmodes

    base = MakeStructured3DMesh(
        hexes=True,
        nx=10,
        ny=4,
        nz=1,
        mapping=lambda x, y, z: (0.10 * (x - 0.5), 0.06 * (y - 0.5),
                                 0.002 * (z - 1.0)),
    )
    for index in range(len(set(base.GetBoundaries()))):
        base.ngmesh.SetBCName(index, "outer")
    mesh = ngsolve.Mesh(base.ngmesh)
    with ngsolve.TaskManager():
        lam, vecs, _mass, free, _fes, _volume = _dirichlet_eigenmodes(
            mesh, 20, "outer"
        )
        with pytest.warns(RuntimeWarning):
            result = lorentz.compute_lorentz_force_via_foster_verified(
                mesh,
                lam,
                vecs,
                free,
                3.5e7,
                lorentz.MU_0,
                2.0j * np.pi * 5000.0,
                10.0,
                0.020,
                0.0,
                rtol=1.0e-6,
            )
    assert result["backend"] == "direct_scalar_fallback"
    assert result["foster_converged"] is False
    assert result["relative_difference"] > 1.0e-6
