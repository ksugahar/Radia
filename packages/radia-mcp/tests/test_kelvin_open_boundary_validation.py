"""Fast policy and argument validation; no NGSolve solve is permitted."""
import pytest
from radia_mcp.radia_ngsolve import kelvin_open_boundary_validation as validation


@pytest.mark.parametrize(
    "case, message",
    [
        ({"method": "pml"}, "kelvin_transform"),
        ({"pml": True}, "PML"),
        ({"physics_regime": "time_harmonic_wave"}, "static-only"),
        ({"wave_boundary_inference": "reuse_for_waves"}, "forbidden"),
        ({"inner_radius": 1.0, "outer_radius": 1.0}, "smaller"),
        ({"offset": 1.5}, "disjoint"),
    ],
)
def test_invalid_static_policy_or_geometry_fails_before_any_solve(case, message, monkeypatch):
    def forbidden_solve(**kwargs):
        pytest.fail("invalid request reached a numerical solver")
    monkeypatch.setattr(validation, "kelvin_dtn_eigenvalue", forbidden_solve)
    monkeypatch.setattr(validation, "kelvin_twosphere_shell_dipole", forbidden_solve)
    with pytest.raises(ValueError, match=message):
        validation.run_kelvin_open_boundary_validation(case)
