"""Fast policy and argument validation; no NGSolve solve is permitted."""
import builtins
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
    real_import = builtins.__import__
    def guarded_import(name, *args, **kwargs):
        if name.rsplit(".", 1)[-1] in {"fem_bem_coupling", "ngsolve"}:
            pytest.fail("invalid request imported a numerical solver")
        return real_import(name, *args, **kwargs)
    monkeypatch.setattr(builtins, "__import__", guarded_import)
    with pytest.raises(ValueError, match=message):
        validation.run_kelvin_open_boundary_validation(case)
