import math

import pytest

from radia.vim import (
    ConductiveSlab,
    parallel_slab_exact_response,
    refine_parallel_slab_reduced,
    slab_surface_model,
    solve_parallel_slab_reduced,
)


SLAB = ConductiveSlab(
    thickness_m=5.0e-4,
    relative_permeability=100.0,
    conductivity_s_per_m=1.0e7,
    area_m2=1.0e-3,
)


@pytest.mark.parametrize("frequency_hz", [100.0, 1.0e4, 1.0e6])
def test_parallel_slab_p1_converges_to_closed_form(frequency_hz):
    result = solve_parallel_slab_reduced(SLAB, frequency_hz, elements=128)
    expected = parallel_slab_exact_response(SLAB, frequency_hz)
    relative_error = abs(result["effective_relative_permeability"] - expected) / abs(expected)
    assert relative_error < 0.012
    assert result["normalized_algebraic_residual"] < 1.0e-10
    assert result["joule_loss_w"] >= 0.0


def test_refinement_ledger_resolves_high_frequency_skin_depth():
    ledger = refine_parallel_slab_reduced(SLAB, [100.0, 1.0e4, 1.0e6])
    assert ledger["levels"] == [64, 96, 128]
    assert ledger["final_max_relative_change"] < 0.01
    assert min(row["elements_per_skin_depth"] for row in ledger["final_rows"]) >= 4.0


def test_surface_model_rejects_half_space_sibc_for_thin_regime():
    low = slab_surface_model(SLAB, 100.0)
    assert low["selected"] == "volumetric"
    assert not low["finite_thickness_separation_verified"]
    with pytest.raises(ValueError, match="thickness/skin_depth"):
        slab_surface_model(SLAB, 100.0, requested="half_space_sibc")


def test_surface_model_allows_half_space_sibc_when_surfaces_are_separated():
    high = slab_surface_model(SLAB, 1.0e6, requested="half_space_sibc")
    assert high["selected"] == "half_space_sibc"
    assert high["finite_thickness_separation_verified"]
    assert high["thickness_to_skin_depth_ratio"] > 30.0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"thickness_m": 0.0, "relative_permeability": 1.0, "conductivity_s_per_m": 1.0},
        {"thickness_m": 1.0, "relative_permeability": math.nan, "conductivity_s_per_m": 1.0},
    ],
)
def test_conductive_slab_rejects_invalid_physical_contract(kwargs):
    with pytest.raises(ValueError):
        ConductiveSlab(**kwargs)


def test_reduced_solve_rejects_nonphysical_frequency_and_bad_mesh():
    with pytest.raises(ValueError, match="frequency_hz"):
        solve_parallel_slab_reduced(SLAB, 0.0)
    with pytest.raises(ValueError, match="elements"):
        solve_parallel_slab_reduced(SLAB, 100.0, elements=1)
