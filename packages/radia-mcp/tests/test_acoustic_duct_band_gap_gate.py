from __future__ import annotations

import copy
import json

from radia_mcp.radia_ngsolve.acoustic_duct_band_gap_gate import acoustic_duct_band_gap_gate
from radia_mcp.radia_ngsolve.server import acoustic_duct_band_gap_gate as mcp_gate


def _summary() -> dict:
    return {
        "model_contract": {
            "cell_pitch": 1.5,
            "inclusion_radius": 0.3,
            "duct_side": 1.0,
            "finite_crystal_cell_count": 5,
            "maximum_wavenumber": 2.9,
            "inclusion": "sound_soft_sphere",
            "confined_geometry": "rigid_duct_periodic_cell",
            "free_space_geometry": "finite_linear_chain",
            "same_inclusion_family": True,
        },
        "duct_result": {
            "empty_lattice_relative_error": 2.04e-4,
            "empty_duct_transparency_error": 4.6e-7,
            "gap_low": 2.519,
            "gap_high": 3.647,
            "first_band_minimum": 2.314,
            "maximum_passband_transmission": 0.267,
            "maximum_gap_transmission": 0.00191,
            "maximum_below_band_transmission": 0.000292,
            "pass_to_gap_contrast": 140.0,
            "replay_relative_error": 8.1e-9,
        },
        "free_space_control": {
            "wavenumbers": [0.6, 1.4, 2.2, 3.0],
            "minimum_insertion_loss_db": 3.47,
            "maximum_insertion_loss_db": 3.94,
            "insertion_loss_spread_db": 0.47,
        },
        "timing_breakdown_s": {"reference": 23.0, "duct": 1.0, "free": 22.0, "verify": 1.0},
    }


def test_accepts_confined_gap_and_free_space_negative_control():
    result = acoustic_duct_band_gap_gate(_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())
    assert json.loads(mcp_gate(json.dumps(_summary())))["status"] == "ok"


def test_rejects_false_gap_and_missing_confinement_contrast():
    bad = copy.deepcopy(_summary())
    bad["duct_result"]["maximum_gap_transmission"] = 0.2
    bad["duct_result"]["pass_to_gap_contrast"] = 1.3
    bad["free_space_control"]["maximum_insertion_loss_db"] = 12.0
    bad["free_space_control"]["insertion_loss_spread_db"] = 8.0
    result = acoustic_duct_band_gap_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["finite_crystal_attenuates_in_gap"] is False
    assert result["checks"]["pass_to_gap_contrast_is_resolved"] is False
    assert result["checks"]["free_space_chain_has_attenuation_but_no_deep_stop_band"] is False


def test_rejects_multimode_sweep_replay_and_timing_drift():
    bad = _summary()
    bad["model_contract"]["maximum_wavenumber"] = 4.0
    bad["duct_result"]["replay_relative_error"] = 1.0e-3
    bad["timing_breakdown_s"].pop("verify")
    result = acoustic_duct_band_gap_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["duct_sweep_stays_below_first_transverse_cutoff"] is False
    assert result["checks"]["fresh_reference_replays_saved_observables"] is False
    assert result["checks"]["exactly_four_timing_stages"] is False
