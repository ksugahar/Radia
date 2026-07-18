"""Regression coverage for identification-conditioning MCP knowledge."""

from radia_mcp.magnetic_materials.hysteresis_models_knowledge import (
    IDENTIFICATION_CONDITIONING,
    get_hysteresis_models_knowledge,
)


def test_identification_conditioning_dispatch_and_aliases():
    expected = get_hysteresis_models_knowledge("identification_conditioning")

    assert expected == IDENTIFICATION_CONDITIONING
    assert len(expected) > 5_000
    for alias in ("conditioning", "basis_sizing", "effective_dof", "linear_in_theta"):
        assert get_hysteresis_models_knowledge(alias) == expected


def test_conditioning_guidance_keeps_uncertainty_and_validation_honest():
    text = get_hysteresis_models_knowledge("identification_conditioning")

    assert "TRUNCATED Gaussian" in text
    assert "APPROXIMATION" in text
    assert "does NOT bound or" in text
    assert "not a universal setting" in text


def test_all_and_unknown_topic_advertise_conditioning():
    assert IDENTIFICATION_CONDITIONING in get_hysteresis_models_knowledge("all")
    assert "identification_conditioning" in get_hysteresis_models_knowledge("unknown")
