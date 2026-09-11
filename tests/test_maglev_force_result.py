from __future__ import annotations

from radia_mcp.force.gates import force_action_reaction_gate

from radia.maglev.ecb import lorentz


def test_foster_lorentz_result_records_action_reaction_and_phasor_contract(monkeypatch):
    """Plumbing only: the stub replaces the solve, so no physics is checked here.

    The physical correctness of compute_lorentz_force_via_foster is held by
    validation_test/maglev/test_ecb_foster_lorentz_reference_evidence.py.  The
    action-reaction residual below is zero by construction (the source force is
    the negated conductor force), so it confirms the bookkeeping, not the force.
    """
    monkeypatch.setattr(
        lorentz,
        "compute_lorentz_force_via_foster",
        lambda *args, **kwargs: (0.125, -0.5, -3.25),
    )

    result = lorentz.compute_lorentz_force_result_via_foster(
        None,
        None,
        None,
        None,
        1.0,
        1.0,
        1j,
        1.0,
        1.0,
        0.0,
        frame="plate",
    )

    assert result["schema"] == "radia.maglev-force-pair/v1"
    assert result["conductor"]["force_N"] == [0.125, -0.5, -3.25]
    assert result["source_pm"]["force_N"] == [-0.125, 0.5, 3.25]
    assert result["action_reaction_residual_N"] == [0.0, 0.0, 0.0]
    for side in ("conductor", "source_pm"):
        assert result[side]["schema"] == "radia.force-result/v1"
        assert result[side]["field_convention"] == "time_average_phasor"
        assert result[side]["phasor_amplitude"] == "peak"
        assert result[side]["frame"] == "plate"

    gate = force_action_reaction_gate(
        result["conductor"]["force_N"],
        result["source_pm"]["force_N"],
    )
    assert gate["status"] == "ok"
