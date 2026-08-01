"""Fast MCP plumbing regressions for CoilBuilder pre-solve audits."""

import radia.coil_builder as coil_builder
from radia_mcp.electromagnet import server


def test_mcp_field_audit_forwards_contract(monkeypatch):
    calls = {}
    monkeypatch.setattr(server, "_load_coils_for_audit", lambda path: "coils")

    def fake_audit(coils, sample_points, **kwargs):
        calls.update(coils=coils, sample_points=sample_points, **kwargs)
        return {"passed": True, "schema": "radia.coil-field-audit/v1"}

    monkeypatch.setattr(
        coil_builder, "audit_coil_field_consistency", fake_audit
    )
    report = server.electromagnet_coil_field_audit(
        "coil.py",
        [[0.0, 0.0, 0.0]],
        n_arc=400,
        relative_tolerance=1.0e-3,
    )

    assert report["passed"]
    assert report["schema"] == "radia.coil-field-audit/v1"
    assert calls == {
        "coils": "coils",
        "sample_points": [[0.0, 0.0, 0.0]],
        "n_arc": 400,
        "arc_max_segment_length": None,
        "relative_tolerance": 1.0e-3,
        "absolute_tolerance_T": 1.0e-9,
        "closure_tolerance": 1.0e-9,
    }


def test_mcp_clearance_audit_forwards_contract_and_fails_loudly(
    tmp_path, monkeypatch
):
    calls = {}
    yoke_step = tmp_path / "yoke.step"
    yoke_step.write_text("contract fixture", encoding="ascii")
    monkeypatch.setattr(server, "_load_coils_for_audit", lambda path: "coils")

    def fake_audit(coils, yoke_path, **kwargs):
        calls.update(coils=coils, yoke_path=yoke_path, **kwargs)
        return {
            "passed": False,
            "no_overlap": False,
            "intersection_volume": 2.0e-6,
            "measured_clearance": 0.0,
            "minimum_clearance": kwargs["minimum_clearance"],
        }

    monkeypatch.setattr(
        coil_builder, "audit_coil_yoke_clearance", fake_audit
    )
    report = server.electromagnet_coil_yoke_clearance_audit(
        "coil.py", str(yoke_step), minimum_clearance=1.0e-3
    )
    assert not report["passed"]
    assert calls["coils"] == "coils"
    assert calls["yoke_path"] == yoke_step.resolve()
    assert calls["minimum_clearance"] == 1.0e-3
    assert calls["intersection_volume_tolerance"] == 1.0e-15

    try:
        server.electromagnet_coil_yoke_clearance_audit(
            "coil.py", str(yoke_step), fail_on_error=True
        )
    except RuntimeError as exc:
        assert "clearance audit failed" in str(exc)
    else:
        raise AssertionError("fail_on_error must reject an overlapping coil")
