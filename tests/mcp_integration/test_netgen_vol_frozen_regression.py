"""Lock calculation outputs to frozen prior records; analytic acceptance is separate."""
import pytest
import json
import hashlib
from contextlib import contextmanager
from types import SimpleNamespace
import subprocess
from validation_test.radia_mcp.netgen_vol_examples import CASES, example_module
from validation_test.radia_mcp.netgen_vol_regression import assert_summary_regression

@pytest.mark.parametrize("module_name", CASES, ids=CASES)
def test_frozen_summary(module_name, tmp_path, monkeypatch):
    with example_module(module_name) as module:
        def forbidden_acceptance(*args, **kwargs):
            raise AssertionError("Analytic acceptance belongs to validation_test")
        monkeypatch.setattr(module, "validate_summary", forbidden_acceptance)
        summary = module.build_summary()
    assert_summary_regression(module_name, summary)


def test_p1_cli_does_not_bypass_validation(tmp_path, monkeypatch):
    """The CLI must dispatch validation before writing, independent of its criteria."""
    import sys
    output = tmp_path / "p1.json"
    with example_module("validation_p1_surface_triangle_element_summary") as module:
        def refused(_summary):
            raise RuntimeError("validation dispatch reached")
        monkeypatch.setattr(module, "validate_summary", refused)
        monkeypatch.setattr(sys, "argv", [module.__file__, "--out", str(output)])
        with pytest.raises(RuntimeError, match="validation dispatch reached"):
            module.main()
    assert not output.exists()


def test_summary_adapter_dispatches_zero_tolerance_and_preserves_metadata(monkeypatch):
    from validation_test.radia_mcp import netgen_vol_regression as adapter
    reference = {"nested": {"x": [1.0, -0.0]}, "flag": True}
    provenance = {"reference_source_sha": "a" * 40, "fixture_sha256": "b" * 64, "git_blob": "c" * 40}
    monkeypatch.setattr(adapter, "approved_case", lambda _: (None, reference, provenance))
    captured = []
    original = adapter.assert_regression
    def inspect(actual, baseline, *, conditions):
        captured.append(baseline)
        return original(actual, baseline, conditions=conditions)
    monkeypatch.setattr(adapter, "assert_regression", inspect)
    adapter.assert_summary_regression("synthetic", {"nested": {"x": (1.0, 0.0)}, "flag": True})
    assert all(q["atol"] == q["rtol"] == 0 for q in captured[0]["quantities"].values())
    assert captured[0]["conditions"]["provenance"] == provenance
    with pytest.raises(AssertionError, match="actual=.*expected="):
        adapter.assert_summary_regression("synthetic", {"nested": {"x": [1.0 + 1e-12, 0.0]}, "flag": True})
    with pytest.raises(AssertionError, match="conditions"):
        adapter.assert_summary_regression("synthetic", {"nested": {"x": [1.0, 0.0]}, "flag": False})


def test_generator_delegates_to_common_proposal_without_analytic_acceptance(tmp_path, monkeypatch):
    import tools.regression_lock as lock
    from validation_test.radia_mcp import regenerate_netgen_vol_regression as generator
    source = tmp_path / "source"
    source.mkdir()
    subprocess.run(["git", "init", "-q", str(source)], check=True)
    (source / "marker").write_text("synthetic", encoding="utf-8")
    subprocess.run(["git", "add", "marker"], cwd=source, check=True)
    subprocess.run(["git", "-c", "user.name=Synthetic", "-c", "user.email=synthetic@example.invalid",
                    "commit", "-qm", "Synthetic"], cwd=source, check=True)
    monkeypatch.setattr(lock, "SOURCE_ROOT", source)
    monkeypatch.setattr(generator, "ROOT", source)
    monkeypatch.setattr(generator, "CASES", ["synthetic"])
    frozen = tmp_path / "frozen.json"
    raw = b'{"numeric": -0.0, "label": "fixture"}\r\n'
    frozen.write_bytes(raw)
    provenance = {"fixture_sha256": hashlib.sha256(raw).hexdigest(), "git_blob": "a" * 40}
    monkeypatch.setattr(generator, "approved_case", lambda _: (frozen, json.loads(raw), provenance))
    @contextmanager
    def example(_):
        def forbidden(_summary):
            raise AssertionError("Analytic acceptance must not run in regeneration")
        yield SimpleNamespace(build_summary=lambda: {"numeric": 0.0, "label": "fixture"}, validate_summary=forbidden)
    monkeypatch.setattr(generator, "example_module", example)
    calls = []
    original = generator.propose_json_regeneration
    def delegate(*args, **kwargs):
        calls.append(kwargs)
        return original(*args, **kwargs)
    monkeypatch.setattr(generator, "propose_json_regeneration", delegate)
    out = tmp_path / "proposal"
    generator.main(["--out", str(out), "--reason", "Synthetic generator contract"])
    assert len(calls) == 1 and calls[0]["conditions"]["reference"] == provenance
    assert (out / "synthetic/baseline.json").read_bytes() == frozen.read_bytes() == raw
    report = json.loads((out / "synthetic/proposal.json").read_text(encoding="utf-8"))
    assert report["generation"]["dirty"] is False and report["diff"] == []
    assert (out / "index.json").is_file()
