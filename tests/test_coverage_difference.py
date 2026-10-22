"""Hand-specified set arithmetic and failure contracts; no solver executions."""
import copy
import json

import pytest

from tools.coverage_difference import compare_exports, load_export, main


@pytest.fixture
def pair():
    heavy = {
        "schema": "radia.coverage-evidence.v1", "role": "heavy",
        "identity": {
            "source_sha": "a" * 40, "artifacts": {"radia.whl": "b" * 64}, "backend": "cpu",
            "route": {"entry_point": "radia.synthetic.solve", "formulation": "scalar", "element": "H1",
                      "order": 2, "nonlinear": True, "boundary": "open", "constitutive": "declared-model"},
            "instrumentation": {"tool": "synthetic", "version": "1", "branch": True,
                                "exclusions": [], "subprocess": {"required": False, "sentinel_verified": False}},
            "production_scope": ["src/radia/synthetic.py"]},
        "run": {"command": ["python", "synthetic-heavy.py"], "conditions": {"mesh": "large", "threads": 1},
                "runtime_s": 9999, "gate": {"status": "passed", "evidence_sha256": "c" * 64}},
        "files": {"src/radia/synthetic.py": {"sha256": "d" * 64, "lines": [1, 2, 3], "arcs": [[-1, 1], [1, 2], [2, 3]]}},
        "native": {}}
    light = copy.deepcopy(heavy)
    light["role"] = "light"
    light["run"].update(command=["python", "synthetic-light.py"], conditions={"mesh": "small", "threads": 1}, runtime_s=31)
    light["files"]["src/radia/synthetic.py"].update(lines=[1, 2, 4], arcs=[[-1, 1], [1, 2], [2, 4]])
    return heavy, light


def test_reports_line_and_branch_differences_without_automatic_downgrade(pair):
    heavy, light = pair
    before = copy.deepcopy(pair)
    report = compare_exports(heavy, light)
    delta = report["files"]["src/radia/synthetic.py"]
    assert delta["lines"] == {"both": [1, 2], "heavy_only": [3], "light_only": [4]}
    assert delta["arcs"] == {"both": [(-1, 1), (1, 2)], "heavy_only": [(2, 3)], "light_only": [(2, 4)]}
    assert report["automatic_downgrade"] is False
    assert report["recommendation"] == "owner-review-required"
    assert report["runs"]["light"]["runtime_s"] == 31  # metadata, not a 30s assertion
    assert pair == before


@pytest.mark.parametrize("key,value", [("source_sha", "e" * 40), ("artifacts", {"radia.whl": "e" * 64}),
                                       ("backend", "gpu"), ("production_scope", ["src/radia/other.py"])])
def test_pair_identity_mismatch_is_loud(pair, key, value):
    pair[1]["identity"][key] = value
    if key == "production_scope":
        pair[1]["files"] = {}
    with pytest.raises(ValueError, match=f"pair identity.{key}.*actual=.*expected="):
        compare_exports(*pair)


@pytest.mark.parametrize("mutation", ["branch", "sentinel", "failed", "partial", "timed-out", "missing-file"])
def test_missing_capture_and_failed_numerical_verdict_stay_unknown(pair, mutation):
    light = pair[1]
    if mutation == "branch":
        light["identity"]["instrumentation"]["branch"] = False
    elif mutation == "sentinel":
        light["identity"]["instrumentation"]["subprocess"]["required"] = True
    elif mutation == "missing-file":
        light["files"] = {}
    else:
        light["run"]["gate"]["status"] = mutation
    report = compare_exports(*pair)
    assert report["unknowns"]
    assert report["recommendation"] == "insufficient-evidence"
    assert report["runs"]["light"]["gate"] == light["run"]["gate"]


def test_equal_lines_do_not_hide_a_missing_branch(pair):
    pair[1]["files"]["src/radia/synthetic.py"]["lines"] = [1, 2, 3]
    delta = compare_exports(*pair)["files"]["src/radia/synthetic.py"]
    assert delta["lines"]["heavy_only"] == []
    assert delta["arcs"]["heavy_only"] == [(2, 3)]


def test_empty_capture_and_missing_gpu_hits_do_not_become_qualification(pair):
    for data in pair:
        data["identity"]["backend"] = "gpu"
        data["files"]["src/radia/synthetic.py"].update(lines=[], arcs=[])
    report = compare_exports(*pair)
    assert report["recommendation"] == "insufficient-evidence"
    assert any("no executed production lines" in item for item in report["unknowns"])
    assert any("GPU entry-point evidence absent" in item for item in report["unknowns"])


@pytest.mark.parametrize("mutation", ["exclusions", "source-hash", "route", "bad-path", "duplicate-lines", "bool-line"])
def test_uncomparable_or_malformed_evidence_is_rejected(pair, mutation):
    light = pair[1]
    if mutation == "exclusions":
        light["identity"]["instrumentation"]["exclusions"] = ["an excluded branch"]
    elif mutation == "source-hash":
        light["files"]["src/radia/synthetic.py"]["sha256"] = "e" * 64
    elif mutation == "route":
        light["identity"]["route"]["nonlinear"] = False
    elif mutation == "bad-path":
        light["identity"]["production_scope"] = ["../outside.py"]
    else:
        light["files"]["src/radia/synthetic.py"]["lines"] = [1, 1] if mutation == "duplicate-lines" else [True]
    with pytest.raises(ValueError, match="actual=.*expected="):
        compare_exports(*pair)


def test_native_hits_never_claim_internal_coverage(pair):
    heavy, light = pair
    heavy["native"]["synthetic_entry"] = {"calls": 30, "artifact": "radia.whl", "configuration": {"order": 2}}
    light["native"]["synthetic_entry"] = {"calls": 2, "artifact": "radia.whl", "configuration": {"order": 2}}
    report = compare_exports(*pair)
    assert report["route_matrix"] == [{"entry_point": "synthetic_entry", "state": "covered-both", "internal_coverage": "UNKNOWN"}]
    light["native"]["synthetic_entry"]["configuration"]["order"] = 3
    report = compare_exports(*pair)
    assert report["route_matrix"][0]["state"] == "unknown" and report["unknowns"]


def test_cli_preserves_inputs_and_refuses_output_overwrite(pair, tmp_path):
    paths = [tmp_path / "heavy.json", tmp_path / "light.json"]
    for path, data in zip(paths, pair):
        path.write_text(json.dumps(data), encoding="utf-8")
    before = [p.read_bytes() for p in paths]
    out = tmp_path / "report.json"
    argv = ["--heavy", str(paths[0]), "--light", str(paths[1]), "--out", str(out)]
    assert main(argv) == 0
    assert load_export(out)["automatic_downgrade"] is False
    with pytest.raises(FileExistsError):
        main(argv)
    assert [p.read_bytes() for p in paths] == before


@pytest.mark.parametrize("text", ['{"a":1,"a":2}', '{"a":NaN}'])
def test_bad_json_never_becomes_evidence(tmp_path, text):
    path = tmp_path / "input.json"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError):
        load_export(path)
