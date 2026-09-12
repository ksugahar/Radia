"""Offline regression coverage; no meshes, native imports or solves."""

import copy
import importlib.util
import json
from pathlib import Path

import pytest

PATH = Path(__file__).resolve().parents[1] / "validation_test/feec/check_wedge_charge_pullback.py"
SPEC = importlib.util.spec_from_file_location("pullback_contract", PATH)
checker = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(checker)


@pytest.fixture
def lanes():
    native = {"lane": "ngsolve-native", "python": "test", "ngsolve": "test", "cases": {}}
    for name, elements in checker.CASES.items():
        case = {"elements": elements, "boundary_faces": 2, "ndof": 10}
        for family, count in (("volume", elements), ("boundary", 2)):
            rows = []
            for host in range(count):
                triangle = host % 2 == 0
                points = checker.VOLUME_POINTS if family == "volume" else (
                    checker.TRI_POINTS if triangle else checker.QUAD_POINTS)
                for point in points:
                    row = {"host": host, "point": list(point),
                           "physical_point": [.01, .02, .03], "reference_charge": 2.0}
                    if family == "boundary":
                        row["triangle"] = triangle
                    rows.append(row)
            case[family] = rows
        native["cases"][name] = case
    radia = copy.deepcopy(native)
    radia["lane"] = "radia-charge-pullback"
    for case in radia["cases"].values():
        case["charge_modes"] = 20
    return native, radia


def test_agreement_is_not_acceptance_or_root_cause(lanes):
    result = checker.compare(*lanes)
    assert result["status"] == "sampled_agreement"
    assert result["sample_checks_passed"]
    for key in ("acceptance_evidence", "campaign_level", "runtime_provenance_verified",
                "basis_rank_verified", "gram_root_cause_established"):
        assert result[key] is False


@pytest.mark.parametrize("field,value", [
    ("reference_charge", float("nan")), ("reference_charge", float("inf")),
    ("reference_charge", True), ("reference_charge", "2"),
    ("physical_point", [0, 0]), ("physical_point", [0, 0, float("nan")]),
    ("point", [.5, .5, .5]), ("host", True), ("host", -1),
])
@pytest.mark.parametrize("lane", [0, 1])
def test_bad_samples_fail_closed(lanes, field, value, lane):
    lanes[lane]["cases"]["two_structured_wedges"]["volume"][0][field] = value
    result = checker.compare(*lanes)
    assert result["status"] == "invalid_input"
    assert not result["sample_checks_passed"]


@pytest.mark.parametrize("damage", ["missing_case", "extra_case", "empty", "duplicate",
                                    "ndof", "version", "face_type", "overflow"])
def test_incomplete_or_inconsistent_inputs(lanes, damage):
    native, radia = lanes
    case = radia["cases"]["two_structured_wedges"]
    if damage == "missing_case":
        del radia["cases"]["explicit_full_four_wedges"]
    elif damage == "extra_case":
        radia["cases"]["unexpected"] = case
    elif damage == "empty":
        case["boundary"] = []
    elif damage == "duplicate":
        case["volume"][1] = copy.deepcopy(case["volume"][0])
    elif damage == "ndof":
        case["ndof"] = 11
    elif damage == "version":
        radia["ngsolve"] = "other"
    elif damage == "face_type":
        case["boundary"][0]["triangle"] = False
    elif damage == "overflow":
        case["volume"][0]["reference_charge"] = 1e308
        native["cases"]["two_structured_wedges"]["volume"][0]["reference_charge"] = -1e308
    assert checker.compare(*lanes)["status"] == "invalid_input"


def test_last_case_mismatch_is_not_hidden(lanes):
    lanes[1]["cases"]["explicit_full_four_wedges"]["boundary"][-1]["reference_charge"] += 1
    result = checker.compare(*lanes)
    assert result["status"] == "sampled_mismatch"
    assert not result["sample_checks_passed"]
    assert not result["gram_root_cause_established"]


@pytest.mark.parametrize("invalid", [False, True])
def test_cli_exit_json_and_input_preservation(lanes, tmp_path, invalid):
    paths = [tmp_path / name for name in ("native.json", "radia.json", "result.json")]
    for path, lane in zip(paths, lanes):
        path.write_text(json.dumps(lane), encoding="utf-8")
    if invalid:
        paths[1].write_text('{"cases": {}, "cases": {}}', encoding="utf-8")
    original = [p.read_bytes() for p in paths[:2]]
    code = checker.main(["--native-json", str(paths[0]), "--radia-json", str(paths[1]),
                         "--output", str(paths[2])])
    assert code == int(invalid)
    result = json.loads(paths[2].read_text(encoding="utf-8"))
    assert result["sample_checks_passed"] is (not invalid)
    assert not result["acceptance_evidence"]
    assert [p.read_bytes() for p in paths[:2]] == original


def test_cli_rejects_output_over_input(tmp_path):
    path = tmp_path / "original.json"
    path.write_text("{}", encoding="utf-8")
    with pytest.raises(SystemExit, match="2"):
        checker.main(["--native-json", str(path), "--radia-json", str(path),
                      "--output", str(path)])
    assert path.read_text(encoding="utf-8") == "{}"


def test_cli_numerical_mismatch_exits_nonzero(lanes, tmp_path):
    lanes[1]["cases"]["explicit_full_four_wedges"]["volume"][-1]["reference_charge"] += 1
    native, radia, output = (tmp_path / name for name in ("n.json", "r.json", "o.json"))
    native.write_text(json.dumps(lanes[0]), encoding="utf-8")
    radia.write_text(json.dumps(lanes[1]), encoding="utf-8")
    assert checker.main(["--native-json", str(native), "--radia-json", str(radia),
                         "--output", str(output)]) == 1
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["status"] == "sampled_mismatch"
    assert set(result["input_sha256"]) == {"native", "radia"}


def test_rows_may_be_reordered_without_losing_identity(lanes):
    for case in lanes[1]["cases"].values():
        case["volume"].reverse()
        case["boundary"].reverse()
    assert checker.compare(*lanes)["sample_checks_passed"]
