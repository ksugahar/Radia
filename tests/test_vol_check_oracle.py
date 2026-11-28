"""check-vol (radia.vol_check) keeps its adopted verdicts on fixed good and bad meshes.

Guards: radia.vol_check, the entry gate every solver-bound .vol passes before
initialization. For each case in tests/fixtures/vol_check/cases.py (committed
fixtures, label-contract violations and topology violations) the normalized
report must equal the adopted baseline (tests/fixtures/vol_check/expected.json,
recorded from radia.vol_check at its source_sha after it matched the former
standalone checker on every case); input errors must match by type. The schema
is radia.vol-check.v1.
Assumes: NGSolve is importable. Measured floats are compared at 1e-9 relative
so a different NGSolve build does not change the verdict; every flag, count,
name and issue list must match exactly.
"""
import importlib.util
import json
import math
from pathlib import Path

import pytest

pytest.importorskip("ngsolve")

HERE = Path(__file__).resolve().parent / "fixtures" / "vol_check"
_spec = importlib.util.spec_from_file_location("vol_check_cases", HERE / "cases.py")
cases = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cases)
EXPECTED = json.loads((HERE / "expected.json").read_text(encoding="utf-8"))
CASES = cases.cases()


def _diff(actual, expected, path="$"):
    if isinstance(expected, float) or isinstance(actual, float):
        if not (isinstance(actual, (int, float)) and isinstance(expected, (int, float))):
            return [f"{path}: {actual!r} != {expected!r}"]
        if not math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-12):
            return [f"{path}: {actual!r} != {expected!r}"]
        return []
    if isinstance(expected, dict) and isinstance(actual, dict):
        out = [f"{path}: keys {sorted(set(actual) ^ set(expected))} differ"] if set(actual) != set(expected) else []
        for key in sorted(set(actual) & set(expected)):
            out += _diff(actual[key], expected[key], f"{path}.{key}")
        return out
    if isinstance(expected, list) and isinstance(actual, list):
        if len(actual) != len(expected):
            return [f"{path}: length {len(actual)} != {len(expected)}"]
        return [d for i, (a, e) in enumerate(zip(actual, expected)) for d in _diff(a, e, f"{path}[{i}]")]
    return [] if actual == expected else [f"{path}: {actual!r} != {expected!r}"]


def test_case_set_matches_the_baseline():
    assert EXPECTED["schema"] == "radia.vol-check-baseline.v1", EXPECTED["schema"]
    assert sorted(CASES) == sorted(EXPECTED["cases"]), (sorted(CASES), sorted(EXPECTED["cases"]))


@pytest.mark.parametrize("name", sorted(CASES))
def test_vol_check_matches_the_baseline(name, tmp_path):
    from radia.vol_check import check_consistency

    mesh, contract = cases.write_case(tmp_path, name, CASES[name])
    kwargs = dict(strict_labels=True, contract=str(contract)) if contract else {}
    expected = EXPECTED["cases"][name]
    if "error" in expected:
        with pytest.raises(Exception) as caught:
            check_consistency(str(mesh), **kwargs)
        assert type(caught.value).__name__ == expected["error"], (type(caught.value).__name__, expected)
        return
    report = json.loads(json.dumps(check_consistency(str(mesh), **kwargs)))
    assert report.pop("schema") == "radia.vol-check.v1"
    actual = cases.normalize(report)
    differences = _diff(actual, expected["report"])
    assert not differences, "\n".join(differences[:20])
