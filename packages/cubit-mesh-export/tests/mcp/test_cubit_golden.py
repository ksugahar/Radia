from pathlib import Path
from cubit_mesh_export.mcp.server import _lint_file


# ---------------------------------------------------------------------------
# Version-robust lint golden: assert (rule, line) pairs, not messages
# ---------------------------------------------------------------------------

_FIXTURE = Path(__file__).parent / "fixtures" / "bad_cubit_script.py"

# Locked 2026-08-05 against the committed fixture.  Assert only
# (rule, line) membership so rule WORDING can change freely; a rule
# deletion or a fixture edit must consciously update this set.
_EXPECTED_RULE_LINES = {
    ("deleted-api-usage", 11),
    ("deleted-api-usage", 46),
    ("missing-mesh-command", 50),
    ("geometry-block-2nd-order", 32),
    ("element-type-before-add", 35),
    ("wrong-connectivity-2nd-order", 39),
    ("nodeset-sideset-usage", 42),
    ("missing-boundary-block", 1),
    ("hardcoded-absolute-path", 23),
    ("wrong-file-extension", 50),
    ("missing-block-names", 1),
}


def test_lint_golden_rule_line_pairs():
    findings = _lint_file(str(_FIXTURE))
    got = {(f["rule"], f["line"]) for f in findings}
    missing = _EXPECTED_RULE_LINES - got
    unexpected = got - _EXPECTED_RULE_LINES
    assert not missing, f"lint findings disappeared: {sorted(missing)}"
    assert not unexpected, (
        f"new lint findings on the golden fixture: {sorted(unexpected)} "
        "-- if intentional, update _EXPECTED_RULE_LINES")


def test_lint_clean_fixture_stays_clean():
    findings = _lint_file(str(_FIXTURE.parent / "clean_cubit_script.py"))
    assert findings == []
