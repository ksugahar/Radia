from pathlib import Path
from radia_mcp.cubit.server import _lint_file


# ---------------------------------------------------------------------------
# Version-robust lint golden: assert (rule, line) pairs, not messages
# ---------------------------------------------------------------------------

_FIXTURE = Path(__file__).parent / "fixtures" / "bad_cubit_script.py"

# Locked 2026-11-27 against the committed fixture.  Assert only
# (rule, line) membership so rule WORDING can change freely; a rule
# deletion or a fixture edit must consciously update this set.
_EXPECTED_RULE_LINES = {
    ("missing-mesh-command", 36),
    ("geometry-block-2nd-order", 25),
    ("element-type-before-add", 28),
    ("wrong-connectivity-2nd-order", 32),
    ("missing-boundary-block", 1),
    ("hardcoded-absolute-path", 19),
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
