"""Tests for Wave-3 cubit MCP features: session-mode triad, --setup CLI
mode, session-status enrichment, and the version-robust lint golden
((rule, line) pairs, not message text)."""

import json
from pathlib import Path

import pytest

from radia_mcp.cubit import server as cubit_server
from radia_mcp.cubit import session as cubit_session
from radia_mcp.cubit.server import _lint_file, cubit_session_status


# ---------------------------------------------------------------------------
# Session-mode triad
# ---------------------------------------------------------------------------

def _bare_session(tmp_path):
    sess = cubit_session.CubitSession.__new__(cubit_session.CubitSession)
    sess._bin_dir = tmp_path
    sess._mode = "gui"
    sess._proc = None
    sess._next_id = 1
    sess._ready_info = None
    sess._drop_dir = None
    sess._outbox = None
    sess._owned = False
    sess._last_license_warmup = {}
    sess._command_history = []
    sess._command_history_max = 10
    return sess


def test_invalid_session_mode_fails_loud(monkeypatch, tmp_path):
    monkeypatch.setenv("RADIA_CUBIT_SESSION_MODE", "sometimes")
    sess = _bare_session(tmp_path)
    with pytest.raises(cubit_session.CubitSessionError) as ei:
        sess._start_gui_bootstrap()
    assert "auto | new | existing" in str(ei.value)


def test_existing_mode_without_daemon_fails_loud(monkeypatch, tmp_path):
    monkeypatch.setenv("RADIA_CUBIT_SESSION_MODE", "existing")
    monkeypatch.setattr(cubit_session, "_user_daemon_dir",
                        lambda: tmp_path / "cubit-session")
    sess = _bare_session(tmp_path)
    with pytest.raises(cubit_session.CubitSessionError) as ei:
        sess._start_gui_bootstrap()
    assert "no live shared" in str(ei.value)


def test_new_mode_uses_private_drop_dir(monkeypatch, tmp_path):
    """Mode 'new' must ignore a live shared daemon and pick a private
    per-process drop dir (never clobber the shared one)."""
    import os

    shared = tmp_path / "cubit-session"
    shared.mkdir()
    (shared / "out").mkdir()
    (shared / "pid.lock").write_text(str(os.getpid()), encoding="utf-8")
    (shared / "ready").write_text(json.dumps(
        {"ready": True, "protocol_version": 2, "pid": os.getpid()}),
        encoding="utf-8")

    monkeypatch.setenv("RADIA_CUBIT_SESSION_MODE", "new")
    monkeypatch.setattr(cubit_session, "_user_daemon_dir", lambda: shared)
    sess = _bare_session(tmp_path)
    # Spawn path fails at launcher discovery (no Cubit under tmp_path) --
    # but by then the private drop dir must already exist and the shared
    # daemon's markers must be untouched.
    with pytest.raises((cubit_session.CubitSessionError, FileNotFoundError)):
        sess._start_gui_bootstrap()
    private = shared.parent / f"cubit-session-{os.getpid()}"
    assert private.is_dir()
    assert (shared / "pid.lock").is_file()          # shared daemon untouched
    assert (shared / "ready").is_file()


def test_attached_clients_use_distinct_request_filenames(tmp_path):
    first = _bare_session(tmp_path)
    second = _bare_session(tmp_path)

    first_stem = first._request_stem(1)
    second_stem = second._request_stem(1)

    assert first_stem != second_stem
    assert first_stem.endswith("-00000001")
    assert second_stem.endswith("-00000001")


# ---------------------------------------------------------------------------
# --setup CLI mode + status enrichment
# ---------------------------------------------------------------------------

def test_setup_mode_runs_and_reports(monkeypatch, capsys):
    # Avoid a real license warmup: stub it.
    import radia_mcp.cubit.license_warmup as lw
    monkeypatch.setattr(lw, "warmup_license",
                        lambda *a, **k: {"status": "skipped",
                                         "reason": "test stub"})
    rc = cubit_server._setup_mode()
    out = capsys.readouterr().out
    assert "Doctor report" in out
    assert rc in (0, 1)
    if rc == 1:
        assert "SETUP FOUND PROBLEMS" in out


def test_session_status_reports_mode_and_journal(monkeypatch):
    monkeypatch.setenv("RADIA_CUBIT_SESSION_MODE", "auto")
    out = json.loads(cubit_session_status())
    assert out["session_mode"] == "auto"
    # keys present whenever a singleton exists are optional here; the
    # baseline contract is bin_dir/alive/session_mode
    assert "alive" in out and "bin_dir" in out


# ---------------------------------------------------------------------------
# Version-robust lint golden: assert (rule, line) pairs, not messages
# ---------------------------------------------------------------------------

_FIXTURE = Path(__file__).resolve().parents[2].parent \
    / "tests" / "mcp_server" / "fixtures" / "bad_cubit_script.py"

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


@pytest.mark.skipif(not _FIXTURE.is_file(),
                    reason="repo-level fixture not present")
def test_lint_golden_rule_line_pairs():
    findings = _lint_file(str(_FIXTURE))
    got = {(f["rule"], f["line"]) for f in findings}
    missing = _EXPECTED_RULE_LINES - got
    unexpected = got - _EXPECTED_RULE_LINES
    assert not missing, f"lint findings disappeared: {sorted(missing)}"
    assert not unexpected, (
        f"new lint findings on the golden fixture: {sorted(unexpected)} "
        "-- if intentional, update _EXPECTED_RULE_LINES")


@pytest.mark.skipif(
    not (_FIXTURE.parent / "clean_cubit_script.py").is_file(),
    reason="repo-level fixture not present")
def test_lint_clean_fixture_stays_clean():
    findings = _lint_file(str(_FIXTURE.parent / "clean_cubit_script.py"))
    assert findings == []
