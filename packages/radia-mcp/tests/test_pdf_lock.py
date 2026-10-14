"""PDF lock release stops only dedicated viewers and never splices the file name."""

import subprocess

import pytest

from radia_mcp.common import pdf_lock


def test_browsers_are_never_stop_targets():
    for browser in ("msedge", "chrome", "firefox"):
        assert browser not in pdf_lock.PDF_VIEWER_PROCESSES
        assert browser not in pdf_lock._PROBE


@pytest.mark.parametrize("name", ["O'Neil.pdf", "a'; Stop-Process -Name explorer; '.pdf"])
def test_file_name_travels_through_environment_only(monkeypatch, tmp_path, name):
    pdf = tmp_path / name
    pdf.write_bytes(b"%PDF-1.4\n")
    calls = []

    def fake_run(argv, **kwargs):
        calls.append((argv, kwargs))
        if argv[0] == "powershell":
            return subprocess.CompletedProcess(argv, 0, stdout="4242\n", stderr="")
        return subprocess.CompletedProcess(argv, 1, stdout="", stderr="denied")

    monkeypatch.setattr(pdf_lock.os, "name", "nt")
    monkeypatch.setattr(pdf_lock, "_is_locked", lambda path: True)
    monkeypatch.setattr(pdf_lock.subprocess, "run", fake_run)

    stopped = pdf_lock.release_pdf_lock(pdf)

    probe_argv, probe_kwargs = calls[0]
    assert name not in " ".join(probe_argv)
    assert probe_kwargs["env"]["RADIA_MCP_PDF_LOCK_NAME"] == name
    assert calls[1][0] == ["taskkill", "/PID", "4242", "/F"]
    assert stopped == []  # a failed taskkill is not reported as stopped


def test_unlocked_file_starts_no_process(monkeypatch, tmp_path):
    pdf = tmp_path / "free.pdf"
    pdf.write_bytes(b"%PDF-1.4\n")
    monkeypatch.setattr(pdf_lock.os, "name", "nt")
    monkeypatch.setattr(pdf_lock, "_is_locked", lambda path: False)
    monkeypatch.setattr(pdf_lock.subprocess, "run",
                        lambda *a, **k: pytest.fail("no process may start"))
    assert pdf_lock.release_pdf_lock(pdf) == []
