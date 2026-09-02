from __future__ import annotations

from radia_mcp.mathematica import server


def test_selftest_does_not_start_wolfram_kernel(monkeypatch, capsys):
    def unexpected_runtime_probe():
        raise AssertionError("selftest must not consume a Mathematica license")

    monkeypatch.setattr(server.sys, "argv", ["mcp-server-mathematica", "--selftest"])
    monkeypatch.setattr(server._tools, "mathematica_status", unexpected_runtime_probe)
    monkeypatch.setattr(server._tools, "_find_wolframscript", lambda: None)

    server.main()

    output = capsys.readouterr().out
    assert "wolframscript: NOT FOUND" in output
    assert "PASSED" in output
