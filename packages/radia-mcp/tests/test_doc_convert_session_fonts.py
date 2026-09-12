"""Mock-only contracts: never call native font APIs or alter a live session."""
import ctypes
from types import SimpleNamespace

import pytest

from radia_mcp.doc_convert.plans import T17_session_fonts as t


class Call:
    def __init__(self, implementation):
        self.implementation = implementation
        self.calls = []

    def __call__(self, *args):
        self.calls.append(args)
        return self.implementation(*args)


@pytest.fixture
def native(monkeypatch):
    dc, font, previous = 0x123456789, 0x223456789, 0x323456789
    def name(_dc, _size, buffer):
        buffer.value = "Arial"
        return 5
    gdi = SimpleNamespace(
        CreateFontW=Call(lambda *a: font),
        SelectObject=Call(lambda *a: previous),
        DeleteObject=Call(lambda *a: 1),
        GetTextFaceW=Call(name),
        GetFontData=Call(lambda *a: 10000),
        GetGlyphOutlineW=Call(lambda *a: 200),
    )
    user = SimpleNamespace(GetDC=Call(lambda *a: dc), ReleaseDC=Call(lambda *a: 1))
    monkeypatch.setattr(t, "_gdi", lambda: t._bind(gdi, user))
    # The mock surfaces deliberately have no registration/removal/notification APIs.
    return gdi, user, dc, font, previous


def test_pointer_width_and_owned_resource_cleanup(native):
    gdi, user, dc, font, previous = native
    result = t.probe_face("Arial")
    assert result["ok"]
    assert ctypes.sizeof(user.GetDC.restype) == ctypes.sizeof(ctypes.c_void_p)
    assert ctypes.sizeof(gdi.CreateFontW.restype) == ctypes.sizeof(ctypes.c_void_p)
    assert ctypes.sizeof(gdi.SelectObject.restype) == ctypes.sizeof(ctypes.c_void_p)
    assert len(gdi.CreateFontW.argtypes) == 14
    assert gdi.SelectObject.calls == [(dc, font), (dc, previous)]
    assert gdi.DeleteObject.calls == [(font,)]
    assert user.ReleaseDC.calls == [(None, dc)]


@pytest.mark.parametrize("error", [-1, 0xFFFFFFFF])
@pytest.mark.parametrize("api", ["GetFontData", "GetGlyphOutlineW"])
def test_query_error_is_not_healthy_or_proof_of_font_host_damage(native, api, error):
    gdi, user, *_ = native
    getattr(gdi, api).implementation = lambda *a: error
    result = t.probe_face("Arial")
    assert result["status"] == "query_failed"
    assert not result["ok"]
    assert "broken" not in result
    assert gdi.DeleteObject.calls and user.ReleaseDC.calls


@pytest.mark.parametrize("stage", ["GetDC", "CreateFontW", "SelectObject", "GetTextFaceW"])
def test_partial_acquisition_failure_cleans_only_acquired_resources(native, stage):
    gdi, user, dc, font, _ = native
    getattr(user if stage == "GetDC" else gdi, stage).implementation = lambda *a: 0
    result = t.probe_face("Arial")
    assert not result["ok"] and result["status"] == "unverified"
    assert bool(user.ReleaseDC.calls) == (stage != "GetDC")
    assert bool(gdi.DeleteObject.calls) == (stage not in ("GetDC", "CreateFontW"))
    assert all(args == (font,) for args in gdi.DeleteObject.calls)


def test_query_exception_still_restores_and_releases(native):
    gdi, user, dc, font, previous = native
    def failure(*args):
        raise OSError("mock query failure")
    gdi.GetFontData.implementation = failure
    result = t.probe_face("Arial")
    assert not result["ok"]
    assert gdi.SelectObject.calls[-1] == (dc, previous)
    assert gdi.DeleteObject.calls == [(font,)]
    assert user.ReleaseDC.calls == [(None, dc)]


def test_cleanup_failure_is_not_a_pass(native):
    gdi, user, *_ = native
    gdi.DeleteObject.implementation = lambda *a: 0
    result = t.probe_face("Arial")
    assert not result["ok"] and result["cleanup_errors"]
    assert user.ReleaseDC.calls


def test_failed_restore_never_deletes_selected_font(native):
    gdi, user, *_ = native
    gdi.SelectObject.implementation = lambda *a: 1 if len(gdi.SelectObject.calls) == 1 else 0
    result = t.probe_face("Arial")
    assert not result["ok"] and result["cleanup_errors"]
    assert not gdi.DeleteObject.calls
    assert user.ReleaseDC.calls


def test_substitution_fails_aggregate_and_faces_are_deduplicated(native):
    result = t.doc_convert_session_font_check("Arial, Missing Face, Arial")
    assert not result["ok"]
    assert len(result["faces"]) == 2
    assert result["faces"][1]["status"] == "substituted"
    assert result["repair_performed"] is False


@pytest.mark.parametrize("faces", ["", ", ,", "x" * 32, "A\x00B", None,
                                   ",".join(str(i) for i in range(33))])
def test_invalid_request_never_touches_native_api(monkeypatch, faces):
    monkeypatch.setattr(t, "_gdi", lambda: pytest.fail("Unexpected native API call"))
    with pytest.raises(ValueError):
        t.doc_convert_session_font_check(faces)


def test_non_windows_reports_unverified_without_loading_libraries(monkeypatch):
    monkeypatch.setattr(t, "sys", SimpleNamespace(platform="linux"))
    result = t.doc_convert_session_font_check("Arial")
    assert not result["ok"]
    assert result["faces"][0]["error"] == "Windows only"


@pytest.fixture
def event_query(monkeypatch):
    monkeypatch.setattr(t, "sys", SimpleNamespace(platform="win32"))
    reply = SimpleNamespace(returncode=0, stdout=b"")
    calls = []
    def run(command, **kwargs):
        calls.append((command, kwargs))
        return reply
    monkeypatch.setattr(t.subprocess, "run", run)
    return reply, calls


def test_event_xml_not_localized_text_and_scope_is_explicit(event_query):
    reply, calls = event_query
    reply.stdout = ('<Events><Event xmlns="http://schemas.microsoft.com/win/2004/08/events/event">'
                    '<EventData><Data Name="AppName">fontdrvhost.exe</Data></EventData>'
                    '</Event></Events>').encode("utf-16")
    result = t.fontdrvhost_crashed_within(60)
    assert result["crashed"] is True and result["events"] == 1
    assert "not proof about this session" in result["scope"]
    assert "/f:xml" in calls[0][0] and "/uni:true" in calls[0][0]
    assert calls[0][1]["timeout"] == 15
    assert not calls[0][1].get("shell")


@pytest.mark.parametrize("output", [b"", "<Events/>".encode("utf-16")])
def test_successful_empty_query_is_no_crash(event_query, output):
    reply, _ = event_query
    reply.stdout = output
    assert t.fontdrvhost_crashed_within(60)["crashed"] is False


@pytest.mark.parametrize("failure", ["exit", "xml", "timeout"])
def test_event_failures_are_unknown(event_query, monkeypatch, failure):
    reply, _ = event_query
    if failure == "exit":
        reply.returncode = 5
    elif failure == "xml":
        reply.stdout = "<broken>".encode("utf-16")
    else:
        def timeout(*args, **kwargs):
            raise t.subprocess.TimeoutExpired("wevtutil", 15)
        monkeypatch.setattr(t.subprocess, "run", timeout)
    result = t.fontdrvhost_crashed_within(60)
    assert result["crashed"] is None and result["status"] == "unverified"


@pytest.mark.parametrize("seconds", [0, -1, True, 1.5, 86401])
def test_bad_window_rejected_without_process_launch(event_query, seconds):
    with pytest.raises(ValueError):
        t.fontdrvhost_crashed_within(seconds)
    assert not event_query[1]


@pytest.mark.parametrize("option", ["--repair", "--force", "--notify"])
def test_unsafe_legacy_cli_options_are_rejected(native, option):
    with pytest.raises(SystemExit) as exc:
        t.main([option])
    assert exc.value.code == 2
    assert not native[1].GetDC.calls


def test_cli_unknown_event_query_fails(native, event_query, capsys):
    event_query[0].returncode = 1
    assert t.main(["--faces", "Arial", "--crash-window", "60"]) == 1
    assert '"crashed": null' in capsys.readouterr().out


def test_registration_exposes_only_readonly_font_tool():
    from radia_mcp.doc_convert import tools
    assert tools.doc_convert_session_font_check is t.doc_convert_session_font_check
    assert not hasattr(tools, "doc_convert_session_font_repair")


@pytest.mark.parametrize("xml", ["garbage", "<Error/>", "<Event/>",
                                  "<Events><Error/></Events>", "<Events>bad</Events>"])
def test_unexpected_event_payload_is_unknown(event_query, xml):
    event_query[0].stdout = xml.encode("utf-16")
    assert t.fontdrvhost_crashed_within(60)["crashed"] is None


def test_truncated_event_coverage_is_unknown(event_query):
    event_query[0].stdout = ("<Event/>" * 1001).encode("utf-16")
    assert t.fontdrvhost_crashed_within(60)["crashed"] is None


def test_mcp_readonly_annotations():
    import asyncio
    from mcp.server.fastmcp import FastMCP
    from radia_mcp.doc_convert import register
    server = FastMCP("mock-font-contract")
    register(server)
    tool = next(tool for tool in asyncio.run(server.list_tools())
                if tool.name == "doc_convert_session_font_check")
    assert tool.annotations.readOnlyHint is True
    assert tool.annotations.destructiveHint is False
