"""Read-only, session-local font diagnostics; never repair shared font state.

Recovered from the shared WIP font investigation. Probing owns only a screen
DC and a temporary logical font. It does not register/unregister fonts, stop
processes, notify sessions, or establish the cause of a failed glyph query.
"""
from __future__ import annotations

import ctypes
from ctypes import wintypes as w
import json
import subprocess
import sys
import xml.etree.ElementTree as ET


def _bind(gdi, user):
    """Declare pointer-width-correct ABI before making any native call."""
    signatures = (
        (user.GetDC, [w.HWND], w.HDC),
        (user.ReleaseDC, [w.HWND, w.HDC], ctypes.c_int),
        (gdi.CreateFontW, [ctypes.c_int] * 5 + [w.DWORD] * 8 + [w.LPCWSTR], w.HFONT),
        (gdi.SelectObject, [w.HDC, w.HANDLE], w.HANDLE),
        (gdi.DeleteObject, [w.HANDLE], w.BOOL),
        (gdi.GetTextFaceW, [w.HDC, ctypes.c_int, w.LPWSTR], ctypes.c_int),
        (gdi.GetFontData, [w.HDC, w.DWORD, w.DWORD, ctypes.c_void_p, w.DWORD], w.DWORD),
        (gdi.GetGlyphOutlineW, [w.HDC, w.UINT, w.UINT, ctypes.c_void_p,
                               w.DWORD, ctypes.c_void_p, ctypes.c_void_p], w.DWORD),
    )
    for function, arguments, result in signatures:
        function.argtypes, function.restype = arguments, result
    return gdi, user


def _gdi():
    if sys.platform != "win32":
        raise OSError("Windows only")
    return _bind(ctypes.WinDLL("gdi32", use_last_error=True),
                 ctypes.WinDLL("user32", use_last_error=True))


def _invalid_handle(value):
    return value in (None, 0, -1, ctypes.c_void_p(-1).value)


def _is_gdi_error(value):
    return int(value) in (-1, 0xFFFFFFFF)


def probe_face(face: str) -> dict:
    """Inspect one realized face; failed queries do not prove font-host damage."""
    result = {"face": face, "ok": False, "status": "unverified"}
    hdc = font = previous = None
    cleanup_errors = []
    try:
        gdi, user = _gdi()
        hdc = user.GetDC(None)
        if _invalid_handle(hdc):
            raise OSError("GetDC failed")
        font = gdi.CreateFontW(-48, 0, 0, 0, 400, 0, 0, 0, 1, 0, 0, 0, 0, face)
        if _invalid_handle(font):
            raise OSError("CreateFontW failed")
        previous = gdi.SelectObject(hdc, font)
        if _invalid_handle(previous):
            raise OSError("SelectObject failed")
        name = ctypes.create_unicode_buffer(64)
        if not gdi.GetTextFaceW(hdc, len(name), name):
            raise OSError("GetTextFaceW failed")
        data = gdi.GetFontData(hdc, 0, 0, None, 0)
        metrics = ctypes.create_string_buffer(64)
        matrix = (ctypes.c_uint32 * 4)(0x10000, 0, 0, 0x10000)
        outline = gdi.GetGlyphOutlineW(hdc, ord("H"), 2, metrics, 0, None, matrix)
        substituted = name.value.casefold() != face.casefold()
        failed = _is_gdi_error(data) or _is_gdi_error(outline)
        status = "substituted" if substituted else "query_failed" if failed else "ok"
        result.update(status=status, ok=status == "ok", realised_as=name.value,
                      substituted=substituted,
                      font_data_bytes=None if _is_gdi_error(data) else int(data),
                      glyph_outline_bytes=None if _is_gdi_error(outline) else int(outline))
    except (OSError, ValueError, TypeError) as exc:
        result["error"] = str(exc)
    finally:
        # Never delete the caller's previous font; restore it before deletion.
        restored = True
        if not _invalid_handle(previous):
            try:
                restored = not _invalid_handle(gdi.SelectObject(hdc, previous))
            except Exception as exc:
                restored = False
                cleanup_errors.append(str(exc))
            if not restored:
                cleanup_errors.append("Could not restore previous font")
        if not _invalid_handle(font) and restored:
            try:
                if not gdi.DeleteObject(font):
                    cleanup_errors.append("DeleteObject failed")
            except Exception as exc:
                cleanup_errors.append(str(exc))
        if not _invalid_handle(hdc):
            try:
                if not user.ReleaseDC(None, hdc):
                    cleanup_errors.append("ReleaseDC failed")
            except Exception as exc:
                cleanup_errors.append(str(exc))
        if cleanup_errors:
            result.update(ok=False, status="unverified", cleanup_errors=cleanup_errors)
    return result


def fontdrvhost_crashed_within(seconds: int) -> dict:
    """Read bounded XML events; query errors remain unknown, never no-crash."""
    if isinstance(seconds, bool) or not isinstance(seconds, int) or not 1 <= seconds <= 86400:
        raise ValueError("seconds must be an integer in [1, 86400]")
    unknown = {"status": "unverified", "crashed": None, "events": None}
    if sys.platform != "win32":
        return {**unknown, "error": "Windows only"}
    query = ("*[System[Provider[@Name='Application Error'] and EventID=1000 and "
             f"TimeCreated[timediff(@SystemTime) <= {seconds * 1000}]]]")
    try:
        run = subprocess.run(["wevtutil", "qe", "Application", f"/q:{query}",
                              "/f:xml", "/uni:true", "/rd:true", "/c:1001"],
                             capture_output=True, timeout=15, check=False)
        if run.returncode:
            return {**unknown, "error": "wevtutil failed", "returncode": run.returncode}
        content = run.stdout.decode("utf-16").lstrip("\ufeff").strip() if run.stdout else ""
        if content.startswith("<?xml"):
            content = content.split("?>", 1)[1]
        root = ET.fromstring("<Root>" + content + "</Root>")
        if (root.text or "").strip() or any(
            node.tag.rsplit("}", 1)[-1] not in ("Events", "Event")
            for node in root
        ):
            raise ValueError("Unexpected event-log response")
        events = [event for event in root.iter() if event.tag.rsplit("}", 1)[-1] == "Event"]
        for container in (root, *[node for node in root.iter()
                                 if node.tag.rsplit("}", 1)[-1] == "Events"]):
            if (container.text or "").strip() or any(
                node.tag.rsplit("}", 1)[-1] not in ("Event", "Events")
                or (node.tail or "").strip() for node in container
            ):
                raise ValueError("Unexpected event container content")
        if len(events) >= 1001:
            return {**unknown, "error": "Event query limit reached; coverage incomplete"}
        matches = 0
        for event in events:
            data = [node for node in event.iter() if node.tag.rsplit("}", 1)[-1] == "Data"]
            if not data:
                raise ValueError("Event has no application data")
            app = next((node.text for node in data if node.get("Name") == "AppName"),
                       data[0].text if data else "")
            matches += (app or "").casefold() == "fontdrvhost.exe"
        return {"status": "observed", "crashed": bool(matches), "events": matches,
                "window_s": seconds, "scope": "Application log; not proof about this session"}
    except (OSError, subprocess.SubprocessError, ValueError, ET.ParseError) as exc:
        return {**unknown, "error": str(exc)}


def doc_convert_session_font_check(faces: str = "Times New Roman, Arial, Calibri") -> dict:
    """Read-only GDI checks for explicitly named faces in this process/session.

    A substituted face or query/cleanup error fails the check. No font repair,
    process restart or notification is performed. Results are not evidence that
    another process can export a document or that fontdrvhost caused a failure.
    """
    if not isinstance(faces, str):
        raise ValueError("faces must be a comma-separated string")
    names = list(dict.fromkeys(name.strip() for name in faces.split(",") if name.strip()))
    if not names or len(names) > 32 or any(len(name) > 31 or "\x00" in name for name in names):
        raise ValueError("Provide 1-32 face names of 1-31 characters without NUL")
    rows = [probe_face(name) for name in names]
    return {"ok": all(row["ok"] for row in rows), "faces": rows,
            "scope": "current process/session only", "repair_performed": False}


def main(argv=None) -> int:
    """Check-only CLI. The legacy --repair/--force/--notify options are rejected."""
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--faces", default="Times New Roman, Arial, Calibri")
    parser.add_argument("--crash-window", type=int)
    args = parser.parse_args(argv)
    result = doc_convert_session_font_check(args.faces)
    if args.crash_window is not None:
        result["crash"] = fontdrvhost_crashed_within(args.crash_window)
        if result["crash"]["status"] == "unverified":
            result["ok"] = False
    print(json.dumps(result, ensure_ascii=False))
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
