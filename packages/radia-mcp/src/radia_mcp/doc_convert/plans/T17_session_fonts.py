"""Tier 1 — doc_convert_session_font_check / doc_convert_session_font_repair.

A Windows session can lose ONE typeface while every other one keeps
working: PowerPoint drops the text set in it (PDF export silently, PNG
export with "an error occurred while saving"), Word and LINE draw blanks,
and at the GDI level ``GetFontData`` / ``GetGlyphOutline`` return
GDI_ERROR for that face only.  Seen on LAB with Calibri (2026-08-31),
Arial (09-01) and Times New Roman (09-04).

What breaks is not the font file and not the FontCache service; it is
the session's win32k entry for that face, poisoned when the session's
user-mode font driver host (``fontdrvhost.exe``) crashed while the face
was being realised.  The file is fine, other sessions are fine, and a
reboot is NOT needed: dropping the session's references to the face's
files (``RemoveFontResourceW`` until it returns FALSE) and adding them
back (``AddFontResourceW`` + ``WM_FONTCHANGE``) makes win32k realise the
face again through the now-healthy driver host.  Verified 2026-09-04 on
Times New Roman: GetFontData -1 -> 1,166,948 bytes, and a fresh
PowerPoint exported the deck with its body text again.

The check probes each face through GDI in THIS process; the repair only
issues public font-resource calls for the session, changes nothing on
disk or in the registry, and is idempotent on a healthy face (the face
is unavailable for the few milliseconds between remove and add, so do
not run it while a document is being exported).
"""
from __future__ import annotations

import os
import pathlib
import sys

_GDI_ERROR = 0xFFFFFFFF
_WM_FONTCHANGE = 0x001D
_HWND_BROADCAST = 0xFFFF
_SMTO_ABORTIFHUNG = 0x0002

# the file names Windows ships for the faces that have broken so far
_KNOWN_FILES = {
    "times new roman": ("times.ttf", "timesbd.ttf", "timesi.ttf", "timesbi.ttf"),
    "arial": ("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"),
    "calibri": ("calibri.ttf", "calibrib.ttf", "calibrii.ttf", "calibriz.ttf", "calibril.ttf", "calibrili.ttf"),
    "cambria": ("cambria.ttc", "cambriab.ttf", "cambriai.ttf", "cambriaz.ttf"),
    "cambria math": ("cambria.ttc",),
    "georgia": ("georgia.ttf", "georgiab.ttf", "georgiai.ttf", "georgiaz.ttf"),
    "segoe ui": ("segoeui.ttf", "segoeuib.ttf", "segoeuii.ttf", "segoeuiz.ttf", "segoeuil.ttf", "segoeuisl.ttf"),
}


def _gdi():
    if sys.platform != "win32":
        return None, None
    import ctypes
    return ctypes.windll.gdi32, ctypes.windll.user32


def probe_face(face: str) -> dict:
    """GDI view of one face in this session: realised name, font data size, glyph outline size."""
    import ctypes
    gdi, user = _gdi()
    if gdi is None:
        return {"face": face, "error": "Windows only"}
    hdc = user.GetDC(0)
    hf = gdi.CreateFontW(-48, 0, 0, 0, 400, 0, 0, 0, 1, 0, 0, 0, 0, face)
    old = gdi.SelectObject(hdc, hf)
    buf = ctypes.create_unicode_buffer(64)
    gdi.GetTextFaceW(hdc, 64, buf)
    data = gdi.GetFontData(hdc, 0, 0, None, 0)
    gm = ctypes.create_string_buffer(64)
    mat = (ctypes.c_uint32 * 4)(0x10000, 0, 0, 0x10000)
    outline = gdi.GetGlyphOutlineW(hdc, ord("H"), 2, gm, 0, None, mat)   # GGO_NATIVE
    gdi.SelectObject(hdc, old)
    gdi.DeleteObject(hf)
    user.ReleaseDC(0, hdc)
    realised = buf.value
    substituted = realised.casefold() != face.casefold()
    broken = (not substituted) and (data == _GDI_ERROR or outline == _GDI_ERROR)
    return {
        "face": face,
        "realised_as": realised,
        "substituted": substituted,
        "font_data_bytes": None if data == _GDI_ERROR else int(data),
        "glyph_outline_bytes": None if outline == _GDI_ERROR else int(outline),
        "ok": not broken and not substituted,
        "broken": broken,
    }


def _font_files(face: str) -> list[pathlib.Path]:
    fonts_dir = pathlib.Path(os.environ.get("WINDIR", r"C:\Windows")) / "Fonts"
    names = _KNOWN_FILES.get(face.casefold())
    if names is None:
        # registry: every file whose registered name starts with the face
        try:
            import winreg
            key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,
                                 r"SOFTWARE\Microsoft\Windows NT\CurrentVersion\Fonts")
            names = []
            i = 0
            while True:
                try:
                    name, value, _t = winreg.EnumValue(key, i)
                except OSError:
                    break
                i += 1
                if name.casefold().startswith(face.casefold()) and isinstance(value, str):
                    names.append(value)
        except OSError:
            names = []
    files = []
    for n in names or ():
        p = pathlib.Path(n)
        if not p.is_absolute():
            p = fonts_dir / n
        if p.exists():
            files.append(p)
    return files


def doc_convert_session_font_check(faces: str = "Times New Roman, Arial, Calibri, Cambria Math, Georgia, Segoe UI") -> dict:
    """このセッションで各書体が GDI から本当に描けるかを確かめる。

    fontdrvhost.exe がクラッシュすると、そのとき読み込み中だった書体だけがセッションの
    フォント表で壊れ、PowerPoint はその書体の文字を PDF では無音で落とし、PNG 書き出しは
    「ファイルの保存中にエラー」で失敗する。列挙は通るので普通の確認では分からない。
    このツールは書体ごとに GetFontData / GetGlyphOutline を呼び、GDI_ERROR を返す書体を
    ``broken`` として報告する（別の書体に置き換わった場合は ``substituted``）。

    Args:
        faces: 点検する書体名のカンマ区切り。

    Returns: ``{"ok", "broken", "substituted", "faces": [...]}``。``broken`` の書体は
    ``doc_convert_session_font_repair`` で直せる（再起動不要）。
    """
    if sys.platform != "win32":
        return {"error": "Windows only: the check reads the session font table through GDI."}
    rows = [probe_face(f.strip()) for f in faces.split(",") if f.strip()]
    broken = [r["face"] for r in rows if r.get("broken")]
    subst = [r["face"] for r in rows if r.get("substituted")]
    return {
        "ok": not broken,
        "broken": broken,
        "substituted": subst,
        "faces": rows,
        "hint": ("" if not broken else
                 "doc_convert_session_font_repair(faces=...) reloads the face for this session; "
                 "no logoff or reboot is needed. Then start a NEW PowerPoint/Word process."),
    }


def doc_convert_session_font_repair(faces: str = "Times New Roman", force: bool = False) -> dict:
    """セッションで壊れた書体を、再起動せずにフォント表へ読み直させる。

    書体のフォントファイルを ``RemoveFontResourceW`` で参照が 0 になるまで外し、
    ``AddFontResourceW`` で読み直し、``WM_FONTCHANGE`` を放送する。レジストリもディスクも
    変えないセッション内の操作で、健全な書体に対しては何も壊さない（外してから戻すまでの
    数ミリ秒だけその書体が使えないので、書き出しの最中には実行しない）。既定では壊れている
    書体だけを直し、``force=True`` で健全な書体も読み直す。

    Args:
        faces: 直す書体名のカンマ区切り（既定 Times New Roman）。
        force: 壊れていない書体も読み直す。

    Returns: 書体ごとの ``before`` / ``after`` の GDI 検査結果と、外した・戻したファイル。
    修復後に新しく起動した Office プロセスから使えるようになる。
    """
    if sys.platform != "win32":
        return {"error": "Windows only."}
    import ctypes
    import ctypes.wintypes as w
    import time
    gdi, user = _gdi()
    results = []
    for face in [f.strip() for f in faces.split(",") if f.strip()]:
        before = probe_face(face)
        if not before.get("broken") and not force:
            results.append({"face": face, "before": before, "action": "healthy, left alone"})
            continue
        files = _font_files(face)
        if not files:
            results.append({"face": face, "before": before, "action": "no font file found for this face; nothing done"})
            continue
        removed = {}
        for p in files:
            n = 0
            while gdi.RemoveFontResourceW(str(p)) and n < 50:
                n += 1
            removed[p.name] = n
        added = {p.name: int(gdi.AddFontResourceW(str(p))) for p in files}
        res = w.DWORD()
        user.SendMessageTimeoutW(_HWND_BROADCAST, _WM_FONTCHANGE, 0, 0, _SMTO_ABORTIFHUNG, 2000, ctypes.byref(res))
        time.sleep(0.5)
        after = probe_face(face)
        results.append({"face": face, "before": before, "after": after,
                        "removed_references": removed, "added": added,
                        "action": "repaired" if after.get("ok") else "reloaded, but the face still fails: log off or reboot"})
    return {
        "ok": all(r.get("action") in ("repaired", "healthy, left alone") for r in results),
        "results": results,
        "note": "Office processes that were already running keep their own broken face references; start a new one.",
    }
