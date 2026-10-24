"""Human Netgen .vol launcher; no meshing, solving, file writes or install patches."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import sys


def validate_mesh_path(filename: str) -> Path:
    """Reject absent/non-mesh inputs before importing a native GUI runtime."""
    path = Path(filename).expanduser().resolve(strict=True)
    if not path.is_file() or path.suffix.lower() != ".vol":
        raise ValueError("Specify one existing .vol file; other extensions are unsupported")
    return path


def _file_identity(filename) -> dict:
    path = Path(filename).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def runtime_identity(path: Path) -> dict:
    """Read identities only; never run a write/read qualification probe."""
    requirements = importlib.metadata.requires('radia') or []
    pins = [match.group(1) for item in requirements
            if (match := re.fullmatch(r'netgen[-_]mesher\s*==\s*([^;\s]+)(?:\s*;.*)?', item, re.I))]
    if len(pins) != 1:
        raise RuntimeError('Installed radia metadata must declare one exact netgen-mesher pin')
    version = importlib.metadata.version("netgen-mesher")
    if version != pins[0]:
        raise RuntimeError(f"Netgen {pins[0]} required by radia metadata; installed {version}")
    import netgen
    if netgen.__version__ != version:
        raise RuntimeError("Netgen module/distribution version mismatch")
    runtime = {"netgen_version": version, "module": _file_identity(netgen.__file__),
               "native_module": _file_identity(netgen.libngpy.__file__)}
    backend = "netgen"
    return {"backend": backend, "runtime": runtime, "python": sys.executable,
            "python_version": sys.version.split()[0],
            "launcher": _file_identity(__file__), "mesh": str(path)}


def open_mesh(path: Path) -> None:
    """Display data in the existing viewer APIs without modifying installations."""
    from netgen.gui import win
    closed = False
    def mark_root_destroyed(event):
        nonlocal closed
        if event.widget is win:
            closed = True
    win.bind("<Destroy>", mark_root_destroyed, add="+")
    def close_owned_window():
        win.destroy()
    win.protocol("WM_DELETE_WINDOW", close_owned_window)
    try:
        # Filename is exactly one Tcl argument, never an interpolated script.
        win.tk.call("Ng_LoadMesh", path.as_posix())
        win.tk.call("set", "selectvisual", "mesh")
        win.tk.call("Ng_SetVisParameters")
        win.tk.call("redraw")
        win.tk.call("Ng_ReadStatus")
        win.title("Mesh viewer - " + path.name)
        win.mainloop()
    finally:
        if not closed:
            win.destroy()


def require_interactive_desktop() -> dict:
    """Fail for services/locked desktops instead of reporting an invisible GUI."""
    if sys.platform != "win32":
        raise RuntimeError("This launcher requires an interactive Windows desktop")
    import ctypes as c
    from ctypes import wintypes as w
    kernel = c.WinDLL("kernel32", use_last_error=True)
    user = c.WinDLL("user32", use_last_error=True)
    session = w.DWORD()
    if not kernel.ProcessIdToSessionId(os.getpid(), c.byref(session)) or not session.value:
        raise RuntimeError("A non-service interactive session is required")
    user.OpenInputDesktop.argtypes = [w.DWORD, w.BOOL, w.DWORD]
    user.OpenInputDesktop.restype = w.HANDLE
    user.GetUserObjectInformationW.argtypes = [w.HANDLE, c.c_int, w.LPVOID, w.DWORD, c.POINTER(w.DWORD)]
    user.CloseDesktop.argtypes = [w.HANDLE]
    desktop = user.OpenInputDesktop(0, False, 1)
    if not desktop:
        raise RuntimeError("Interactive input desktop is unavailable; unlock the desktop")
    try:
        name = c.create_unicode_buffer(256)
        needed = w.DWORD()
        if not user.GetUserObjectInformationW(desktop, 2, name, c.sizeof(name), c.byref(needed)) or name.value.lower() != "default":
            raise RuntimeError("The interactive desktop is not unlocked")
        return {"session_id": session.value, "desktop": name.value}
    finally:
        user.CloseDesktop(desktop)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh")
    parser.add_argument("--check", action="store_true", help="read-only preflight; no GUI")
    args = parser.parse_args(argv)
    try:
        path = validate_mesh_path(args.mesh)
        desktop = require_interactive_desktop()
        identity = {**runtime_identity(path), **desktop}
        if args.check:
            print(json.dumps(identity, ensure_ascii=True))
        else:
            open_mesh(path)
        return 0
    except Exception as error:
        if args.check or sys.platform != "win32":
            print(str(error), file=sys.stderr)
        else:
            import ctypes
            ctypes.windll.user32.MessageBoxW(None, str(error), "Mesh viewer error", 0x10)
        return 1
