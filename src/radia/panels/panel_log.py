"""Shared panel debug log writer.

All Radia GUI components — Cubit-side ``register_toolbar.py``,
external-Python ``radia_gui_base.py`` (the IH/EM/PCB analysis windows),
and the ``calc_*.py`` subprocess scripts — write to the same file:

    Windows: C:/radia_panel_log.txt
    Other:   $TMPDIR/radia_panel_log.txt

This gives the user (and Claude) **one place** to look when something
goes wrong in the GUI:

  - Did the panel dialog open?
  - Did the model labels validate?
  - Was the subprocess command line correct?
  - Did calc_inductance.py print an error?
  - What was the JSON result?

Every line is tagged with:

  - millisecond timestamp
  - source component (cubit / ih-window / inductance / fem_kelvin / ...)
  - **user@host** so logs from multiple machines (LAB, 100号機, mdx)
    or multiple users on the same machine can be told apart

Example::

    [2026-04-12 14:30:12.345] [ksugahar@LAB         ] (cubit       ) register_toolbar.py loaded
    [2026-04-12 14:30:18.892] [ksugahar@LAB         ] (cubit       ) _launch_radia_ngsolve: ENTER
    [2026-04-12 14:30:25.103] [ksugahar@LAB         ] (ih-window   ) _on_run: cmd=...
    [2026-04-12 14:30:25.567] [ksugahar@LAB         ] (inductance  ) MESH:loaded radia_model.vol
    [2026-04-12 14:30:37.842] [ksugahar@LAB         ] (inductance  ) SOLVE_DONE 12.3s
    [2026-04-12 14:30:37.901] [ksugahar@LAB         ] (ih-window   ) result: L=87.81 nH

The user@host tag is captured **once at process start** (in
``init_panel_log``) so it does not change mid-session even if
environment variables shift. This makes a multi-user log immediately
self-explanatory: Kubota's runs on 100号機 will be tagged
``[kubota@KUBOTA-PC      ]``, Sugahara's lab runs ``[ksugahar@LAB         ]``.

The log is **NOT truncated** by individual processes — only the
top-level Cubit-side ``register_toolbar.py`` truncates it on each
Cubit session start. Subprocess writes append. This way one Cubit
session produces one continuous log file across all the processes
that run during it.

Usage in any module::

    from radia.panels.panel_log import panel_log, panel_log_exception
    panel_log("hello")  # auto-tagged with the importing module name
    try:
        ...
    except Exception:
        panel_log_exception("solve failed")
"""

from __future__ import annotations

import getpass
import os
import platform
import socket
import sys
import time
import traceback


# ============================================================
# Log file path
# ============================================================
if sys.platform == "win32":
    PANEL_LOG_PATH = "C:/radia_panel_log.txt"
else:
    PANEL_LOG_PATH = os.path.join(
        os.environ.get("TMPDIR", "/tmp"), "radia_panel_log.txt")


# Component tag — set per-process via init_panel_log() or auto-detected
# from the calling module name. Width is fixed so columns line up.
_COMPONENT_TAG = "radia       "

# user@host tag — captured once at process start so it survives env
# changes mid-session. 20 chars wide (8 user + @ + 11 host) so the
# columns stay aligned in the log file. If the real user/host strings
# are longer they are truncated.
_USERHOST_TAG = ""


def _detect_user():
    """Best-effort current user name.

    Tries getpass.getuser() (which honors USER, LOGNAME, USERNAME etc.)
    and falls back through environment variables. Returns 'unknown'
    if every method fails.
    """
    try:
        u = getpass.getuser()
        if u:
            return u
    except Exception:
        pass
    for var in ("USER", "USERNAME", "LOGNAME"):
        u = os.environ.get(var)
        if u:
            return u
    return "unknown"


def _detect_host():
    """Best-effort short host name.

    Strips domain suffix (foo.example.com -> foo) and truncates to
    12 chars. Returns 'unknown' if every method fails.
    """
    try:
        h = platform.node() or socket.gethostname()
    except Exception:
        h = ""
    if not h:
        h = os.environ.get("COMPUTERNAME", "") or os.environ.get(
            "HOSTNAME", "") or "unknown"
    h = h.split(".", 1)[0]  # strip FQDN suffix
    return h


def _format_userhost():
    """Return ``user@host`` truncated/padded to 20 chars (fixed-width)."""
    u = _detect_user()
    h = _detect_host()
    full = f"{u}@{h}"
    return f"{full:<20.20s}"


def init_panel_log(component, *, truncate=False, banner=True):
    """Set the source-component tag for this process.

    Call once near the top of any Radia GUI / panel script.

    Args:
        component: short string (max 12 chars). e.g. "cubit",
            "ih-window", "calc_induct", "calc_fem".
        truncate: if True, the log file is wiped first. Use this only
            in the top-level Cubit ``register_toolbar.py`` to start a
            fresh session log. Subprocesses must NEVER truncate.
        banner: if True, write a separator + a "process started" line.
    """
    global _COMPONENT_TAG, _USERHOST_TAG
    _COMPONENT_TAG = f"{component:<12.12s}"
    _USERHOST_TAG = _format_userhost()

    if truncate:
        try:
            open(PANEL_LOG_PATH, "w", encoding="utf-8").close()
        except Exception:
            pass

    if banner:
        panel_log("=" * 70)
        panel_log(
            f"{component} started "
            f"(pid={os.getpid()}, python={sys.version.split()[0]}, "
            f"user={_detect_user()}, host={_detect_host()}, "
            f"platform={sys.platform})")


def panel_log(msg):
    """Append a single line to the panel debug log.

    Format::

        [YYYY-MM-DD HH:MM:SS.mmm] [user@host          ] (component  ) msg

    Never raises — all errors are swallowed (we cannot log a logging
    failure without infinite recursion).

    If ``init_panel_log`` was never called the user@host tag is filled
    in lazily on the first ``panel_log`` call so logs from forgotten
    init paths are still attributable.
    """
    global _USERHOST_TAG
    if not _USERHOST_TAG:
        _USERHOST_TAG = _format_userhost()
    ts = (time.strftime("%Y-%m-%d %H:%M:%S")
          + f".{int((time.time() % 1) * 1000):03d}")
    line = f"[{ts}] [{_USERHOST_TAG}] ({_COMPONENT_TAG}) {msg}\n"
    try:
        with open(PANEL_LOG_PATH, "a", encoding="utf-8") as f:
            f.write(line)
    except Exception:
        pass


def panel_log_exception(prefix=""):
    """Append the current exception's traceback to the panel debug log.

    Multi-line tracebacks are written line-by-line so they line up with
    the timestamp / component column.
    """
    tb = traceback.format_exc()
    panel_log(f"{prefix} EXCEPTION")
    for line in tb.rstrip().splitlines():
        panel_log(f"  {line}")


def panel_log_path():
    """Return the absolute path of the log file (for printing in UI)."""
    return PANEL_LOG_PATH
