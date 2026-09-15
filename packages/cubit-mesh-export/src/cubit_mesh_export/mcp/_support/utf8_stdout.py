"""Force stdout/stderr to UTF-8 so --selftest output doesn't crash on cp932.

Default Windows consoles in Japanese locale use cp932, which cannot encode
em-dashes (U+2014), arrows (U+2192), or CJK characters.  Knowledge-server
--selftest paths print documentation text that frequently contains such
characters, so without this helper `mcp-server-radia-ngsolve --selftest`
and friends crash with UnicodeEncodeError in Japanese Windows consoles or
service-hosted selftest runs.

Call ``use_utf8_stdout()`` once at the top of any ``main()`` path that
writes user-facing text.  MCP stdio (``mcp.run(transport="stdio")``) is
independent of this and must not be touched — JSON-RPC framing uses raw
bytes over the FastMCP stdio transport.
"""
from __future__ import annotations

import sys


def use_utf8_stdout() -> None:
    """Reconfigure ``sys.stdout``/``sys.stderr`` to UTF-8 with replace errors.

    No-op if the streams don't expose ``reconfigure`` (e.g., redirected to
    a file with a fixed encoding, or a non-TextIOWrapper).  Never raises.
    """
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is None:
            continue
        try:
            reconfigure(encoding="utf-8", errors="replace")
        except (OSError, ValueError):
            pass
