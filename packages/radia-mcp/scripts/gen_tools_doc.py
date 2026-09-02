"""Generate docs/TOOLS.md -- single-file inventory of all radia-mcp tools.

radia-mcp ships one MCP server per subpackage listed in
radia_mcp.meta.catalog (cubit / build123d / radia-ngsolve / gmsh / ih /
peec / electromagnet / ... -- the count grows over time and is derived
dynamically below, never hardcoded).  Each has its own FastMCP instance
and console-script entry point.  This script imports each server module,
enumerates its tools, and writes a unified markdown inventory so an AI
agent can grep one file to find any tool across all servers.

Usage:
    python scripts/gen_tools_doc.py            # write docs/TOOLS.md
    python scripts/gen_tools_doc.py --check    # exit 1 if file is stale
"""
from __future__ import annotations

import asyncio
import importlib
import os
import sys
from pathlib import Path

PKG_ROOT = Path(__file__).resolve().parent.parent
DOC_PATH = PKG_ROOT / "docs" / "TOOLS.md"

# Auto-derive the server list from radia_mcp.meta.catalog so this
# script stays in sync with the canonical catalog. Previously this
# was a hardcoded list of 9 entries and silently fell behind the
# 27-subpackage batch additions (caught by the 2026-05-24 review).
sys.path.insert(0, str(PKG_ROOT / "src"))
# The checked inventory documents the production surface. Compatibility and
# debugging may expose legacy individual gates with the ``full`` profile, but
# that deliberately larger list is not a public API catalog.
os.environ["RADIA_MCP_TOOL_PROFILE"] = "core"
from radia_mcp.meta.catalog import CATALOG  # noqa: E402

# (subpackage_name, console_script, header_blurb)
SERVERS = [
    (info["subpackage"].replace("radia_mcp.", ""),
     info["entry_point"],
     info["description"])
    for _name, info in CATALOG.items()
]


def _first_line(description: str | None) -> str:
    if not description:
        return "(no description)"
    for line in description.splitlines():
        s = line.strip()
        if s:
            return s
    return "(no description)"


def _list_tools_for(subpkg: str) -> list:
    """Import the server module and enumerate its tools.

    Calls `mcp.list_tools()` TWICE and returns the second result.  The
    first call has a lazy-registration side effect in some FastMCP
    versions (some servers register meta tools only on the first
    list_tools() invocation), which would otherwise make the output
    non-deterministic between cold and warm runs.  Calling twice
    guarantees the warm-cache result every time.
    """
    module = importlib.import_module(f"radia_mcp.{subpkg}.server")
    asyncio.run(module.mcp.list_tools())   # warm-up; discard
    return asyncio.run(module.mcp.list_tools())


def render() -> str:
    out: list[str] = []
    out.append("# radia-mcp Tools Inventory")
    out.append("")
    out.append(
        "Auto-generated from each server's production `core` "
        "`mcp.list_tools()` via "
        "`scripts/gen_tools_doc.py`. **Do not edit by hand** — "
        "regenerate after adding/renaming tools."
    )
    out.append("")
    out.append(
        "Fine-grained validation and identity operations are discovered with "
        "each server's `*_validation_catalog` tool and invoked through "
        "`*_validation_run`; they are not repeated as top-level schemas."
    )
    out.append("")

    # Collect tools per server first so we can compute totals
    per_server: list[tuple[str, str, str, list]] = []
    grand_total = 0
    for subpkg, script, blurb in SERVERS:
        try:
            tools = _list_tools_for(subpkg)
        except Exception as e:  # noqa: BLE001
            tools = []
            blurb = f"{blurb}  _(import failed: {type(e).__name__}: {e})_"
        per_server.append((subpkg, script, blurb, tools))
        grand_total += len(tools)

    out.append(f"Total: **{grand_total} tools** across {len(SERVERS)} MCP servers.")
    out.append("")

    # TOC
    out.append("| Server (console-script) | Subpackage | Tools |")
    out.append("|---|---|---:|")
    for subpkg, script, _blurb, tools in per_server:
        anchor = script.replace("_", "-")
        out.append(f"| [`{script}`](#{anchor}) | `radia_mcp.{subpkg}` | {len(tools)} |")
    out.append("")

    # Per-server section
    for subpkg, script, blurb, tools in per_server:
        out.append(f"## `{script}`")
        out.append("")
        out.append(f"_{blurb}_")
        out.append("")
        out.append(f"Module: `radia_mcp.{subpkg}.server`")
        out.append("")
        if not tools:
            out.append("_(no tools registered or import failed)_")
            out.append("")
            continue
        out.append("| Tool | Description |")
        out.append("|---|---|")
        for t in sorted(tools, key=lambda t: t.name):
            desc = _first_line(t.description).replace("|", r"\|")
            if len(desc) > 200:
                desc = desc[:197] + "..."
            out.append(f"| `{t.name}` | {desc} |")
        out.append("")

    return "\n".join(out).rstrip() + "\n"


def write() -> Path:
    DOC_PATH.parent.mkdir(parents=True, exist_ok=True)
    # newline="" + manual \n preserves LF endings on Windows.  Without
    # this, Windows write_text translates \n -> \r\n on disk, then CI
    # (Ubuntu) reads \r\n literally and compares against render() (\n),
    # causing a spurious drift-gate failure that only fires on CI.
    with open(DOC_PATH, "w", encoding="utf-8", newline="") as fh:
        fh.write(render())
    return DOC_PATH


def check() -> bool:
    if not DOC_PATH.exists():
        return False
    # Read raw bytes, then decode -- matches what write() puts on disk
    # (LF-only). Avoids platform-dependent newline translation at READ
    # time. But the file on disk MAY have CRLF anyway because of:
    #   (a) Windows `core.autocrlf=true` (LAB + self-hosted runner default)
    #   (b) `actions/checkout@v4` on persistent runners doesn't re-write
    #       files that git considers "clean" via autocrlf, even after
    #       `.gitattributes` adds `text eol=lf` rule.
    # So normalize CRLF -> LF on BOTH sides before comparing -- the
    # drift check is about CONTENT, not about line-ending bytes.
    with open(DOC_PATH, encoding="utf-8", newline="") as fh:
        on_disk = fh.read()
    return on_disk.replace("\r\n", "\n") == render().replace("\r\n", "\n")


def main(argv: list[str]) -> int:
    if "--check" in argv:
        if check():
            print(f"OK: {DOC_PATH.relative_to(PKG_ROOT)} is up-to-date.")
            return 0
        print(
            f"STALE: {DOC_PATH.relative_to(PKG_ROOT)} differs from generated. "
            "Run `python scripts/gen_tools_doc.py` to refresh.",
            file=sys.stderr,
        )
        return 1

    p = write()
    print(f"wrote {p.relative_to(PKG_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
