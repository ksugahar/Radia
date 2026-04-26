"""Generate docs/TOOLS.md — single-file inventory of all radia-mcp tools.

radia-mcp ships 9 separate MCP servers (cubit / build123d / radia-ngsolve /
gmsh / elf / electromagnet / ih / peec / interop), each with its own
FastMCP instance and console-script entry point. This script imports
each server module, enumerates its tools, and writes a unified markdown
inventory so an AI agent can grep one file to find any tool across the
9 servers.

Usage:
    python scripts/gen_tools_doc.py            # write docs/TOOLS.md
    python scripts/gen_tools_doc.py --check    # exit 1 if file is stale
"""
from __future__ import annotations

import asyncio
import importlib
import sys
from pathlib import Path

PKG_ROOT = Path(__file__).resolve().parent.parent
DOC_PATH = PKG_ROOT / "docs" / "TOOLS.md"

# (subpackage_name, console_script, header_blurb)
SERVERS = [
    ("cubit",          "mcp-server-cubit",          "Cubit hex-mesh export, Netgen/NGSolve curving, scripting + API reference"),
    ("build123d",      "mcp-server-build123d",      "build123d (Pythonic OCCT) + STEP/XCAF labels + Cubit pipeline interop"),
    ("radia_ngsolve",  "mcp-server-radia-ngsolve",  "Radia + NGSolve coupled magnetostatics, Kelvin transformation, sparse solver"),
    ("gmsh",           "mcp-server-gmsh",           "Gmsh script linting + post-processing spec helpers"),
    ("elf",            "mcp-server-elf",            "ELF (Electromagnetic Loss/Field) postprocessing knowledge"),
    ("electromagnet",  "mcp-server-electromagnet",  "Electromagnet design (symmetry reductions, BC choices)"),
    ("ih",             "mcp-server-ih",             "IH (induction-heating) coil + load workflow"),
    ("peec",           "mcp-server-peec",           "PEEC (partial element equivalent circuit) inductance modeling"),
    ("interop",        "mcp-server-radia-interop",  "Cross-tool interop (CadQuery / build123d / Cubit STEP boundary)"),
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
    module = importlib.import_module(f"radia_mcp.{subpkg}.server")
    return asyncio.run(module.mcp.list_tools())


def render() -> str:
    out: list[str] = []
    out.append("# radia-mcp Tools Inventory")
    out.append("")
    out.append(
        "Auto-generated from each server's `mcp.list_tools()` via "
        "`scripts/gen_tools_doc.py`. **Do not edit by hand** — "
        "regenerate after adding/renaming tools."
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

    return "\n".join(out) + "\n"


def write() -> Path:
    DOC_PATH.parent.mkdir(parents=True, exist_ok=True)
    DOC_PATH.write_text(render(), encoding="utf-8")
    return DOC_PATH


def check() -> bool:
    if not DOC_PATH.exists():
        return False
    return DOC_PATH.read_text(encoding="utf-8") == render()


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
