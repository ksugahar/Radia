"""MCP server: radia_mcp.md2html.

Exposes md2html_convert(md_file, output_file=None, title=None) for
turning a Markdown file into a self-contained HTML file with full
MathJax / table / code-highlight / image-embed support.

Usage:
    mcp-server-md2html              # stdio
    mcp-server-md2html --selftest   # smoke test (no Mathematica needed)
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool

from . import tools as _tools

mcp = FastMCP("mcp-server-md2html")


# ============================================================
# Tool registration (auto-pick md2html_* callables)
# ============================================================

_REGISTERED: list[str] = []
for _name in dir(_tools):
    if _name.startswith("md2html_") and callable(getattr(_tools, _name)):
        mcp.tool()(getattr(_tools, _name))
        _REGISTERED.append(_name)


# ============================================================
# Self-introspection
# ============================================================

register_status_tool(
    mcp,
    server_name="mcp-server-md2html",
    description=(
        "Markdown -> self-contained HTML with MathJax v3 + tables + "
        "fenced code + codehilite + base64 image embed + [N] reference"
        " linking.  Promoted from mcp-server-document."
    ),
    subpackage="radia_mcp.md2html",
    related_servers=["mathematica", "figure"],
    optional_deps=["markdown"],
)


# ============================================================
# Entry point
# ============================================================


def main():
    if "--selftest" in sys.argv:
        print("md2html MCP server self-test:")
        print(f"  Registered tools ({len(_REGISTERED)}):")
        for n in sorted(_REGISTERED):
            print(f"    - {n}")

        # Try a real conversion if `markdown` is importable; otherwise
        # just report that the module loaded.
        try:
            import markdown  # noqa: F401
            with tempfile.TemporaryDirectory() as td:
                src = Path(td) / "selftest.md"
                src.write_text(
                    "# md2html selftest\n\n"
                    "Inline math: $E = m c^2$.\n\n"
                    "Block math:\n\n"
                    "$$\\nabla \\times \\mathbf{H} = \\mathbf{J}\\,.$$\n\n"
                    "Reference [1].\n\n"
                    "## References\n\n"
                    "1. Ampere et al, 1820.\n",
                    encoding="utf-8")
                out = Path(td) / "selftest.html"
                result = _tools.md2html_convert(
                    str(src), output_file=str(out))
                print(f"  Test convert:")
                for line in result.splitlines():
                    print(f"    {line}")
                assert out.exists(), "output HTML not produced"
                html = out.read_text(encoding="utf-8-sig")
                assert "MathJax" in html, "MathJax block missing"
                assert "$$\\nabla" in html or "%%MATH_BLOCK_" not in html, (
                    "math placeholders leaked into HTML")
                assert "ref1" in html, "[N] reference linker did not fire"
                print(f"    output size: {out.stat().st_size} bytes")
        except ImportError:
            print("  [skipped real convert: 'markdown' package not "
                  "installed]")

        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
