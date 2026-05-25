"""Tool functions for the md2html domain of mcp-server-document.

All functions are plain callables; they get wrapped by ``@mcp.tool()``
in the sub-module's ``register()`` (called from the top-level
``server.py``). No FastMCP import here (that would accidentally
register a second MCP server).
"""
from __future__ import annotations

import pathlib

from .converter import md_to_html as _md_to_html


def md2html_convert(md_file: str,
                    output_file: str | None = None,
                    title: str | None = None) -> str:
    """Convert a Markdown file to a self-contained HTML file.

    The output HTML embeds MathJax (CDN), inlines local images as
    base64 data URIs, and preserves Markdown tables, fenced code,
    code-highlighting, and ``[N]`` reference links.

    Math support:
        - inline ``$..$``
        - display ``$$..$$``
        - GitHub-style ``` ```math ``` ``` blocks
        - ``||x||`` is auto-converted to ``\\Vert x \\Vert`` so it does
          not collide with Markdown table column separators.

    Args:
        md_file: Path of the input ``.md`` file (absolute or relative
                 to the current working directory).
        output_file: Optional output ``.html`` path. Defaults to
                     ``<md_file>.html`` next to the input.
        title: Optional ``<title>`` for the HTML head. Defaults to
               the first ``<h1>`` text in the document, or the input
               file's basename if no ``<h1>`` is present.

    Returns:
        Status string with output path, number of math blocks, number
        of embedded images, and any embedding warnings.
    """
    src = pathlib.Path(md_file)
    if not src.is_absolute():
        src = pathlib.Path.cwd() / src
    if not src.exists():
        return f"Error: file not found: {src}"
    if src.suffix.lower() != '.md':
        return f"Error: not a .md file: {src}"

    out: pathlib.Path | None = None
    if output_file:
        out = pathlib.Path(output_file)
        if not out.is_absolute():
            out = pathlib.Path.cwd() / out

    result = _md_to_html(str(src), str(out) if out else None, title)

    lines = [
        f"Converted: {src} -> {result['output_file']}",
        f"  math blocks: {result['math_blocks']}",
        f"  images embedded: {result['images_embedded']}",
    ]
    if result['log'].strip():
        lines.append("  details:")
        for ln in result['log'].splitlines():
            lines.append(f"    {ln}")
    return "\n".join(lines)
