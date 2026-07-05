"""radia_mcp.md2html — Markdown -> self-contained HTML with MathJax.

Promoted on 2026-05-26 from
  legacy private source tree
into radia-mcp.  Source files (converter.py, tools.py) carried over
byte-identical; only the registration layer (server.py) is new and
follows the standard radia_mcp.<topic>.server pattern (statusable,
--selftest-able, discoverable via the meta catalog).

One MCP tool:
    md2html_convert(md_file, output_file=None, title=None)

Features (preserved verbatim from the originally-promoted skill):
    - MathJax v3 CDN for LaTeX math: inline $..$, display $$..$$,
      GitHub-style ```math``` blocks
    - Markdown tables, fenced code, codehilite syntax highlighting
    - ||x|| -> \\Vert x \\Vert to avoid the Markdown table column-
      separator | clash
    - [N] reference linking + auto <li id="refN"> on References / 参考文献
    - Local images embedded as base64 data URIs (single-file HTML)
    - UTF-8 first, cp932 fallback for legacy Japanese files
"""

from .converter import md_to_html  # noqa: F401  -- public helper
from .tools import md2html_convert  # noqa: F401
