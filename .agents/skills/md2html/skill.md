---
name: md2html
description: Compatibility guidance for Markdown to HTML conversion through radia-mcp.
---

# Markdown to HTML

Use the `radia-mcp` tool `md2html_convert` or its public Python helper:

```python
from radia_mcp.md2html import md_to_html
result = md_to_html("input.md", "output.html", "Document title")
print(result["output_file"])
```

The helper returns a dictionary containing the output path, math/image counts,
and conversion log. Read warnings about missing local images before accepting
the output. MathJax rendering requires access to its CDN.

The canonical implementation is
`packages/radia-mcp/src/radia_mcp/md2html/converter.py`.
Do not use old workstation script paths or assume an `md2html` executable exists.

For callers of the old repository CLI, `md2html.py` in this directory delegates
to the same installed package and retains its path-string return value:

```powershell
python .agents/skills/md2html/md2html.py input.md output.html "Document title"
```

Install `radia-mcp[md2html]` in the chosen interpreter if needed. Missing
dependencies fail visibly; do not restore a second converter implementation or
repoint a shared editable installation merely to run this command.
