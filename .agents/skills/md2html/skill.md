---
name: md2html
description: DEPRECATED — moved to mcp-server-document as md2html_convert tool
allowed-tools: Bash(python *), Read
---

# md2html — moved to mcp-server-document

This skill has been promoted to a tool in **mcp-server-document** on
2026-05-03 so it is reachable from any MCP-aware client without going
through the slash-command skill mechanism.

## New invocation

Through MCP (preferred):

- Server: `mcp-server-document`
- Tool: `md2html_convert(md_file, output_file=None, title=None)`

The tool returns a status string with output path, math-block count,
embedded image count, and per-image embedding warnings.

Source: `S:\mcp-server\src\mcp_server_document\md2html\`
- `tools.py` — MCP tool wrapper (`md2html_convert`)
- `converter.py` — pure conversion core (`md_to_html`)

## Legacy CLI (still works)

The original standalone script remains at this folder and can be run
directly when MCP is unavailable:

```bash
python md2html.py <input.md> [<output.html>] [<title>]
```

## Why moved

- Single MCP server (`mcp-server-document`) already handles
  paper_writing / presentation / diagram / ocr / etc., and md2html is
  a natural fit for the document family.
- MCP tools are usable from any client (Codex Desktop, Codex,
  custom integrations); skills are Codex only.
- Pure-function `converter.md_to_html` is now importable as a regular
  Python API for other tools.
