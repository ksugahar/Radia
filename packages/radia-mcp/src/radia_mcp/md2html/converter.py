"""Markdown -> self-contained HTML with MathJax (pure-function core).

Ported from .claude/skills/md2html/md2html.py (2026-05-03). Behaviour is
preserved verbatim; only the entry point changes (callable function +
module-level template instead of a CLI script).

Features:
    - MathJax v3 (CDN) for LaTeX math: inline ``$..$``, display ``$$..$$``,
      and GitHub-style ``` ```math ``` ``` blocks
    - Markdown tables, fenced code, codehilite syntax highlighting
    - Auto-conversion of ``||x||`` to ``\\Vert x \\Vert`` to avoid the
      Markdown table column-separator ``|`` clash
    - ``[N]`` citations turned into ``<a href="#refN">`` clickable links;
      ``<li>`` items under "References" / "参考文献" get matching IDs
    - Local images embedded as base64 data URIs (single-file HTML)
    - UTF-8 input with cp932 fallback for legacy Japanese files
    - UTF-8 BOM (``utf-8-sig``) output
"""
from __future__ import annotations

import base64
import io
import mimetypes
import os
import re

# `markdown` is lazy-imported inside md_to_html(); it's the only
# consumer.  Keeping it out of the module-level import-chain lets
# `import radia_mcp.md2html.server` succeed even when the optional
# `markdown` package isn't installed (the [md2html] extra), which
# matters for the radia-mcp meta-health selftest (CI).  Mirrors the
# build123d/server.py pattern for optional heavy deps.


# ---------------------------------------------------------------------------
# Math protection
# ---------------------------------------------------------------------------

def _convert_norm_notation(math_content: str) -> str:
    """``||x||`` -> ``\\Vert x \\Vert`` to avoid Markdown table conflict."""
    result = re.sub(r'\\\|\\\|([^|]+)\\\|\\\|', r'\\Vert \1 \\Vert', math_content)
    result = re.sub(r'\|\|([^|]+)\|\|', r'\\Vert \1 \\Vert', result)
    return result


def _protect_math(md_content: str) -> tuple[str, list[str]]:
    """Replace math expressions with placeholders so Markdown leaves them alone."""
    math_blocks: list[str] = []

    def save_math(match: re.Match) -> str:
        converted = _convert_norm_notation(match.group(0))
        math_blocks.append(converted)
        return f"%%MATH_BLOCK_{len(math_blocks)-1}%%"

    def save_math_codeblock(match: re.Match) -> str:
        math_content = match.group(1).strip()
        converted = _convert_norm_notation(f"$$\n{math_content}\n$$")
        math_blocks.append(converted)
        return f"%%MATH_BLOCK_{len(math_blocks)-1}%%"

    # ```math ``` (GitHub style) first, then $$..$$, then $..$
    content = re.sub(r'```math\s*\n([\s\S]*?)\n```', save_math_codeblock, md_content)
    content = re.sub(r'\$\$[\s\S]*?\$\$', save_math, content)
    content = re.sub(r'\$[^\$\n]+?\$', save_math, content)
    return content, math_blocks


def _restore_math(html_content: str, math_blocks: list[str]) -> str:
    for i, math in enumerate(math_blocks):
        html_content = html_content.replace(f"%%MATH_BLOCK_{i}%%", math)
    return html_content


# ---------------------------------------------------------------------------
# Reference linking
# ---------------------------------------------------------------------------

def _convert_reference_links(html_content: str) -> str:
    """``[N]`` -> ``<a href="#refN" class="ref-link">[N]</a>``."""
    html_content = re.sub(
        r'((?:<[^>]*>|[^<])*)(\[\d+\])',
        lambda m: m.group(1) + re.sub(
            r'\[(\d+)\]', r'<a href="#ref\1" class="ref-link">[\1]</a>', m.group(2)
        ) if 'href=' not in m.group(0) and '#ref' not in m.group(0) else m.group(0),
        html_content,
    )
    html_content = re.sub(
        r'(?<!href="#ref)(?<!">)\[(\d+)\](?!</a>)',
        r'<a href="#ref\1" class="ref-link">[\1]</a>',
        html_content,
    )
    return html_content


def _add_reference_ids(html_content: str) -> str:
    """``id="refN"`` on each <li> in the References / 参考文献 section."""
    ref_section_match = re.search(r'(参考文献|References)', html_content, re.IGNORECASE)
    if not ref_section_match:
        return html_content

    ref_pos = ref_section_match.end()
    before_ref = html_content[:ref_pos]
    after_ref = html_content[ref_pos:]

    ol_match = re.search(r'<ol>([\s\S]*?)</ol>', after_ref)
    if not ol_match:
        return html_content

    ol_start = ol_match.start()
    ol_end = ol_match.end()
    ol_content = ol_match.group(1)

    li_count = [0]
    def add_li_id(_: re.Match) -> str:
        li_count[0] += 1
        return f'<li id="ref{li_count[0]}">'

    ol_with_ids = re.sub(r'<li>', add_li_id, ol_content)
    after_ref = (after_ref[:ol_start] +
                 '<ol class="references">\n' + ol_with_ids + '</ol>' +
                 after_ref[ol_end:])
    return before_ref + after_ref


# ---------------------------------------------------------------------------
# Image embedding
# ---------------------------------------------------------------------------

def _embed_images(html_content: str, base_dir: str, log: io.StringIO) -> tuple[str, int]:
    """Replace local <img src="..."> with base64 data URIs."""
    img_count = 0

    def replace_img_src(match: re.Match) -> str:
        nonlocal img_count
        full_tag = match.group(0)
        src = match.group(1)

        if src.startswith(('http://', 'https://', 'data:')):
            return full_tag

        img_path = os.path.join(base_dir, src)
        if not os.path.isfile(img_path):
            log.write(f"  Warning: image not found: {img_path}\n")
            return full_tag

        mime_type, _ = mimetypes.guess_type(img_path)
        if mime_type is None:
            mime_type = 'application/octet-stream'

        with open(img_path, 'rb') as f:
            img_data = base64.b64encode(f.read()).decode('ascii')

        img_count += 1
        size_kb = os.path.getsize(img_path) / 1024
        log.write(f"  Embedded: {src} ({size_kb:.0f} KB)\n")
        return full_tag.replace(f'src="{src}"', f'src="data:{mime_type};base64,{img_data}"')

    html_content = re.sub(r'<img\s[^>]*src="([^"]+)"[^>]*>', replace_img_src, html_content)
    return html_content, img_count


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _read_file_with_fallback(file_path: str) -> str:
    """UTF-8 first, fall back to cp932 (Shift-JIS) for legacy JP files."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return f.read()
    except UnicodeDecodeError:
        with open(file_path, 'r', encoding='cp932') as f:
            return f.read()


# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------

_HTML_TEMPLATE = '''<!DOCTYPE html>
<html lang="ja">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <script>
    MathJax = {{
      tex: {{
        inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
        displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
        processEscapes: true
      }},
      options: {{
        skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre']
      }}
    }};
    </script>
    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                line-height: 1.8; max-width: 900px; margin: 0 auto;
                padding: 20px; background-color: #fafafa; color: #333; }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db;
              padding-bottom: 10px; }}
        h2 {{ color: #34495e; border-bottom: 1px solid #bdc3c7;
              padding-bottom: 5px; margin-top: 40px; }}
        h3 {{ color: #555; margin-top: 30px; }}
        h4 {{ color: #666; margin-top: 20px; }}
        code {{ background-color: #ecf0f1; padding: 2px 6px;
                border-radius: 3px;
                font-family: 'Consolas', 'Monaco', monospace; }}
        pre {{ background-color: #2c3e50; color: #ecf0f1;
               padding: 15px; border-radius: 5px; overflow-x: auto; }}
        pre code {{ background-color: transparent; padding: 0; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #bdc3c7; padding: 10px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
        tr:nth-child(even) {{ background-color: #ecf0f1; }}
        blockquote {{ border-left: 4px solid #3498db; margin: 20px 0;
                      padding: 10px 20px; background-color: #ecf0f1; }}
        .highlight {{ background-color: #fffacd; padding: 2px 5px;
                      border-radius: 3px; }}
        ul, ol {{ padding-left: 30px; }}
        li {{ margin: 8px 0; }}
        a {{ color: #3498db; text-decoration: none; }}
        a:hover {{ text-decoration: underline; }}
        hr {{ border: none; border-top: 1px solid #ddd; margin: 30px 0; }}
        .ref-link {{ color: #e74c3c; font-weight: bold; text-decoration: none; }}
        .ref-link:hover {{ text-decoration: underline; }}
        ol.references {{ background-color: #f5f5f5;
                          padding: 20px 20px 20px 45px; border-radius: 8px; }}
        ol.references li {{ padding: 10px; background-color: #fff;
                            border-radius: 5px; border-left: 3px solid #3498db;
                            margin: 10px 0; }}
        ol.references li:target {{ background-color: #fffacd;
                                    border-left-color: #e74c3c; }}
        .MathJax_Display {{ background-color: #f8f9fa; padding: 10px;
                             margin: 15px 0; border-radius: 5px; overflow-x: auto; }}
    </style>
</head>
<body>
{body}
</body>
</html>
'''


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def md_to_html(md_file: str, output_file: str | None = None,
               title: str | None = None) -> dict:
    """Convert a Markdown file to a self-contained HTML file.

    Parameters
    ----------
    md_file
        Path to the input ``.md`` file.
    output_file
        Path to the output ``.html`` file. Defaults to ``<md_file>.html``.
    title
        ``<title>`` text. Defaults to the first ``<h1>`` text, or the
        input file's basename if no ``<h1>`` is present.

    Returns
    -------
    dict with keys ``output_file``, ``math_blocks``, ``images_embedded``,
    and ``log`` (string with embedding warnings / per-image notes).
    """
    # Lazy import: the `markdown` package is an optional dep of
    # radia-mcp ([md2html] extra).  Importing at module-top would
    # break the meta-health selftest that runs without optional deps.
    import markdown

    if output_file is None:
        output_file = md_file.replace('.md', '.html')

    md_content = _read_file_with_fallback(md_file)

    protected_content, math_blocks = _protect_math(md_content)

    html_body = markdown.markdown(
        protected_content,
        extensions=['tables', 'fenced_code', 'codehilite'],
    )

    html_body = _restore_math(html_body, math_blocks)
    html_body = _add_reference_ids(html_body)
    html_body = _convert_reference_links(html_body)

    log = io.StringIO()
    md_dir = os.path.dirname(os.path.abspath(md_file))
    html_body, img_count = _embed_images(html_body, md_dir, log)

    if title is None:
        title = os.path.basename(md_file).replace('.md', '')
        if '<h1>' in html_body:
            start = html_body.find('<h1>') + 4
            end = html_body.find('</h1>')
            title = html_body[start:end]

    html_content = _HTML_TEMPLATE.format(title=title, body=html_body)
    with open(output_file, 'w', encoding='utf-8-sig') as f:
        f.write(html_content)

    return {
        'output_file': output_file,
        'math_blocks': len(math_blocks),
        'images_embedded': img_count,
        'log': log.getvalue(),
    }
