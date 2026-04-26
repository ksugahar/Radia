"""
md2html.py knowledge base - Markdown to HTML converter with MathJax support.

Location: C:/Program Files/Python312/Scripts/md2html.py
Command:  md2html <input.md> [output.html] [title]

This module stores the full implementation details so the MCP server can
guide users on usage, customization, and troubleshooting.
"""

MD2HTML_USAGE = """
# md2html - Markdown to HTML Converter

## Location
- Script: `C:/Program Files/Python312/Scripts/md2html.py`
- Wrapper command: `md2html` (available on PATH)

## Usage

```bash
# Auto-generate output filename and title from input
md2html README.md
# -> README.html (title extracted from first <h1>)

# Explicit output file and title
md2html input.md output.html "My Document Title"

# From Python
python "C:/Program Files/Python312/Scripts/md2html.py" input.md
```

## Features
1. **MathJax support** - LaTeX inline `$...$` and display `$$...$$` math
2. **GitHub-style ```math blocks** - Converted to `$$...$$` automatically
3. **Table support** - `markdown.extensions.tables`
4. **Code highlighting** - `markdown.extensions.codehilite`
5. **Reference links** - `[1]` notation becomes clickable links to `#ref1`
6. **Reference list anchors** - `<ol>` items in References section get `id="refN"`
7. **Norm notation** - `||H||` auto-converted to `\\Vert H \\Vert` in math
8. **Image embedding** - Local images (PNG, JPG, etc.) are base64-encoded into data URIs, producing self-contained HTML
9. **Fallback encoding** - UTF-8 primary, cp932 fallback (Japanese Windows)
10. **BOM output** - Writes UTF-8 with BOM (`utf-8-sig`) for Windows compatibility

## Dependencies
- `markdown` (Python package: `pip install markdown`)
- MathJax 3 loaded from CDN (requires internet for math rendering)
- `base64`, `mimetypes` (stdlib, for image embedding)
"""

MD2HTML_IMPLEMENTATION = r"""
# md2html Implementation Details

## Architecture

```
Input .md -> protect_math() -> markdown.markdown() -> restore_math()
          -> add_reference_ids() -> convert_reference_links()
          -> embed_images() -> wrap in HTML template with MathJax + CSS -> output .html
```

## Key Functions

### `protect_math(md_content) -> (content, math_blocks)`
Protects math expressions from markdown processing by replacing them
with `%%MATH_BLOCK_N%%` placeholders before markdown conversion.

Processing order (important!):
1. ````math` code blocks (GitHub style) -> converted to `$$...$$`
2. Display math `$$...$$`
3. Inline math `$...$`

Also converts `||...||` to `\Vert ... \Vert` inside math (avoids table conflicts).

### `restore_math(html_content, math_blocks) -> html`
Restores `%%MATH_BLOCK_N%%` placeholders back to original math expressions.

### `convert_reference_links(html_content) -> html`
Converts `[N]` notation (where N is a number) to clickable `<a href="#refN">`
links. Skips `[N]` already inside anchor tags.

### `add_reference_ids(html_content) -> html`
Finds "References" or "参考文献" section header, then adds `id="refN"` to
each `<li>` in the first `<ol>` after it. Enables anchor scrolling from
reference links.

### `read_file_with_fallback(file_path) -> str`
Reads file with UTF-8, falls back to cp932 for Japanese files on Windows.

### `embed_images(html_content, base_dir) -> html`
Finds all `<img src="...">` tags in the HTML. For each local image
(not http/https/data URI), reads the file, base64-encodes it, and
replaces the `src` with a `data:` URI. Skips missing files with a warning.
Uses `mimetypes.guess_type()` for MIME detection.

This makes the output HTML **self-contained** -- no external image files needed.
Prints each embedded image filename and size for progress tracking.

### `md_to_html(md_file, output_file=None, title=None) -> output_path`
Main conversion function. If `output_file` is None, derives from input
(`.md` -> `.html`). If `title` is None, extracts from first `<h1>` tag.
Automatically embeds all local images referenced by `![alt](path)` markdown.

## HTML Template

- MathJax 3 configuration: inline `$...$`, display `$$...$$`, processEscapes
- Styled with Segoe UI font, max-width 900px, responsive
- Dark code blocks (`#2c3e50` background)
- Blue-themed headings and table headers
- Reference links in red (`#e74c3c`), highlighted on `:target`
- Output encoding: `utf-8-sig` (BOM)

## MathJax Configuration

```javascript
MathJax = {
  tex: {
    inlineMath: [['$', '$'], ['\\(', '\\)']],
    displayMath: [['$$', '$$'], ['\\[', '\\]']],
    processEscapes: true
  },
  options: {
    skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre']
  }
};
```

## CSS Styling Summary

| Element | Style |
|---------|-------|
| Body | Segoe UI, max-width 900px, #fafafa background |
| h1 | #2c3e50, 3px blue bottom border |
| h2 | #34495e, 1px gray bottom border, margin-top 40px |
| Code (inline) | #ecf0f1 background, Consolas font |
| Code (block) | #2c3e50 dark background, #ecf0f1 text |
| Table header | #3498db blue background, white text |
| Table even rows | #ecf0f1 alternating |
| Blockquote | 4px blue left border, #ecf0f1 background |
| Reference links | #e74c3c red, bold |
| Reference items | White cards with 3px blue left border |
| Reference :target | #fffacd yellow highlight, red border |

## Markdown Extensions Used

```python
markdown.markdown(content, extensions=['tables', 'fenced_code', 'codehilite'])
```
"""

MD2HTML_TIPS = """
# md2html Tips & Troubleshooting

## Common Patterns

### LaTeX Math in Markdown
```markdown
Inline: The impedance is $Z = R + j\\omega L$.

Display:
$$
Z_{eff}(f) = R_{DC} \\cdot F_R(\\delta) + j\\omega L
$$

GitHub-style math block:
```math
\\nabla \\times H = J + \\frac{\\partial D}{\\partial t}
```
```

### Norms in Tables
Use `||H||` notation - automatically converted to `\\Vert H \\Vert`:
```markdown
| Quantity | Formula |
|----------|---------|
| Error | $||H_1 - H_2||$ |
```

### Reference Links
```markdown
The PEEC method [1] uses partial inductance [2].

## References

1. Ruehli, A. E., "Equivalent circuit models..." (1974)
2. Rosa, E. B., "The self and mutual inductances..." (1908)
```
`[1]` and `[2]` become clickable, scrolling to the numbered reference.

## Troubleshooting

- **Math not rendering**: Check internet connection (MathJax CDN).
  For offline use, download MathJax and change the `<script src=...>`.
- **Table conflicts with ||**: The `||` in markdown creates table cells.
  In math blocks, `||...||` is auto-converted to `\\Vert`. Outside math,
  avoid `||` on lines that could be parsed as table rows.
- **Encoding errors**: File is read with UTF-8 fallback to cp932.
  Output is always UTF-8 with BOM.
- **Missing codehilite styles**: `codehilite` extension generates CSS
  classes but the template uses inline dark-theme styling for `<pre>`.
- **Large HTML output**: Images are base64-encoded (33% size increase).
  A README with 7 PNGs (~800 KB total) produces ~1.1 MB HTML.
  This is intentional for self-contained distribution.
- **Image not found warning**: If a `![](path)` references a missing file,
  the original `src` is kept as-is (not embedded). Check the path is
  relative to the `.md` file location.
"""


def get_md2html_documentation(topic: str = "all") -> str:
    """Return md2html documentation for the given topic."""
    topic = topic.lower().strip()

    sections = {
        "usage": MD2HTML_USAGE,
        "implementation": MD2HTML_IMPLEMENTATION,
        "tips": MD2HTML_TIPS,
    }

    if topic == "all":
        return "\n\n---\n\n".join([MD2HTML_USAGE, MD2HTML_IMPLEMENTATION, MD2HTML_TIPS])
    elif topic in sections:
        return sections[topic]
    else:
        return (
            f"Unknown topic: '{topic}'. "
            "Available topics: all, usage, implementation, tips"
        )
