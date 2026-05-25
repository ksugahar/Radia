"""Golden tests for radia_mcp.md2html (Markdown -> HTML conversion).

Promoted from mcp-server-document; locks the public-API surface so
later edits to converter.py / tools.py cannot silently drop behaviour.

Tests are skipped (NOT failed) if the `markdown` package is absent.
"""
from __future__ import annotations

import pytest

markdown = pytest.importorskip("markdown")

from radia_mcp.md2html import md2html_convert, md_to_html  # noqa: E402


def _write(p, txt: str):
    p.write_text(txt, encoding="utf-8")
    return str(p)


def test_basic_convert_creates_html(tmp_path):
    """Round-trip: .md -> .html, output file exists + non-empty."""
    src = _write(tmp_path / "hello.md", "# Hello\n\nworld\n")
    msg = md2html_convert(src)
    out = tmp_path / "hello.html"
    assert out.exists()
    assert out.stat().st_size > 500     # template alone is > 1 KB
    assert "Converted:" in msg


def test_mathjax_block_emitted(tmp_path):
    """MathJax CDN script must be in the output (math support promise)."""
    src = _write(tmp_path / "m.md",
                  "# title\n\n$E = m c^2$ and $$x = 1$$\n")
    md2html_convert(src)
    html = (tmp_path / "m.html").read_text(encoding="utf-8-sig")
    assert "MathJax" in html
    assert "tex-mml-chtml" in html


def test_double_dollar_block_preserved(tmp_path):
    """$$..$$ block must survive verbatim (not converted by Markdown)."""
    src = _write(tmp_path / "b.md",
                  "$$\n\\nabla \\times \\mathbf{H} = \\mathbf{J}\n$$\n")
    md2html_convert(src)
    html = (tmp_path / "b.html").read_text(encoding="utf-8-sig")
    assert "\\nabla" in html
    assert "\\mathbf{H}" in html


def test_norm_notation_converted(tmp_path):
    """||x|| -> \\Vert x \\Vert (Markdown table-column-separator
    collision avoidance)."""
    src = _write(tmp_path / "n.md", "$||v||^2$\n")
    md2html_convert(src)
    html = (tmp_path / "n.html").read_text(encoding="utf-8-sig")
    # original ||v|| must NOT survive (would corrupt table parsing)
    assert "||v||" not in html
    assert "\\Vert" in html


def test_reference_linking(tmp_path):
    """[1] becomes <a href="#ref1"> + <li id="ref1"> under References."""
    src = _write(tmp_path / "r.md",
                  "intro [1]\n\n## References\n\n1. Author 2024.\n")
    md2html_convert(src)
    html = (tmp_path / "r.html").read_text(encoding="utf-8-sig")
    assert 'href="#ref1"' in html
    assert 'id="ref1"' in html


def test_japanese_references_section(tmp_path):
    """参考文献 heading also triggers reference-id rewriting."""
    src = _write(tmp_path / "jp.md",
                  "本文 [1]\n\n## 参考文献\n\n1. 菅原 2026.\n")
    md2html_convert(src)
    html = (tmp_path / "jp.html").read_text(encoding="utf-8-sig")
    assert 'id="ref1"' in html


def test_local_image_embedded_as_base64(tmp_path):
    """Local <img src='foo.png'> -> data:image/png;base64,..."""
    # Minimal 1x1 PNG (transparent)
    png_bytes = (b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01"
                  b"\x00\x00\x00\x01\x08\x06\x00\x00\x00\x1f\x15\xc4\x89"
                  b"\x00\x00\x00\rIDAT\x08\x99c\xf8\x0f\x00\x00\x01\x01"
                  b"\x01\x00\x1bV\xf9\xb1\x00\x00\x00\x00IEND\xaeB`\x82")
    img = tmp_path / "pix.png"
    img.write_bytes(png_bytes)
    # converter's regex is r'<img\s[^>]*src="([^"]+)"[^>]*>' -- it
    # ONLY accepts double-quoted src (matches the standard HTML
    # serialisation Python-Markdown produces from ![alt](pix.png)).
    src = _write(tmp_path / "img.md",
                  '<img src="pix.png" alt="pix">\n')
    md2html_convert(src)
    html = (tmp_path / "img.html").read_text(encoding="utf-8-sig")
    assert "data:image/png;base64" in html


def test_cp932_fallback(tmp_path):
    """A legacy cp932-encoded .md must read correctly via fallback."""
    src = tmp_path / "legacy.md"
    src.write_text("# 古い日本語\n\nテスト\n", encoding="cp932")
    md2html_convert(str(src))
    html = (tmp_path / "legacy.html").read_text(encoding="utf-8-sig")
    # The CJK content must appear in the output
    assert "古い日本語" in html


def test_unknown_file_returns_error_string(tmp_path):
    """Missing file -> 'Error: file not found' string, NOT raise."""
    msg = md2html_convert(str(tmp_path / "does_not_exist.md"))
    assert msg.startswith("Error:")
    assert "not found" in msg


def test_non_md_extension_rejected(tmp_path):
    """Non-.md extension -> 'Error: not a .md file' string."""
    src = _write(tmp_path / "bogus.txt", "not markdown\n")
    msg = md2html_convert(src)
    assert msg.startswith("Error:")
    assert ".md" in msg


def test_md_to_html_returns_full_dict(tmp_path):
    """The lower-level md_to_html() returns dict with required keys."""
    src = _write(tmp_path / "dict.md", "# x\n\n$y=1$\n")
    out = md_to_html(src)
    for k in ("output_file", "math_blocks", "images_embedded", "log"):
        assert k in out, f"md_to_html missing key {k}"
    assert out["math_blocks"] == 1
    assert out["images_embedded"] == 0


def test_title_explicit_override(tmp_path):
    src = _write(tmp_path / "t.md", "# Auto Title\n\nbody\n")
    md2html_convert(src, title="Override Title")
    html = (tmp_path / "t.html").read_text(encoding="utf-8-sig")
    assert "<title>Override Title</title>" in html
    # The H1 still says the original
    assert "<h1>Auto Title</h1>" in html
