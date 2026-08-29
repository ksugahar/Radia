"""Repository-level contract checks for the browser equation editor."""
from __future__ import annotations

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (ROOT / "web" / "equation-editor.js").read_text(encoding="utf-8")
FRAGMENT = (ROOT / "web" / "equation-editor.fragment.html").read_text(
    encoding="utf-8"
)
WEB_README = (ROOT / "web" / "README.md").read_text(encoding="utf-8")


def test_web_source_is_radia_owned_tex_only() -> None:
    lowered = (SOURCE + FRAGMENT).lower()
    assert "mtef" not in lowered
    assert ".eqn" not in lowered
    assert 'class="eqed-source"' in FRAGMENT
    assert "TeXソース" in FRAGMENT


def test_palette_and_learning_contract() -> None:
    for label in ["基本", "解析", "集合・記号", "幾何", "ギリシャ"]:
        assert f'label: "{label}"' in SOURCE
    assert "showRecentInsertion(snippet)" in SOURCE
    assert 'event.key !== "Tab"' in SOURCE
    assert "nextHole(input.value" in SOURCE


def test_common_math_alphabets_are_always_visible() -> None:
    assert "var MATH_ALPHABETS" in SOURCE
    for command in [r"\\mathrm{}", r"\\mathit{}", r"\\mathbf{}"]:
        assert f'snippet: "{command}"' in SOURCE
    common = SOURCE.split("var MATH_ALPHABETS", 1)[1].split("];", 1)[0]
    assert 'name: "立体"' in common
    assert 'name: "変数（斜体）"' in common
    assert 'name: "ベクトル"' in common
    assert "paletteHead.appendChild(styleBar)" in SOURCE
    assert "paletteHost.appendChild(paletteHead)" in SOURCE
    assert 'button.addEventListener("click", function () { insert(alphabet.snippet); })' in SOURCE


def test_extended_math_alphabets_remain_in_decoration_palette() -> None:
    for command in [r"\\mathsf{}", r"\\mathtt{}", r"\\mathcal{}",
                    r"\\mathbb{}", r"\\mathfrak{}",
                    r"\\bm{}", r"\\mathnormal{}"]:
        assert command in SOURCE
    assert 'function mathJaxTex(tex)' in SOURCE
    assert 'doc.convert(mathJaxTex(tex)' in SOURCE
    assert '"\\\\[" + mathJaxTex(tex) + "\\\\]"' in SOURCE
    assert '"\\\\displaystyle " + mathJaxTex(tex)' in SOURCE


def test_office_copy_is_editable_mathml_without_png_competition() -> None:
    office = SOURCE.split(
        'root.querySelector(".eqed-copy-office")', 1
    )[1].split(
        'root.querySelector(".eqed-copy-display")', 1
    )[0]
    assert "officeMathMl(tex)" in office
    assert '"image/png"' not in office
    assert "officeOmml(mml)" in office
    assert "writeOfficeClipboard(html, tex)" in office
    assert "<!--[if gte msEquation 12]>" in office
    assert "<![if !msEquation]>" in office
    assert 'text-align:left;font-size:24pt' in office
    assert '<m:oMathParaPr><m:jc m:val=\\\"left\\\"/>' in office
    assert "<!DOCTYPE html>" not in office


def test_office_mathml_is_canonical_inline_24pt() -> None:
    canonical = SOURCE.split("function officeMathMl", 1)[1].split(
        "function getSvgConverter", 1)[0]
    assert "MathJax.tex2mml" in canonical
    assert '{ display: false }' in canonical
    assert 'setAttribute("display", "inline")' in canonical
    assert 'setAttribute("mathsize", "24pt")' in canonical
    assert 'querySelectorAll("mstyle")' in canonical
    assert 'attribute.name.indexOf("data-") === 0' in canonical
    assert 'setAttribute("largeop", "true")' in canonical


def test_office_omml_transport_covers_structured_math() -> None:
    for tag in ["m:f", "m:rad", "m:sSup", "m:sSub", "m:sSubSup",
                "m:nary", "m:acc", "m:bar", "m:d", "m:m", "m:borderBox"]:
        assert f"<{tag}>" in SOURCE
    assert "font-size:24.0pt" in SOURCE
    assert "Cambria Math" in SOURCE


def test_office_copy_prefers_exact_cf_html_fragment() -> None:
    transport = SOURCE.split("function writeOfficeClipboard", 1)[1].split(
        "function getSvgConverter", 1
    )[0]
    assert 'document.addEventListener("copy", onCopy, true)' in transport
    assert 'event.clipboardData.setData("text/html", html)' in transport
    assert 'event.clipboardData.setData("text/plain", tex)' in transport
    assert '"text/html": new Blob([html]' in transport
    assert '"text/plain": new Blob([tex]' in transport
    assert 'document.execCommand("copy")' in transport
    assert 'return Promise.resolve("copy-event")' in transport
    assert "window.ClipboardItem" in transport


def test_png_is_a_separate_user_action() -> None:
    assert 'root.querySelector(".eqed-copy-png")' in SOURCE
    assert 'new window.ClipboardItem({ "image/png"' in SOURCE
    assert re.search(r"var\s+PNG_SCALE\s*=\s*\d+\s*;", SOURCE)


def test_fragment_matches_script_mount_contract() -> None:
    required = {
        "eqed-palettes",
        "eqed-source",
        "eqed-preview",
        "eqed-actions",
        "eqed-copy-office",
        "eqed-copy-png",
        "eqed-save-svg",
        "eqed-copy-display",
        "eqed-copy-equation",
        "eqed-clear",
        "eqed-status",
    }
    for class_name in required:
        assert class_name in FRAGMENT
        assert class_name in SOURCE
    assert "data-equation-editor" in FRAGMENT


def test_homepage_publication_imports_radia_source() -> None:
    normalized = " ".join(WEB_README.split())
    assert "RADIA_REPOSITORY" in WEB_README
    assert "SHA-256" in WEB_README
    assert "Do not retain or edit an independent homepage source copy" in normalized
    assert "run_eqnedit64_release_qa.ps1" in WEB_README
    assert "does not require the Mathematica, NGSolve" in normalized
