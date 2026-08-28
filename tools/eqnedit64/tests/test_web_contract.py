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


def test_office_copy_is_editable_mathml_without_png_competition() -> None:
    office = SOURCE.split(
        'root.querySelector(".eqed-copy-office")', 1
    )[1].split(
        'root.querySelector(".eqed-copy-display")', 1
    )[0]
    assert "MathJax.tex2mml" in office
    assert '{ display: false }' in office
    assert '"text/html"' in office
    assert '"text/plain"' in office
    assert '"image/png"' not in office
    assert "&#160;</body></html>" in office


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
