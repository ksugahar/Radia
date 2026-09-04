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


def test_source_surfaces_share_a_cjk_safe_font_stack() -> None:
    assert (
        "--eqed-source-font: ui-monospace, 'Cascadia Mono', "
        "'Yu Gothic UI', Meiryo, monospace"
    ) in SOURCE
    assert re.search(
        r"\.eqed-source\s*\{[^}]*font-family:\s*var\(--eqed-source-font\)",
        SOURCE,
    )
    assert re.search(
        r"\.eqed-recent\s*\{[^}]*font-family:\s*var\(--eqed-source-font\)",
        SOURCE,
    )
    assert "font-family: Consolas" not in SOURCE


def test_palette_and_learning_contract() -> None:
    for label in ["基本", "解析", "集合・記号", "幾何", "ギリシャ"]:
        assert f'label: "{label}"' in SOURCE
    assert "showRecentInsertion(snippet)" in SOURCE
    assert 'event.key !== "Tab"' in SOURCE
    assert "nextHole(input.value" in SOURCE


def test_structural_row_break_and_undo_contract() -> None:
    """Enter is the native row break, and every edit stays undoable."""
    assert 'event.key !== "Enter" || event.shiftKey' in SOURCE
    assert "event.isComposing || event.keyCode === 229" in SOURCE
    assert "composeRowBreak(" in SOURCE
    assert "function prettyTex(raw)" in SOURCE
    assert r'snippet.indexOf("\\begin{") >= 0) snippet = prettyTex(snippet)' in SOURCE
    # The palette insert and the row break go through the same undoable edit.
    assert "applyEdit(input, edit)" in SOURCE
    assert "input.value = edit.value" not in SOURCE
    assert 'document.execCommand(text ? "insertText" : "delete", false, text)' in SOURCE
    assert "input.setRangeText(text, start, end, \"end\")" in SOURCE
    # Row environments get {} cells so Tab reaches every one of them.
    for template in [r"\\begin{pmatrix} {} & {} \\\\ {} & {} \\end{pmatrix}",
                     r"\\begin{cases} {} & {} \\\\ {} & {} \\end{cases}"]:
        assert template in SOURCE
    assert "ROW_ENVIRONMENT" in SOURCE


def test_prime_never_creates_a_double_superscript() -> None:
    """A prime is a superscript, so it needs a group after another one.

    The structural editor attaches the prime to its base; a source pane
    inserts at the caret, so a palette prime pressed after `a^{2}` produced
    `a^{2}'`, which MathJax rejects with "Prime causes double exponent".
    """
    assert "function endsWithSuperscript(text)" in SOURCE
    assert "/^'+$/.test(snippet)" in SOURCE
    assert 'endsWithSuperscript(before) ? "{}" : ""' in SOURCE


def test_palette_keys_explain_themselves_without_hover() -> None:
    """A two-letter face such as "sf" or "tt" says nothing on its own.

    The native status bar describes the cell under the cursor.  A title
    tooltip cannot do that on a touch screen or for a keyboard user, so the
    Web editor mirrors the native status line on hover and on focus.
    """
    assert "function showKeyHelp(text)" in SOURCE
    assert "function clearKeyHelp(text)" in SOURCE
    assert 'button.setAttribute("aria-label", description)' in SOURCE
    assert 'button.addEventListener("mouseenter"' in SOURCE
    assert 'button.addEventListener("focus"' in SOURCE
    assert 'status.setAttribute("data-tex-literal-ok", "true")' in SOURCE


def test_teaching_surfaces_are_exempt_from_the_bare_tex_scan() -> None:
    """The hint line and recent-insertion display show TeX on purpose."""
    assert 'recent.setAttribute("data-tex-literal-ok", "true")' in SOURCE
    assert 'hint.setAttribute("data-tex-literal-ok", "true")' in SOURCE
    assert '"eqed-hint"' in SOURCE
    assert "Shift+Enter" in SOURCE


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
    assert "officeOmml(mml)" not in office
    assert "writeOfficeClipboard(html, tex)" in office
    assert "<!--[if gte msEquation 12]>" not in office
    assert "var html = mml + '<span style=\"font-size:18pt\">&#160;</span>'" in office
    assert '<m:oMath' not in office
    assert '<m:oMathPara' not in office
    assert "<!DOCTYPE html>" not in office


def test_office_mathml_is_canonical_inline_18pt() -> None:
    canonical = SOURCE.split("function officeMathMl", 1)[1].split(
        "function getSvgConverter", 1)[0]
    assert "MathJax.tex2mml" in canonical
    assert '{ display: false }' in canonical
    assert 'setAttribute("display", "inline")' in canonical
    assert 'setAttribute("mathsize", "18pt")' in canonical
    assert "splitUnanchoredOfficeRows(math)" in canonical
    assert 'cells[0].textContent.trim() === ""' in SOURCE
    assert "onlyNestedSyntheticTable(content)" in SOURCE
    assert "collectLeafCells(nested, output)" in SOURCE
    assert "leafCells.length < 2" in SOURCE
    assert "officeRows.join(" in canonical
    assert 'font-size:18pt' in canonical
    assert "maligngroup" not in canonical
    assert "malignmark" not in canonical
    assert 'querySelectorAll("mstyle")' in canonical
    assert 'attribute.name.indexOf("data-") === 0' in canonical
    assert 'setAttribute("largeop", "true")' in canonical
    assert 'node.textContent === "―"' in canonical
    assert 'parent === "mover" ? "¯" : "_"' in canonical
    assert 'node.setAttribute("stretchy", "true")' in canonical


def test_autoloaded_macros_are_warmed_before_the_first_office_copy() -> None:
    """The Office copy is synchronous, so no package may still be in flight.

    MathJax loads \\boldsymbol and \\cancel on demand and makes the
    synchronous tex2mml throw "MathJax retry" until the package arrives, which
    failed the first copy of a `\\bm` equation on a freshly opened page.
    """
    assert "function warmAutoloadedMacros()" in SOURCE
    assert r'tex2mmlPromise("\\boldsymbol{x}+\\cancel{x}"' in SOURCE
    assert "window.MathJax.startup.promise.then(warmAutoloadedMacros)" in SOURCE


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
