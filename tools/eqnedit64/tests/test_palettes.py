# -*- coding: utf-8 -*-
"""The toolbar palettes must be a complete, non-overlapping catalogue.

Eqnedt32 reaches every one of its symbols and templates from a toolbar
palette, so nothing is available only to someone who already knows the
command.  Eqnedit64 drifted away from that: at one point 117 of its 160
symbols had neither a shortcut nor a mouse route.

These are rules, not examples.  Adding a symbol to the table or a kind to
insert_template() with no palette home fails here, which is what stops the
gap from reopening one entry at a time.
"""
import json
import os
from pathlib import Path
import re
import sys
import unicodedata

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "build"))
import eqnedit_core as E  # noqa: E402


def palette_items():
    for title, face, columns, items in E.palettes():
        for command, cell, label in items:
            yield title, command, cell, label


def test_every_symbol_has_a_palette_home():
    placed = {c[len("symbol."):] for _, c, _, _ in palette_items()
              if c.startswith("symbol.")}
    missing = sorted(set(E.symbol_commands()) - placed)
    assert not missing, (
        "%d symbol(s) reachable only by typing the command: %s"
        % (len(missing), " ".join(missing)))


def test_every_template_has_a_palette_home():
    placed = {c[len("template."):] for _, c, _, _ in palette_items()
              if c.startswith("template.")}
    missing = sorted(set(E.Equation.templates()) - placed)
    assert not missing, (
        "%d template(s) with no button: %s" % (len(missing), " ".join(missing)))


def test_no_command_sits_on_two_palettes():
    seen, twice = set(), []
    for _, command, _, _ in palette_items():
        if command in seen:
            twice.append(command)
        seen.add(command)
    assert not twice, "on more than one palette: %s" % " ".join(sorted(set(twice)))


def test_no_palette_names_something_that_does_not_exist():
    symbols = set(E.symbol_commands())
    templates = set(E.Equation.templates())
    matrix_actions = {
        "matrix.add_row", "matrix.remove_row",
        "matrix.add_column", "matrix.remove_column",
    }
    unknown = []
    for _, command, _, _ in palette_items():
        if command.startswith("symbol.") and command[7:] not in symbols:
            unknown.append(command)
        elif command.startswith("template.") and command[9:] not in templates:
            unknown.append(command)
        elif command.startswith("matrix.") and command not in matrix_actions:
            unknown.append(command)
    assert not unknown, "no such command: %s" % " ".join(unknown)


def test_every_palette_cell_is_labelled():
    for title, command, cell, label in palette_items():
        assert cell, "%s: %s has no face" % (title, command)
        assert label, "%s: %s has no description" % (title, command)


def test_every_palette_cell_produces_a_tex_lesson():
    """The shortcut coach derives its lesson by executing each real command.

    A parser hardening change once made the supported ``\\flat`` and
    ``\\sharp`` geometry entries disappear as if they were unknown commands.
    Catalogue completeness alone could not see that regression: the cells
    still existed, but their generated lesson and insertion were empty.
    """
    empty = []
    for title, command, _, _ in palette_items():
        equation = E.Equation()
        if command.startswith("template."):
            changed = equation.insert_template(command[9:])
        elif command.startswith("symbol."):
            changed = equation.insert_symbol(command[7:])
        elif command.startswith("latex."):
            changed = equation.insert_latex(command[6:])
        elif command.startswith("style."):
            equation.insert_text("x")
            equation.select_all()
            changed = equation.restyle_selection(command[6:])
        elif command.startswith("matrix."):
            changed = (equation.insert_template("matrix2x2") and
                       equation.command(command))
        else:
            changed = False
        if not changed or not equation.latex():
            empty.append(f"{title}: {command}")
    assert not empty, "palette cells with no TeX lesson: " + ", ".join(empty)


def test_symbol_palettes_come_first_then_template_palettes():
    """The shape of Eqnedt32's bar: symbols on the first row, templates on
    the second, so the eye knows which row to look at.  Eqnedit64 adds a
    browser-shared differential-geometry symbol palette.  The dedicated
    over/underline palette remains on the template row, while labelled arrows
    (\\xrightarrow and friends) have no node type yet."""
    entries = E.palettes()
    split = E.symbol_palette_count()
    assert split == 11
    assert len(entries) > split, "no template palettes"
    for title, _, _, items in entries[:split]:
        kinds = {c.split(".")[0] for c, _, _ in items}
        assert kinds <= {"symbol", "latex", "template", "style"}, title
    for title, _, _, items in entries[split:]:
        kinds = {c.split(".")[0] for c, _, _ in items}
        assert kinds <= {"template", "matrix"} and "template" in kinds, (
            "%s is on the template row but holds %s" % (title, kinds))


def test_matrix_palette_offers_rectangles_and_structural_resize():
    commands = {
        command for title, command, _, _ in palette_items() if title == "行列"
    }
    assert {
        "template.matrix1x2", "template.matrix2x1",
        "template.matrix2x3", "template.matrix3x2",
        "template.matrix4x4", "template.matrix5x5", "template.matrix6x6",
        "matrix.add_row", "matrix.remove_row",
        "matrix.add_column", "matrix.remove_column",
    } <= commands


def test_palette_faces_use_unambiguous_owned_glyphs():
    by_title = {
        title: {command: cell for command, cell, _ in items}
        for title, _, _, items in E.palettes()
    }
    scripts = by_title["上下付き"]
    for command in ("template.sup", "template.sub", "template.subsup"):
        assert not any(character.isdigit() for character in scripts[command]), (
            f"{command} uses a real digit as an empty-slot marker: "
            f"{scripts[command]!r}")

    all_faces = [palette_face for _, palette_face, _, _ in E.palettes()]
    all_faces.extend(
        face
        for _, _, _, items in E.palettes()
        for _, face, _ in items)
    assert all("▯" not in face for face in all_faces), (
        "U+25AF is absent from the embedded Latin Modern Math cmap; "
        "use the owned U+25A1 empty slot instead")
    assert any("□" in face for face in all_faces)
    assert not [
        (face, character)
        for face in all_faces
        for character in face
        if unicodedata.combining(character)
    ], "GDI owner-draw palette faces must not rely on combining-mark shaping"

    matrix = by_title["行列"]
    assert matrix["matrix.add_row"] == "+R"
    assert matrix["matrix.remove_row"] == "−R"
    assert matrix["matrix.add_column"] == "+C"
    assert matrix["matrix.remove_column"] == "−C"


def test_five_categories_cover_every_palette_once():
    """The native tabs use the same five-way vocabulary as the web editor.

    Keeping this map in the tested catalogue prevents a menu-only refactor
    from hiding a palette or showing it in two places.
    """
    categories = E.palette_categories()
    assert [title for title, _ in categories] == [
        "基本", "解析", "集合・記号", "幾何", "ギリシャ"]
    indices = [index for _, members in categories for index in members]
    assert sorted(indices) == list(range(len(E.palettes())))
    assert len(indices) == len(set(indices)), "a palette appears on two tabs"
    assert all(members for _, members in categories), "an empty tab is not useful"
    assert max(len(members) for _, members in categories) <= 6


def test_categories_follow_the_web_learning_order():
    """The web editor starts with structures, then brackets, then decoration.

    The native catalogue is larger, but the first left-to-right choices should
    teach the same path instead of putting decoration before basic structure.
    """
    palettes = E.palettes()
    by_category = {
        title: [palettes[index][0] for index in members]
        for title, members in E.palette_categories()
    }
    assert by_category["基本"] == [
        "分数と根号", "上下付き", "行列", "括弧", "装飾", "上線と下線"]
    assert by_category["集合・記号"] == [
        "矢印", "集合記号", "論理記号", "総乗と集合演算", "空白と点"]
    assert by_category["幾何"] == ["微分幾何"]


def test_native_geometry_palette_matches_the_web_editor():
    """Geometry is a shared product contract, not a spare category.

    It once pointed at the unrelated over/underline palette while the Web
    product taught differential geometry.  Keep faces, insertion TeX, labels,
    and order identical at the product boundary.
    """
    expected = [
        ("∧", r"\wedge ", "ウェッジ積"),
        ("★", r"\star ", "Hodge star作用素"),
        ("dω", r"\mathrm{d} ", "外微分（立体のd）"),
        ("ι", r"\iota_{} ", "内部積（縮約）ι_X"),
        ("ℒ", r"\mathcal{L}_{} ", "Lie微分"),
        ("f^*", "^{*} ", "引き戻し（pullback）"),
        ("f_*", "_{*} ", "押し出し（pushforward）"),
        ("♭", r"^{\flat} ", "フラット（添字を下げる）"),
        ("♯", r"^{\sharp} ", "シャープ（添字を上げる）"),
        ("⊗", r"\otimes ", "テンソル積"),
        ("⊕", r"\oplus ", "直和"),
    ]
    palette = next(items for title, _, _, items in E.palettes()
                   if title == "微分幾何")
    actual = []
    for command, face, label in palette:
        if command.startswith("symbol."):
            payload = command[len("symbol."):]
        else:
            assert command.startswith("latex.")
            payload = command[len("latex."):]
        actual.append((face, payload, label))
    assert actual == expected

    source = (Path(__file__).resolve().parents[1] / "web" /
              "equation-editor.js").read_text(encoding="utf-8")
    block = source.split('label: "微分幾何"', 1)[1].split("]\n    },", 1)[0]
    web = []
    for match in re.finditer(r'^\s*(\[".*"\]),?$', block, re.MULTILINE):
        face, payload, label = json.loads(match.group(1))
        web.append((face, payload, label))
    assert web == expected


def test_every_palette_command_round_trips():
    """A button that inserts something the parser cannot read back is worse
    than no button: the equation survives until it is saved and reopened."""
    broken = []
    for _, command, _, _ in palette_items():
        equation = E.Equation()
        equation.load_latex("")
        if command.startswith("template."):
            kind = command[9:]
            assert equation.insert_template(kind), command
            equation.insert_text("x")
            if kind.startswith("matrix"):
                rows, columns = equation.matrix_dimensions()
                for _ in range(rows * columns - 1):
                    assert equation.next_slot(), command
                    equation.insert_text("x")
        elif command.startswith("symbol."):
            assert equation.insert_symbol(command[7:]), command
        elif command.startswith("matrix."):
            assert equation.insert_template("matrix2x2"), command
            for position, value in enumerate("1234"):
                equation.insert_text(value)
                if position != 3:
                    equation.next_slot()
            assert equation.command(command), command
        elif command.startswith("style."):
            equation.insert_text("x")
            equation.select_all()
            assert equation.restyle_selection(command[6:]), command
        else:
            equation.insert_latex(command[6:])
            equation.insert_text("x")
        tex = equation.latex()
        once = E.tex_normalize(tex)
        twice = E.tex_normalize(once)
        if once != twice or once != tex:
            broken.append("%s: %r -> %r -> %r" % (command, tex, once, twice))
        if not tex.isascii():
            broken.append("%s: not ASCII: %r" % (command, tex))
    assert not broken, "\n".join(broken)


if __name__ == "__main__":
    failures = 0
    for name, function in sorted(globals().items()):
        if not name.startswith("test_"):
            continue
        try:
            function()
            print("[OK]   %s" % name)
        except AssertionError as error:
            failures += 1
            print("[FAIL] %s\n       %s" % (name, error))
    total = sum(len(items) for _, _, _, items in E.palettes())
    print("\n%d palettes, %d cells, %d failure(s)"
          % (len(E.palettes()), total, failures))
    sys.exit(1 if failures else 0)
