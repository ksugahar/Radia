"""The editing model, exercised the way a person types.

Every case is a keystroke sequence, not an API call sequence: what is being
checked is that Equation Editor's *feel* survived -- Ctrl+L after `x` gives `x`
a subscript rather than an empty box, Tab walks the holes of a template,
backspace at the start of an empty slot unwraps the template instead of
swallowing it, and undo puts everything back.

Run:  python tests\\test_edit.py
"""
from __future__ import annotations

import re
import os
import sys
import xml.etree.ElementTree as ET

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "build"))

from eqnedit_core import (  # noqa: E402
    Equation, MAX_NESTING_DEPTH, SvgStyle, tex_to_mathml,
    tex_to_office_mathml_fragment, tex_to_svg,
)

# Each case: (name, keystrokes, expected LaTeX).
# A keystroke is a literal string to type, or a command name to dispatch.
CASES = [
    ("plain text",
     ["x", "+", "y"],
     "x+y"),

    ("Ctrl+L makes the typed letter the base",
     ["x", "template.sub", "i"],
     "x_{i}"),

    ("Ctrl+H after a letter",
     ["a", "template.sup", "2"],
     "a^{2}"),

    ("Ctrl+J gives both scripts, Tab moves between them",
     ["x", "template.subsup", "i", "caret.next_slot", "2"],
     "x_{i}^{2}"),

    ("fraction: type numerator, Tab, type denominator",
     ["template.frac", "a", "caret.next_slot", "b"],
     "\\frac{a}{b}"),

    ("fraction then leave it and keep typing",
     ["template.frac", "1", "caret.next_slot", "2",
      "caret.next_slot", "+", "3"],
     "\\frac{1}{2}+3"),

    ("square root",
     ["template.sqrt", "2"],
     "\\sqrt{2}"),

    ("vector accent wraps the previous item",
     ["x", "template.vec"],
     "\\vec{x}"),

    ("hat template serializes structurally",
     ["template.hat", "x"],
     "\\hat{x}"),

    ("nested: a fraction inside a square root",
     ["template.sqrt", "template.frac", "a", "caret.next_slot", "b"],
     "\\sqrt{\\frac{a}{b}}"),

    ("parentheses",
     ["template.paren", "x", "+", "1"],
     "\\left( x+1 \\right)"),

    # LaTeX has no place to record "limits are stacked" on a \sum -- it stacks
    # them in display style anyway -- so that flag does not survive the trip.
    ("summation with both limits",
     ["template.sum", "i", "caret.next_slot", "n", "caret.next_slot", "a"],
     "\\sum _{i}^{n} a"),

    ("Greek by name",
     ["\\alpha", "+", "\\beta"],
     "\\alpha +\\beta "),

    ("backspace deletes the character before the caret",
     ["a", "b", "edit.backspace"],
     "a"),

    ("backspace deletes a populated superscript character",
     ["E", "=", "m", "c", "template.sup", "2", "caret.next_slot",
      "H", "template.sup", "2", "edit.backspace"],
     "E = mc^{2}H^{}"),

    ("backspace in an empty slot unwraps the template, keeping the content",
     ["template.frac", "a", "caret.next_slot", "caret.prev_slot",
      "caret.home", "edit.backspace"],
     "a"),

    ("undo restores what the last edit changed",
     ["a", "template.frac", "b", "edit.undo", "edit.undo"],
     "a"),

    ("redo after undo",
     ["a", "b", "edit.undo", "edit.redo"],
     "ab"),
]


def type_in(eq: "Equation", keys) -> None:
    for k in keys:
        if k.startswith("\\"):
            if not eq.insert_symbol(k):
                raise AssertionError(f"unknown symbol {k}")
        elif "." in k and not k.isspace():
            if not eq.command(k):
                raise AssertionError(f"command {k} did nothing")
        else:
            eq.insert_text(k)


def main() -> int:
    failures = []
    selection_contract_checks = 0

    for name, keys, expect in CASES:
        eq = Equation()
        try:
            type_in(eq, keys)
        except AssertionError as exc:
            failures.append(f"{name}: {exc}")
            continue
        got = eq.latex()
        if got.strip() != expect.strip():
            failures.append(f"{name}:\n    want {expect!r}\n    got  {got!r}")

    # The caret has to stay inside the tree no matter how it is driven.
    eq = Equation()
    type_in(eq, ["template.frac", "a", "caret.next_slot", "b"])
    for _ in range(40):
        eq.move_left()
    for _ in range(40):
        eq.move_right()
    if eq.latex().strip() != "\\frac{a}{b}":
        failures.append(f"walking the caret changed the equation: {eq.latex()!r}")

    # Mouse hit-testing lands on structural insertion sites, not character
    # offsets in a hidden LaTeX string.
    eq = Equation()
    eq.insert_text("ab")
    width, _height, baseline = eq.metrics()
    if not eq.hit_test(0.0, baseline) or eq.caret() != ":0":
        failures.append(f"left-edge hit did not reach the first caret: {eq.caret()}")
    if not eq.hit_test(width, baseline) or eq.caret() != ":2":
        failures.append(f"right-edge hit did not reach the last caret: {eq.caret()}")
    caret_geometry = eq.caret_geometry()
    if (caret_geometry is None or
            not (0 <= caret_geometry[0] <= width) or
            not (caret_geometry[1] < caret_geometry[2])):
        failures.append(
            f"caret geometry was invalid for IME placement: {caret_geometry!r}")

    def caret_point(equation):
        geometry = equation.caret_geometry()
        if geometry is None:
            return None
        return geometry[0], 0.5 * (geometry[1] + geometry[2])

    def compact(tex):
        return tex.replace(" ", "").replace("\r", "").replace("\n", "")

    # Cross-slot pointer selection is a tree rule, not a special case for
    # fractions.  Exercise every editable multi-slot container class through
    # the same public begin/extend/end lifecycle.
    multi_slot_templates = (
        ("frac", 2), ("nthroot", 2),
        ("sub", 2), ("sup", 2), ("subsup", 3),
        ("dirac", 2), ("int", 3), ("sum", 3),
        ("over", 2), ("under", 2),
        ("overbrace", 2), ("underbrace", 2),
        ("matrix2x2", 4), ("cases", 4),
    )
    for template, slot_count in multi_slot_templates:
        structured = Equation()
        if not structured.insert_template(template):
            failures.append(f"{template}: could not build cross-slot fixture")
            continue
        built = True
        for slot_number in range(slot_count):
            structured.insert_text(chr(ord("a") + slot_number))
            if slot_number + 1 < slot_count and not structured.next_slot():
                failures.append(f"{template}: could not reach slot {slot_number + 1}")
                built = False
                break
        if not built:
            continue
        points = []
        for slot_number in range(slot_count - 1, -1, -1):
            structured.move_end()
            point = caret_point(structured)
            if point is None:
                failures.append(f"{template}: slot {slot_number} has no caret geometry")
                break
            points.append(point)
            if slot_number > 0 and not structured.prev_slot():
                failures.append(f"{template}: could not return to slot {slot_number - 1}")
                break
        points.reverse()
        if len(points) != slot_count:
            continue
        whole = structured.latex()
        began = structured.begin_pointer_selection(*points[0])
        extended = structured.extend_pointer_selection(*points[-1])
        selected = structured.selection_latex()
        structured.end_pointer_selection()
        if not began or not extended or compact(selected) != compact(whole):
            failures.append(
                f"{template}: cross-slot drag selected {selected!r}, want {whole!r}")
        selection_contract_checks += 1

    # PileNode is produced by gathered rather than by a palette template.
    gathered = Equation()
    gathered.load_latex(r"\begin{gathered}a\\b\end{gathered}")
    gathered.move_home()
    if not gathered.move_right():
        failures.append("gathered: could not enter first row")
    else:
        gathered.move_end()
        first_row = caret_point(gathered)
        if not gathered.next_slot():
            failures.append("gathered: could not enter second row")
        else:
            gathered.move_end()
            second_row = caret_point(gathered)
            whole = gathered.latex()
            if (first_row is None or second_row is None or
                    not gathered.begin_pointer_selection(*first_row) or
                    not gathered.extend_pointer_selection(*second_row) or
                    compact(gathered.selection_latex()) != compact(whole)):
                failures.append("gathered: cross-row drag did not select the pile")
            gathered.end_pointer_selection()
    selection_contract_checks += 1

    # Single-slot containers are crossed by dragging from their content to a
    # neighbour in the parent slot.  This covers radical, fence, limit,
    # stretch decoration, and character embellishment nodes.
    for template in ("sqrt", "paren", "lim", "overline", "hat"):
        structured = Equation()
        structured.insert_template(template)
        structured.insert_text("x")
        inside = caret_point(structured)
        structured.move_out()
        structured.insert_text("z")
        outside = caret_point(structured)
        whole = structured.latex()
        if (inside is None or outside is None or
                not structured.begin_pointer_selection(*inside) or
                not structured.extend_pointer_selection(*outside) or
                compact(structured.selection_latex()) != compact(whole)):
            failures.append(
                f"{template}: inside-to-parent drag did not include the container")
        structured.end_pointer_selection()
        selection_contract_checks += 1

    # Parser-owned grouping/math-alphabet containers obey the same rule.
    for name, tex in (("group", "{x}z"), ("math alphabet", r"\mathbf{x}z")):
        structured = Equation()
        structured.load_latex(tex)
        structured.move_home()
        entered = structured.move_right()
        structured.move_end()
        inside = caret_point(structured)
        structured.move_out()
        structured.move_end()
        outside = caret_point(structured)
        whole = structured.latex()
        if (not entered or inside is None or outside is None or
                not structured.begin_pointer_selection(*inside) or
                not structured.extend_pointer_selection(*outside) or
                compact(structured.selection_latex()) != compact(whole)):
            failures.append(f"{name}: inside-to-parent drag lost its container")
        structured.end_pointer_selection()
        selection_contract_checks += 1

    def promoted_fraction_selection():
        equation = Equation()
        equation.load_latex(r"a\frac{x}{y}b")
        equation.move_home()
        equation.move_right()  # after a
        equation.move_right()  # numerator slot
        equation.move_end()
        inside = caret_point(equation)
        equation.move_out()
        equation.move_end()
        outside = caret_point(equation)
        if (inside is None or outside is None or
                not equation.begin_pointer_selection(*inside) or
                not equation.extend_pointer_selection(*outside)):
            return equation, None, equation.latex()
        selected = equation.selection_latex()
        original = equation.latex()
        equation.end_pointer_selection()
        return equation, selected, original

    # Every consumer must use exactly the promoted structural range.  This
    # catches a future fix that updates painting/copy but leaves deletion,
    # replacement, wrapping, or recursive styling on the old shallow range.
    probe, selected, _ = promoted_fraction_selection()
    if selected is None or compact(selected) != r"\frac{x}{y}b":
        failures.append(f"promoted copy range was {selected!r}")
    selection_contract_checks += 1

    for name, operation in (
            ("delete", lambda e: e.delete_selection()),
            ("backspace", lambda e: e.backspace()),
            ("forward delete", lambda e: e.erase())):
        probe, _selected, _original = promoted_fraction_selection()
        if not operation(probe) or compact(probe.latex()) != "a":
            failures.append(f"promoted {name} left {probe.latex()!r}")
        selection_contract_checks += 1

    # Promotion stops at the child reached by the pointer.  An adjacent
    # operator which the drag did not reach is deliberately not absorbed;
    # Word and MathType use the same structural boundary.
    probe = Equation()
    probe.load_latex(r"a\frac{x}{y}+")
    probe.move_home()
    probe.move_right()  # after a
    probe.move_right()  # numerator slot
    probe.move_end()
    inside = caret_point(probe)
    probe.move_out()    # root slot immediately after the fraction
    after_fraction = caret_point(probe)
    if (inside is None or after_fraction is None or
            not probe.begin_pointer_selection(*inside) or
            not probe.extend_pointer_selection(*after_fraction)):
        failures.append("adjacent-operator selection fixture could not drag")
    else:
        selected = compact(probe.selection_latex())
        probe.end_pointer_selection()
        if selected != r"\frac{x}{y}":
            failures.append(
                f"adjacent-operator selection promoted to {selected!r}")
        elif not probe.delete_selection() or compact(probe.latex()) != "a+":
            failures.append(
                f"promoted deletion absorbed an adjacent operator: {probe.latex()!r}")
    selection_contract_checks += 1

    probe, _selected, _original = promoted_fraction_selection()
    probe.insert_text("q")
    if compact(probe.latex()) != "aq":
        failures.append(f"promoted typing replacement left {probe.latex()!r}")
    selection_contract_checks += 1

    probe, _selected, _original = promoted_fraction_selection()
    if not probe.insert_latex(r"\sqrt{q}") or compact(probe.latex()) != r"a\sqrt{q}":
        failures.append(f"promoted TeX replacement left {probe.latex()!r}")
    selection_contract_checks += 1

    probe, _selected, _original = promoted_fraction_selection()
    if not probe.insert_symbol(r"\alpha") or compact(probe.latex()) != r"a\alpha":
        failures.append(f"promoted symbol replacement left {probe.latex()!r}")
    selection_contract_checks += 1

    probe, _selected, _original = promoted_fraction_selection()
    if (not probe.insert_template("paren") or
            compact(probe.latex()) != r"a\left(\frac{x}{y}b\right)"):
        failures.append(f"promoted template wrapping left {probe.latex()!r}")
    selection_contract_checks += 1

    probe, _selected, original = promoted_fraction_selection()
    if not probe.restyle_selection("vector"):
        failures.append("promoted recursive style change was rejected")
    else:
        styled = probe.latex()
        for letter in ("x", "y", "b"):
            if rf"\mathbf{{{letter}}}" not in styled:
                failures.append(
                    f"promoted style missed {letter!r} inside {styled!r}")
        if not probe.undo() or probe.latex() != original:
            failures.append("Undo did not restore a promoted style change")
        elif not probe.redo() or probe.latex() != styled:
            failures.append("Redo did not restore a promoted style change")
    selection_contract_checks += 1

    # Rebuilding the tree and completing a drag both invalidate its deep
    # origin.  An extend call without a live begin must fail rather than reuse
    # a path into an old tree.
    lifecycle = Equation()
    lifecycle.load_latex("xy")
    lifecycle.move_home()
    old_point = caret_point(lifecycle)
    if old_point is None or not lifecycle.begin_pointer_selection(*old_point):
        failures.append("pointer lifecycle fixture could not begin")
    lifecycle.load_latex("ab")
    if lifecycle.extend_pointer_selection(0.0, lifecycle.metrics()[2]):
        failures.append("document load retained a stale pointer origin")
    lifecycle.begin_pointer_selection(0.0, lifecycle.metrics()[2])
    lifecycle.end_pointer_selection()
    if lifecycle.extend_pointer_selection(0.0, lifecycle.metrics()[2]):
        failures.append("pointer end retained a reusable origin")
    selection_contract_checks += 2

    # A visual selection is structural: copy emits valid LaTeX, typing
    # replaces it, and a template wraps it instead of discarding it.
    eq.select_all()
    if eq.selection_latex().strip() != "ab":
        failures.append(f"selection copied the wrong structure: {eq.selection_latex()!r}")
    eq.insert_text("π")
    if "\\pi" not in eq.latex():
        failures.append(f"UTF-8 typing was split into bytes: {eq.latex()!r}")
    eq.load_latex("a+b")
    eq.select_all()
    eq.insert_template("frac")
    if eq.latex().replace(" ", "") != "\\frac{a+b}{}":
        failures.append(f"fraction did not wrap the selection: {eq.latex()!r}")

    # Eqnedit32 double-clicked the current structural slot, not the whole
    # equation.  The outermost slot still means the full equation.
    eq = Equation()
    eq.insert_template("frac")
    eq.insert_text("numerator")
    if not eq.select_current_slot() or eq.selection_latex() != "numerator":
        failures.append(
            "nested double-click selection escaped its slot: "
            f"{eq.selection_latex()!r}")
    eq.clear_selection()
    eq.next_slot()
    eq.insert_text("denominator")
    if not eq.select_current_slot() or eq.selection_latex() != "denominator":
        failures.append(
            "denominator slot selection was incorrect: "
            f"{eq.selection_latex()!r}")
    eq.load_latex("a+b")
    if not eq.select_current_slot() or eq.selection_latex().replace(" ", "") != "a+b":
        failures.append(
            "outer-slot double-click did not select the equation: "
            f"{eq.selection_latex()!r}")

    # Ctrl+click on a template selects the complete innermost structure in its
    # parent slot.  Cutting it must never strand one half of a fraction.
    eq = Equation()
    eq.insert_template("frac")
    eq.insert_text("a")
    eq.next_slot()
    eq.insert_text("b")
    if (not eq.select_containing_structure() or
            eq.selection_latex().replace(" ", "") != r"\frac{a}{b}"):
        failures.append(
            "containing-structure selection tore the fraction apart: "
            f"{eq.selection_latex()!r}")
    elif not eq.delete_selection() or eq.latex() != "":
        failures.append(
            "deleting a containing structure left an invalid remainder: "
            f"{eq.latex()!r}")

    for template, contents in (
            ("sqrt", "x"), ("paren", "x"), ("matrix2x2", "x")):
        eq = Equation()
        eq.insert_template(template)
        eq.insert_text(contents)
        complete = eq.latex()
        if (not eq.select_containing_structure() or
                eq.selection_latex() != complete):
            failures.append(
                f"{template} was not selected as one structure: "
                f"{eq.selection_latex()!r} vs {complete!r}")
        elif not eq.delete_selection() or eq.latex() != "":
            failures.append(
                f"deleting {template} left a remainder: {eq.latex()!r}")

    # Creating a script checkpoints before moving its base into the template.
    eq = Equation()
    eq.insert_text("x")
    eq.insert_template("sub")
    eq.undo()
    if eq.latex().strip() != "x":
        failures.append(f"undo lost the script base: {eq.latex()!r}")

    # The same action name drives both the Undo and Redo menu labels/logs.
    eq = Equation()
    eq.insert_text("x")
    if not eq.can_undo() or eq.undo_name() != "Typing" or eq.can_redo():
        failures.append(
            "typing history was not named/enabled correctly: "
            f"undo={eq.undo_name()!r} redo={eq.redo_name()!r}")
    eq.undo()
    if not eq.can_redo() or eq.redo_name() != "Typing":
        failures.append(f"redo lost the Typing name: {eq.redo_name()!r}")
    eq.redo()
    eq.insert_template("frac")
    if eq.undo_name() != "Template":
        failures.append(f"template history name was {eq.undo_name()!r}")
    eq.select_all()
    eq.restyle_selection("vector")
    if eq.undo_name() != "Style Change":
        failures.append(f"style history name was {eq.undo_name()!r}")
    eq.replace_latex("a+b", True)
    if eq.undo_name() != "TeX Edit":
        failures.append(f"source history name was {eq.undo_name()!r}")

    # Eqnedit32's persistent style chords are represented by actual TeX
    # typefaces, so the mode is not merely a visual flag in the native app.
    styled_cases = {
        "text": ("abc", r"\text{abc}"),
        "function": ("sin", r"\sin"),
        "variable": ("x", "x"),
        "roman": ("abc", r"\mathrm{abc}"),
        "italic": ("abc", r"\mathit{abc}"),
        "sans": ("abc", r"\mathsf{abc}"),
        "mono": ("abc", r"\mathtt{abc}"),
        "script": ("ABC", r"\mathcal{ABC}"),
        "double": ("ABC", r"\mathbb{ABC}"),
        "fraktur": ("ABC", r"\mathfrak{ABC}"),
        "boldsymbol": ("α", r"\bm{\alpha }"),
        "vector": ("x", r"\mathbf{x}"),
    }
    for style, (typed, expected) in styled_cases.items():
        eq = Equation()
        if not eq.insert_styled_text(typed, style):
            failures.append(f"styled input {style}: rejected")
        elif eq.latex().strip() != expected:
            failures.append(
                f"styled input {style}: got {eq.latex()!r}, want {expected!r}")
    eq = Equation()
    if eq.insert_styled_text("x", "unknown"):
        failures.append("unknown styled input mode was accepted")
    eq = Equation()
    eq.insert_styled_text("foo", "function")
    custom_function = eq.latex()
    again = Equation()
    again.load_latex(custom_function)
    if custom_function != r"\operatorname{foo}" or again.latex() != custom_function:
        failures.append(
            "custom function style was not stable: "
            f"{custom_function!r} -> {again.latex()!r}")

    # Exact Eqnedit32 automatic-Function dictionary: PE string resources
    # 11000--11038, loaded by x86 routine 0x00430670 as 0x27 entries.
    # Check both bulk model insertion and the character-by-character sequence
    # delivered by the real Win32 input path.
    auto_functions = {
        "Im": r"\operatorname{Im}",
        "Pr": r"\Pr",
        "Re": r"\operatorname{Re}",
        "arg": r"\arg",
        "arcsin": r"\arcsin",
        "arccos": r"\arccos",
        "arctan": r"\arctan",
        "cosh": r"\cosh",
        "cos": r"\cos",
        "coth": r"\coth",
        "cot": r"\cot",
        "cov": r"\operatorname{cov}",
        "csc": r"\csc",
        "deg": r"\deg",
        "det": r"\det",
        "dim": r"\dim",
        "exp": r"\exp",
        "gcd": r"\gcd",
        "glb": r"\operatorname{glb}",
        "hom": r"\hom",
        "inf": r"\inf",
        "int": r"\operatorname{int}",
        "ker": r"\ker",
        "ln": r"\ln",
        "lg": r"\lg",
        "lim": r"\lim",
        "log": r"\log",
        "lub": r"\operatorname{lub}",
        "max": r"\max",
        "min": r"\min",
        "mod": r"\operatorname{mod}",
        "sec": r"\sec",
        "sgn": r"\operatorname{sgn}",
        "sinh": r"\sinh",
        "sin": r"\sin",
        "sup": r"\sup",
        "tanh": r"\tanh",
        "tan": r"\tan",
        "var": r"\operatorname{var}",
    }
    eqnedit64_function_extensions = {
        "curl": r"\operatorname{curl}",
        "div": r"\operatorname{div}",
        "grad": r"\operatorname{grad}",
        "rot": r"\operatorname{rot}",
        "tr": r"\operatorname{tr}",
        "diag": r"\operatorname{diag}",
        "Res": r"\operatorname{Res}",
        "const": r"\operatorname{const}",
    }
    for typed, expected in {
            **auto_functions, **eqnedit64_function_extensions}.items():
        eq = Equation()
        eq.insert_text(typed)
        if eq.latex().strip() != expected:
            failures.append(
                f"automatic function {typed}: got {eq.latex()!r}, "
                f"want {expected!r}")
        sequential = Equation()
        for char in typed:
            sequential.insert_text(char)
        if sequential.latex().strip() != expected:
            failures.append(
                f"sequential automatic function {typed}: "
                f"got {sequential.latex()!r}, want {expected!r}")
    eq = Equation()
    eq.insert_text("single")
    if eq.latex() != "single":
        failures.append(
            f"ordinary word was partially recognised as a function: {eq.latex()!r}")
    eq = Equation()
    eq.insert_text("sinx")
    if eq.latex() != "sinx":
        failures.append(f"function prefix was recognised inside a word: {eq.latex()!r}")
    eq.backspace()
    if eq.latex().strip() != r"\sin":
        failures.append(
            f"deleting to a complete function did not restyle it: {eq.latex()!r}")

    # Real keyboard input arrives one character at a time.  Once `sin` has
    # become upright, the immediately following x is its italic argument and
    # must not turn the recognised function back into variables.
    eq = Equation()
    for char in "sinx":
        eq.insert_text(char)
    if eq.latex().strip() != r"\sin x":
        failures.append(
            f"typing sinx did not retain upright sin plus variable x: {eq.latex()!r}")

    eq = Equation()
    for char in "sinhx":
        eq.insert_text(char)
    if eq.latex().strip() != r"\sinh x":
        failures.append(
            f"typing sinhx did not prefer the longer function: {eq.latex()!r}")

    # Slots such as a fraction numerator are serialized through NodeList rather
    # than a top-level LineNode.  Style grouping must be identical there.
    nested_styles = {
        "function": ("sin", r"\frac{\sin}{}"),
        "text": ("速度", r"\frac{\text{速度}}{}"),
        "roman": ("abc", r"\frac{\mathrm{abc}}{}"),
        "italic": ("abc", r"\frac{\mathit{abc}}{}"),
        "vector": ("xy", r"\frac{\mathbf{x}\mathbf{y}}{}"),
    }
    for style, (typed, expected) in nested_styles.items():
        eq = Equation()
        eq.insert_template("frac")
        eq.insert_styled_text(typed, style)
        saved = eq.latex()
        again = Equation()
        again.load_latex(saved)
        if saved.replace(" ", "") != expected or again.latex() != saved:
            failures.append(
                f"nested styled input {style} was lost or unstable: "
                f"{saved!r} -> {again.latex()!r}")

    # A script inserted before any base uses an explicit empty TeX atom.  The
    # old emitter wrote bare `^{}`, which changed again on the second reload.
    eq = Equation()
    eq.insert_template("paren")
    eq.insert_template("sup")
    empty_base_script = eq.latex()
    once = Equation()
    once.load_latex(empty_base_script)
    twice = Equation()
    twice.load_latex(once.latex())
    if once.latex().strip() != twice.latex().strip() or "{}^{" not in empty_base_script:
        failures.append(
            "empty-base script was not stable through one reload: "
            f"{empty_base_script!r} -> {once.latex()!r} -> {twice.latex()!r}")

    # Deleting an empty outer script shell must not delete its visible base.
    # This is the exact regression reported for ``E = m{c^{2}}^{}``.
    eq = Equation()
    eq.load_latex(r"E = m{c^{2}}^{}")
    before_empty_sup = eq.latex()
    if not eq.backspace() or eq.latex().replace(" ", "") != r"E=mc^{2}":
        failures.append(
            "backspace deleted the base with an empty outer superscript: "
            f"{before_empty_sup!r} -> {eq.latex()!r}")
    if not eq.undo() or eq.latex() != before_empty_sup:
        failures.append("undo did not restore the removed empty superscript shell")

    eq.move_home()
    for _ in range(3):
        eq.move_right()
    if not eq.erase() or eq.latex().replace(" ", "") != r"E=mc^{2}":
        failures.append(
            "forward delete deleted the base with an empty outer superscript: "
            f"{eq.latex()!r}")

    mathml = tex_to_mathml(r"\frac{x_{1}+\alpha}{\sqrt{y}}")
    for required in ('display="inline"', 'mathsize="24pt"',
                     "<mfrac>", "<msub>", "<msqrt>"):
        if required not in mathml:
            failures.append(f"Office MathML is missing {required!r}: {mathml!r}")
    operator_mathml = tex_to_mathml(
        r"\sum_{n=1}^{m} a^3 \int_{a}^{b} f(x)\, dx^3")
    if not re.search(r"<munderover><mo[^>]*>&#x2211;</mo>", operator_mathml):
        failures.append(
            "Office MathML no longer puts display sum limits above/below: "
            f"{operator_mathml!r}")
    if not re.search(r"<msubsup><mo[^>]*>&#x222B;</mo>", operator_mathml):
        failures.append(
            "Office MathML does not match MathJax side limits for integrals: "
            f"{operator_mathml!r}")
    if "</munderover><msup>" not in operator_mathml or \
            "</msup><msubsup>" not in operator_mathml:
        failures.append(
            "Office MathML re-nested a large-operator body instead of matching "
            f"the browser's flat MathML import contract: {operator_mathml!r}")
    decoration_mathml = tex_to_mathml(r"\overline{u}+\underline{v}")
    if ("<mo stretchy=\"true\">&#x00AF;</mo>" not in decoration_mathml or
            "<mo stretchy=\"true\">_</mo>" not in decoration_mathml):
        failures.append(
            "Office MathML lost the shared overline/underline marks: "
            f"{decoration_mathml!r}")
    alphabet_mathml = tex_to_mathml(
        r"\mathrm{r}\mathit{i}\mathbf{v}\mathsf{s}\mathtt{t}"
        r"\mathcal{C}\mathbb{R}\mathfrak{F}\bm{\alpha}")
    for variant in ('mathvariant="normal"', 'mathvariant="italic"',
                    'mathvariant="bold"', 'mathvariant="sans-serif"',
                    'mathvariant="monospace"', 'mathvariant="script"',
                    'mathvariant="double-struck"',
                    'mathvariant="fraktur"'):
        if variant not in alphabet_mathml:
            failures.append(
                f"Office MathML lost {variant!r}: {alphabet_mathml!r}")

    unanchored_mathml = tex_to_mathml(
        r"\begin{aligned}a\\cccccc\end{aligned}")
    mathml_ns = {"m": "http://www.w3.org/1998/Math/MathML"}
    unanchored_root = ET.fromstring(unanchored_mathml)
    unanchored_rows = unanchored_root.findall("./m:mtable/m:mtr", mathml_ns)
    unanchored_cells = [row.findall("./m:mtd", mathml_ns)
                        for row in unanchored_rows]
    unanchored_office = tex_to_office_mathml_fragment(
        r"\begin{aligned}a\\cccccc\end{aligned}")
    if ('columnalign="left"' not in unanchored_mathml or
            len(unanchored_rows) != 2 or
            unanchored_office.count("<math ") != 2 or
            unanchored_office.count("<br>") != 1 or
            "mtable" in unanchored_office or
            "malign" in unanchored_office or "&amp;" in unanchored_office):
        failures.append(
            "one-column aligned Office fragment is not two clean math rows: "
            f"{unanchored_office!r}")
    nested_unanchored_office = tex_to_office_mathml_fragment(
        r"\begin{aligned}\begin{aligned}a^2\\\text{first}"
        r"\end{aligned}\\\text{second}\end{aligned}")
    if (nested_unanchored_office.count("<math ") != 3 or
            nested_unanchored_office.count("<br>") != 2 or
            "mtable" in nested_unanchored_office or
            "malign" in nested_unanchored_office or
            "&amp;" in nested_unanchored_office):
        failures.append(
            "nested one-column aligned Office fragment was not flattened: "
            f"{nested_unanchored_office!r}")
    anchored_mathml = tex_to_mathml(
        r"\begin{aligned}F&=ma\\E&=mc^2\end{aligned}")
    anchored_root = ET.fromstring(anchored_mathml)
    anchored_rows = anchored_root.findall("./m:mtable/m:mtr", mathml_ns)
    anchored_cells = [row.findall("./m:mtd", mathml_ns)
                      for row in anchored_rows]
    if ('columnalign="right left"' not in anchored_mathml or
            len(anchored_rows) != 2 or
            any(len(cells) != 2 for cells in anchored_cells) or
            "malign" in anchored_mathml):
        failures.append(
            "explicit alignment tabs lost invisible right/left cells: "
            f"{anchored_mathml!r}")
    anchored_office = tex_to_office_mathml_fragment(
        r"\begin{aligned}F&=ma\\E&=mc^2\end{aligned}")
    if anchored_office.count("<math ") != 1 or "<br>" in anchored_office:
        failures.append(
            "explicit aligned Office fragment was split into independent rows: "
            f"{anchored_office!r}")

    # Enter promotes a one-line equation to an aligned multi-line structure;
    # subsequent Enter and Up/Down operate on rows, not on raw TeX newlines.
    eq = Equation()
    eq.insert_latex("a=b")
    eq.new_line()
    eq.insert_latex("c=d")
    multi = eq.latex()
    if "\\begin{aligned}" not in multi or "a = b" not in multi or "c = d" not in multi:
        failures.append(f"multi-line equation was not structural: {multi!r}")
    if not eq.move_up() or not eq.move_down():
        failures.append("Up/Down did not traverse aligned rows")

    # A line break belongs at the caret, just as it does in the raw TeX pane.
    # The original implementation always moved the complete equation to row
    # one and appended a blank row, so Enter in the middle appeared to ignore
    # the insertion point.  A selection is replaced by the same row break.
    for selected, expected in (
            (False, r"\begin{aligned}abc\\def\end{aligned}"),
            (True, r"\begin{aligned}ab\\ef\end{aligned}")):
        eq = Equation()
        eq.insert_text("abcdef")
        eq.move_home()
        for _ in range(2 if selected else 3):
            eq.move_right()
        if selected:
            eq.begin_selection()
            eq.select_step_right()
            eq.select_step_right()
        if not eq.new_line():
            failures.append("Enter rejected a direct row split")
            continue
        if re.sub(r"\s+", "", eq.latex()) != expected:
            failures.append(
                f"Enter split at the wrong place: {eq.latex()!r}, "
                f"want {expected!r}")
        if eq.undo_name() != "Line Break":
            failures.append(
                f"row split has the wrong Undo name: {eq.undo_name()!r}")

    # Backspace at the start of the later row and Delete at the end of the
    # earlier row remove exactly the visible row separator.  Once only one
    # unanchored row remains, its redundant aligned wrapper disappears too.
    for key in ("backspace", "erase"):
        eq = Equation()
        eq.insert_text("abcdef")
        eq.move_home()
        for _ in range(3):
            eq.move_right()
        eq.new_line()
        if key == "erase":
            if not eq.move_up():
                failures.append("could not reach the row before Delete join")
                continue
            eq.move_end()
        if not getattr(eq, key)() or eq.latex() != "abcdef":
            failures.append(
                f"{key} did not join rows at the caret: {eq.latex()!r}")
        elif eq.caret() != ":3":
            failures.append(
                f"{key} row join lost its horizontal caret: {eq.caret()!r}")
        elif eq.undo_name() != "Join Lines":
            failures.append(
                f"{key} row join has wrong Undo name: {eq.undo_name()!r}")

    # An alignment tab is also a structural boundary.  Enter before it moves
    # both the current-cell suffix and the cells on its right to the new row,
    # matching insertion of `\\` at the same point in the TeX pane.
    eq = Equation()
    eq.load_latex(r"\begin{aligned}ab&=cd\\ef&=gh\end{aligned}")
    eq.move_home()
    eq.move_right()                    # enter the first cell
    eq.move_right()                    # between a and b
    eq.new_line()
    aligned_split = re.sub(r"\s+", "", eq.latex())
    if aligned_split != (
            r"\begin{aligned}a&\\b&=cd\\ef&=gh\end{aligned}"):
        failures.append(
            "Enter did not carry right-hand aligned cells to the new row: "
            f"{eq.latex()!r}")

    # Joining an explicitly aligned row preserves the column schema.  The
    # corresponding cells join pairwise; flattening all four cells into the
    # root would discard the user's `&` contract.
    for key in ("backspace", "erase"):
        eq = Equation()
        eq.load_latex(r"\begin{aligned}a&b\\c&d\end{aligned}")
        eq.move_home()
        eq.move_right()                  # first row, first cell
        if key == "backspace":
            eq.move_down()               # second row, first cell start
        else:
            eq.next_slot()               # first row, final cell
            eq.move_end()
        if not getattr(eq, key)():
            failures.append(f"{key} rejected a two-column row join")
            continue
        joined = re.sub(r"\s+", "", eq.latex())
        if joined != r"\begin{aligned}ac&bd\end{aligned}":
            failures.append(
                f"{key} flattened or lost aligned columns: {eq.latex()!r}")

    # Joining can create a function word just as typing its final letter can.
    # Re-run automatic classification at the old boundary in both directions.
    for key in ("backspace", "erase"):
        eq = Equation()
        eq.load_latex(r"\begin{aligned}s\\in\end{aligned}")
        eq.move_home()
        eq.move_right()
        if key == "backspace":
            eq.move_down()
        else:
            eq.move_end()
        getattr(eq, key)()
        if eq.latex() != r"\sin":
            failures.append(
                f"{key} row join did not refresh function style: "
                f"{eq.latex()!r}")

    # Automatic function styling is a property of a complete word.  Splitting
    # sin after s must reclassify both fragments instead of leaving stale
    # upright letters on either row.
    eq = Equation()
    eq.insert_text("sin")
    eq.move_home()
    eq.move_right()
    eq.new_line()
    function_split = eq.latex()
    if ("operatorname" in function_split or "\\sin" in function_split or
            re.sub(r"\s+", "", function_split) !=
            r"\begin{aligned}s\\in\end{aligned}"):
        failures.append(
            f"function style survived across a row split: {function_split!r}")

    # Enter deliberately leaves a blank row ready for typing.  Saving before
    # filling it must not write a spurious final `\\`, and must reach the same
    # canonical TeX after one load/save cycle.
    eq = Equation()
    eq.insert_latex("a=b")
    eq.new_line()
    trailing_blank = eq.latex()
    once = Equation()
    once.load_latex(trailing_blank)
    if once.latex().strip() != trailing_blank.strip():
        failures.append(
            "trailing aligned edit row was serialized as mathematical content: "
            f"{trailing_blank!r} -> {once.latex()!r}")

    # An ampersand is an alignment command in the editor, and remains a real
    # TeX alignment tab when saved.
    eq = Equation()
    eq.insert_text("F")
    eq.alignment_tab()
    eq.insert_latex("=ma")
    if "F &  = ma" not in eq.latex() and "F & = ma" not in eq.latex():
        failures.append(f"alignment tab was lost: {eq.latex()!r}")

    # Every template must be insertable and must survive a LaTeX round trip
    # once its holes are filled.  An *empty* hole is a different matter: LaTeX
    # has no way to write "there is a row here and it is blank", so a template
    # is only required to be stable in the state a person leaves it in.
    for kind in Equation.templates():
        eq = Equation()
        if not eq.insert_template(kind):
            failures.append(f"template {kind}: insert_template returned False")
            continue
        eq.insert_text("x")
        for _ in range(8):
            if not eq.next_slot():
                break
            eq.insert_text("x")
        latex = eq.latex()
        again = Equation()
        again.load_latex(latex)
        if again.latex().strip() != latex.strip():
            failures.append(
                f"template {kind}: not stable through LaTeX\n"
                f"    {latex!r}\n    {again.latex()!r}")

    # A matrix is an editable table, not a choice between two fixed stamps.
    # Any practical rectangular size can be created by the model, while the
    # four structural operations preserve every cell outside the edited row
    # or column and participate in normal Undo/Redo history.
    eq = Equation()
    if not eq.insert_template("matrix7x9") or eq.matrix_dimensions() != (7, 9):
        failures.append(
            f"arbitrary 7x9 matrix was not created: {eq.matrix_dimensions()!r}")
    invalid = Equation()
    if (invalid.insert_template("matrix0x2") or
            invalid.insert_template("matrix2x0") or
            invalid.insert_template("matrix100x2") or
            invalid.insert_template("matrix2x100") or
            invalid.insert_template("matrix2x3junk")):
        failures.append("an invalid/out-of-range matrix size was accepted")
    limit = Equation()
    if (not limit.insert_template("matrix99x99") or
            limit.matrix_dimensions() != (99, 99) or
            limit.matrix_add_row() or limit.matrix_add_column()):
        failures.append("the documented 99x99 matrix limit was not enforced")

    # Empty edge cells are structural too.  Canonical `{}` placeholders make
    # a 1x3 row and a 3x1 column survive ordinary TeX save/reopen instead of
    # shrinking once on every normalization pass.
    for kind in ("matrix1x3", "matrix3x1"):
        edge = Equation()
        edge.insert_template(kind)
        edge.insert_text("x")
        edge_tex = edge.latex()
        reopened = Equation()
        reopened.load_latex(edge_tex)
        renormalized = reopened.latex()
        reopened_again = Equation()
        reopened_again.load_latex(renormalized)
        if renormalized != edge_tex or reopened_again.latex() != edge_tex:
            failures.append(
                f"{kind} empty edge cells shrank across TeX: "
                f"{edge_tex!r} -> {renormalized!r} -> {reopened_again.latex()!r}")

    eq = Equation()
    eq.insert_template("matrix2x3")
    for position, value in enumerate("abcdef"):
        eq.insert_text(value)
        if position != 5:
            eq.next_slot()
    original_matrix = eq.latex()
    if eq.matrix_dimensions() != (2, 3):
        failures.append(f"2x3 matrix dimensions were lost: {eq.matrix_dimensions()!r}")
    if not eq.move_up() or not eq.move_down():
        failures.append("Up/Down did not navigate between ordinary matrix rows")
    if not eq.matrix_add_row() or eq.matrix_dimensions() != (3, 3):
        failures.append(f"matrix row was not added: {eq.matrix_dimensions()!r}")
    elif eq.undo_name() != "Matrix Row":
        failures.append(f"matrix row Undo name was {eq.undo_name()!r}")
    if not eq.matrix_add_column() or eq.matrix_dimensions() != (3, 4):
        failures.append(f"matrix column was not added: {eq.matrix_dimensions()!r}")
    elif eq.undo_name() != "Matrix Column":
        failures.append(f"matrix column Undo name was {eq.undo_name()!r}")
    if any(value not in eq.latex() for value in "abcdef"):
        failures.append(f"matrix resize lost an existing cell: {eq.latex()!r}")
    if not eq.undo() or eq.matrix_dimensions() != (3, 3):
        failures.append("Undo did not restore the pre-column matrix")
    if not eq.undo() or eq.matrix_dimensions() != (2, 3) or \
            eq.latex() != original_matrix:
        failures.append(
            f"Undo did not restore the original 2x3 matrix: {eq.latex()!r}")
    if not eq.redo() or not eq.redo() or eq.matrix_dimensions() != (3, 4):
        failures.append("Redo did not restore the rectangular matrix resize")

    eq = Equation()
    eq.insert_template("matrix2x2")
    if (not eq.matrix_remove_row() or eq.matrix_dimensions() != (1, 2) or
            eq.matrix_remove_row()):
        failures.append("matrix rows did not stop at the one-row minimum")
    if (not eq.matrix_remove_column() or eq.matrix_dimensions() != (1, 1) or
            eq.matrix_remove_column()):
        failures.append("matrix columns did not stop at the one-column minimum")

    # Shortcut declarations are a contract shared by Help, the native GUI,
    # and this headless model.  Duplicate chords are ambiguous even if both
    # targets happen to work.
    chords = {}
    for chord, cmd, _label in Equation.shortcuts():
        if chord in chords:
            failures.append(f"duplicate shortcut {chord}: {chords[chord]} / {cmd}")
        chords[chord] = cmd
        eq = Equation()
        eq.insert_text("x")
        if cmd.startswith("template.") and cmd[9:] not in Equation.templates():
            failures.append(f"{chord} -> {cmd}: no such template")
        if cmd.startswith("matrix."):
            eq.load_latex("")
            eq.insert_template("matrix2x2")
            if not eq.command(cmd):
                failures.append(f"{chord} -> {cmd}: command did nothing")

    # Every two-stroke sequence must resolve through the same table and its
    # resulting model command must execute.  This catches a Help-only binding
    # and a key-handler-only binding in either direction.
    sequences = Equation.sequence_shortcuts()
    sequence_keys = set()
    for prefix, key, shift, cmd, _label in sequences:
        identity = (prefix, key, shift)
        if identity in sequence_keys:
            failures.append(f"duplicate sequence {identity}")
        sequence_keys.add(identity)
        resolved = Equation.resolve_sequence(prefix, key, shift)
        if resolved != cmd:
            failures.append(f"sequence {identity}: resolves to {resolved!r}, want {cmd!r}")
            continue
        eq = Equation()
        eq.insert_text("x")
        if not eq.command(cmd):
            failures.append(f"sequence {identity} -> {cmd}: command did nothing")

    known_sequences = {
        ("T", "R", False): "template.sqrt",
        ("T", "N", False): "template.nthroot",
        ("T", "M", False): "template.matrix3x3",
        ("T", "/", False): "template.slashfrac",
        ("T", "U", False): "template.under",
        ("T", "C", False): "template.cases",
        ("T", "3", False): "template.matrix3x3",
        ("K", "I", False): "symbol.\\infty",
        ("K", "E", True): "symbol.\\notin",
        ("B", "X", False): "latex.\\mathbf{x}",
        ("B", "X", True): "latex.\\mathbf{X}",
        ("G", "A", False): "symbol.\\alpha",
        ("G", "G", True): "symbol.\\Gamma",
        ("G", "O", False): "symbol.\\omega",
    }
    for identity, expected in known_sequences.items():
        got = Equation.resolve_sequence(*identity)
        if got != expected:
            failures.append(f"known sequence {identity}: got {got!r}, want {expected!r}")

    if Equation.resolve_sequence("T", "Z", False):
        failures.append("undefined Ctrl+T,Z unexpectedly resolved")

    # Shift+arrow selects a script as a unit and, once the caret is inside a
    # script, selects its content.  Movement used to descend into the
    # structure mid-drag, so Shift+arrow across a slot boundary silently
    # selected nothing -- the superscript 2 of a^2 could not be selected.
    eq = Equation()
    eq.load_latex("a^{2}")
    eq.move_end()
    eq.move_left()                         # into the superscript, after the 2
    if not eq.select_step_left() or eq.selection_latex() != "2":
        failures.append(
            f"Shift+Left in the exponent did not select the 2: "
            f"{eq.selection_latex()!r}")

    eq = Equation()
    eq.load_latex("a^{2}")
    eq.move_home()
    if not eq.select_step_right() or eq.selection_latex().replace(" ", "") \
            != "a^{2}":
        failures.append(
            f"Shift+Right before a^2 did not select the whole script: "
            f"{eq.selection_latex()!r}")

    # At a slot boundary the step does nothing and must not discard an
    # existing selection.
    eq = Equation()
    eq.load_latex("xy")
    eq.move_home()
    eq.select_step_right()                  # selects x
    held = eq.selection_latex()
    eq.select_step_left()                   # back to the anchor: empty again
    eq.select_step_right()
    if eq.selection_latex() != held:
        failures.append("selection stepping is not reversible")

    # Applying a sub/superscript to a whole selection and then wanting it off
    # again.  The caret sits at the end of the base with the empty script
    # boxes ahead of it; Delete removes those boxes and keeps the base, the
    # mirror of Backspace unwrapping a template from the start of a slot.
    # Reported as the scripts being unremovable from `{...}_{}^{}` (caret at
    # 0.0:8 in the screenshot).
    def base_end_after_subsup(content):
        eq = Equation()
        for c in content:
            eq.insert_text(c)
        eq.select_all()
        eq.insert_template("subsup")
        eq.command("caret.left")            # sub slot -> end of base
        return eq

    eq = base_end_after_subsup("abc")
    if eq.caret() != "0.0:3":
        failures.append(f"expected the caret at the base end, got {eq.caret()}")
    if not eq.erase() or eq.latex().replace(" ", "") != "abc":
        failures.append(
            "Delete did not remove the empty scripts from the base end: "
            f"{eq.latex()!r}")

    # Backspace there still edits the base -- it must not become a second way
    # to unwrap.
    eq = base_end_after_subsup("abc")
    if not eq.backspace() or "_{}" not in eq.latex():
        failures.append(
            f"Backspace at the base end should edit the base: {eq.latex()!r}")

    # A script that still carries content on one side keeps it.
    eq = Equation()
    eq.load_latex(r"x_{i}^{}")
    eq.command("caret.right")               # into the base
    eq.command("caret.end")
    while not eq.caret().startswith("0.0"):
        eq.command("caret.left")
    eq.erase()
    if "i" not in eq.latex():
        failures.append(
            f"Delete destroyed a non-empty subscript: {eq.latex()!r}")

    # One Backspace removes one thing the reader can see.
    #
    # This is stated as a rule on purpose.  The earlier defect here was
    # reported as `E = m{c^{2}}^{}` losing its base, and both the fix and its
    # regression were written to that one state: a branch for "the previous
    # node is a script with an empty slot", and assertions with the caret
    # either inside a script or in an empty slot.  Every other state kept the
    # general path, which erased the whole node -- so `E = mc^{2}` with the
    # caret after it lost the c together with the 2, and nothing could tell.
    #
    # Counting glyphs in the rendered SVG measures "visible item" directly and
    # does not care which template is involved, so a new template cannot
    # quietly fall outside the rule.
    # An empty template may go as a unit -- a bare "( )" is one object to the
    # reader even though it draws two glyphs -- so the rule binds only while
    # the expression still holds content.
    style = SvgStyle()

    def visible_glyphs(latex):
        return tex_to_svg(latex, style).count("<text")

    def holds_content(latex):
        return bool(re.search(r"[A-Za-z0-9]", re.sub(r"\\[a-zA-Z]+", "", latex)))

    # Both keys, not one.  The rule was first written for Backspace alone and
    # Delete kept the old behaviour -- `c^{2}` lost its c and its 2 to a
    # single press -- which a test that exercised only Backspace could not
    # see.  That is the same mistake this rule exists to prevent, so the two
    # are checked together.
    for key in ("backspace", "delete"):
        for start in (r"E = mc^{2}a", r"\frac{a}{b}", r"\sqrt{x^{2}}",
                      r"\left( xy \right)", r"\sum_{k=1}^{n}x",
                      r"\sqrt{\frac{a}{b}}", r"x_{i}^{2}",
                      r"\begin{cases} x & x>0 \\ -x & x<0 \end{cases}"):
            eq = Equation()
            eq.load_latex(start)
            if key == "backspace":
                eq.move_end()
                press = eq.backspace
            else:
                eq.move_home()
                press = eq.erase
            before = eq.latex()
            seen = visible_glyphs(before)
            for step in range(60):
                if not press():
                    break
                after = eq.latex()
                now = visible_glyphs(after)
                if seen - now > 1 and holds_content(before):
                    failures.append(
                        f"{start!r}: {key} {step + 1} removed {seen - now} "
                        f"visible items at once, {before!r} -> {after!r}")
                    break
                before, seen = after, now
            else:
                failures.append(
                    f"{start!r}: 60 {key} presses did not empty it, "
                    f"still {eq.latex()!r}")

    # Editing the raw TeX pane used to go through load_latex(), which clears
    # the undo stack: one keystroke there discarded everything the canvas had
    # done.  replace_latex() keeps the history and takes one checkpoint per
    # editing burst, so Ctrl+Z returns to the equation as it stood before the
    # pane was touched -- not one character at a time, and not nowhere.
    pane = Equation()
    pane.load_latex("a+b")
    pane.insert_template("frac")          # a canvas edit worth getting back
    canvas_state = pane.latex()
    for burst, text in enumerate(("a+b+", "a+b+c")):
        pane.replace_latex(text, burst == 0)
    if pane.latex() != "a+b+c":
        failures.append(f"source pane edit did not apply: {pane.latex()!r}")
    if not pane.undo():
        failures.append("source pane edit left no undo step")
    elif pane.latex() != canvas_state:
        failures.append(
            f"undo after a source burst gave {pane.latex()!r}, "
            f"want {canvas_state!r}")
    if not pane.undo():
        failures.append("undo did not reach the canvas edit before the pane")

    # --- a deletion leaves the caret where it happened -------------------
    # Reported for x_1^2: Backspace removed the 2 correctly but left the
    # caret outside the script, so typing put the next character after the
    # whole thing (x_{1}^{}9) instead of back in the exponent.  Written for
    # both keys at once, deliberately: the two previous defects in this area
    # were each fixed for one key and left standing for the other.
    for tex, expect_bs, expect_del in (
            (r"x_{1}^{2}", r"x_{1}^{9}", r"9_{1}^{2}"),
            # Forward from the start, the first thing a reader sees in
            # x^{ab} is the base x, not the exponent.
            (r"x^{ab}", r"x^{a9}", r"9^{ab}"),
            (r"\sqrt{ab}", r"\sqrt{a9}", r"\sqrt{9b}"),
            (r"\frac{ab}{c}", r"\frac{ab}{9}", r"\frac{9b}{c}")):
        for key, expect in (("backspace", expect_bs), ("erase", expect_del)):
            eq = Equation()
            eq.load_latex(tex)
            if key == "backspace":
                eq.move_end()
            else:
                eq.move_home()
            getattr(eq, key)()
            eq.insert_text("9")
            if eq.latex() != expect:
                failures.append(
                    f"{tex}: {key} then typing gave {eq.latex()!r}, want "
                    f"{expect!r} -- the caret did not follow the deletion")

    # --- and repeating the key still clears the equation ------------------
    # Moving the caret into the slot broke this the first time: once the slot
    # emptied there was nothing ahead of (or behind) the caret and the key
    # stopped working, leaving half the equation behind.
    for tex in (r"x_{1}^{2}", r"x^{ab}", r"\sqrt{ab}", r"\frac{ab}{c}",
                r"\int _{a}^{b}f(x)dx", r"\sum _{i=1}^{n}x_{i}"):
        for key in ("backspace", "erase"):
            eq = Equation()
            eq.load_latex(tex)
            if key == "backspace":
                eq.move_end()
            else:
                eq.move_home()
            for _ in range(60):
                if not getattr(eq, key)():
                    break
            if eq.latex().strip():
                failures.append(
                    f"{tex}: repeating {key} stopped with {eq.latex()!r} left")

    # --- editor/parser depth is one shared, lossless boundary -------------
    # Undo snapshots are LaTeX.  Before this gate, Ctrl+R could build more
    # radicals than the parser accepted; the next Undo silently restored only
    # the first 200.  The rejected operation must now be a complete no-op,
    # including caret and history, and the largest accepted tree must remain a
    # parse/emit fixed point that layout and SVG can traverse safely.
    deep = Equation()
    for level in range(MAX_NESTING_DEPTH):
        before_level = deep.latex()
        if not deep.insert_template("sqrt"):
            failures.append(f"nesting limit rejected safe level {level + 1}")
            break
        accepted_level = deep.latex()
        probe = Equation()
        if not probe.load_latex(accepted_level) or probe.latex() != accepted_level:
            failures.append(
                f"nesting level {level + 1} is not a parse/emit fixed point")
            break
        probe.metrics()
        if not deep.undo() or deep.latex() != before_level:
            failures.append(
                f"nesting level {level + 1} did not Undo losslessly")
            break
        if not deep.redo() or deep.latex() != accepted_level:
            failures.append(
                f"nesting level {level + 1} did not Redo losslessly")
            break
    accepted = deep.latex()
    accepted_caret = deep.caret()
    accepted_undo = deep.undo_name()
    if accepted.count(r"\sqrt{") != MAX_NESTING_DEPTH:
        failures.append("largest accepted editor tree has the wrong depth")
    if deep.insert_template("sqrt"):
        failures.append("template insertion exceeded the shared nesting limit")
    if (deep.latex(), deep.caret(), deep.undo_name()) != (
            accepted, accepted_caret, accepted_undo):
        failures.append("rejected template mutated equation/caret/Undo state")
    if deep.last_error() != "maximum-nesting-depth":
        failures.append("depth rejection did not expose its diagnostic")

    reopened = Equation()
    if not reopened.load_latex(accepted) or reopened.latex() != accepted:
        failures.append("maximum safe nesting is not a parse/emit fixed point")
    else:
        reopened.metrics()
        reopened.svg()
    if not deep.undo():
        failures.append("rejected template consumed the preceding Undo step")
    elif deep.latex().count(r"\sqrt{") != MAX_NESTING_DEPTH - 1:
        failures.append("Undo after depth rejection did not restore level 199")
    elif not deep.redo() or deep.latex() != accepted:
        failures.append("Redo after depth rejection did not restore level 200")

    too_deep = r"\sqrt{" * (MAX_NESTING_DEPTH + 1) + "x" + \
        "}" * (MAX_NESTING_DEPTH + 1)
    guarded = Equation()
    guarded.insert_text("a")
    state = (guarded.latex(), guarded.caret(), guarded.undo_name())
    if guarded.insert_latex(too_deep):
        failures.append("structural paste accepted parser-truncating input")
    if (guarded.latex(), guarded.caret(), guarded.undo_name()) != state:
        failures.append("rejected structural paste mutated editor state")
    if guarded.replace_latex(too_deep, True):
        failures.append("source edit accepted parser-truncating input")
    if (guarded.latex(), guarded.caret(), guarded.undo_name()) != state:
        failures.append("rejected source edit mutated editor state")
    if guarded.load_latex(too_deep):
        failures.append("document load accepted parser-truncating input")
    if (guarded.latex(), guarded.caret(), guarded.undo_name()) != state:
        failures.append("rejected document load mutated editor state")

    total = (len(CASES) + len(Equation.templates()) + len(sequences) + 39 + 32 +
             selection_contract_checks)
    if failures:
        print(f"FAIL  {len(failures)} of {total}")
        for f in failures:
            print("  " + f)
        return 1
    print(f"ok    {total} checks: typing, caret, templates, undo, shortcuts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
