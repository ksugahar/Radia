"""The editing model, exercised the way a person types.

Every case is a keystroke sequence rather than an API call sequence, because
what is being checked is that Equation Editor's *feel* survived: Ctrl+L after
`x` gives `x` a subscript instead of an empty box, Tab walks the holes of a
template, backspace at the start of an empty slot unwraps the template instead
of swallowing it, and undo puts everything back.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
Equation = equation.Equation


def type_in(eq, keys):
    """A key is a literal string to type, a \\command, or a command name."""
    for k in keys:
        if k.startswith("\\"):
            assert eq.insert_symbol(k), f"unknown symbol {k}"
        elif "." in k:
            assert eq.command(k), f"command {k} did nothing"
        else:
            eq.insert_text(k)
    return eq


TYPING = [
    ("plain text", ["x", "+", "y"], "x+y"),
    ("Ctrl+L makes the letter just typed the base",
     ["x", "template.sub", "i"], "x_{i}"),
    ("Ctrl+H after a letter", ["a", "template.sup", "2"], "a^{2}"),
    ("Ctrl+J gives both scripts, Tab moves between them",
     ["x", "template.subsup", "i", "caret.next_slot", "2"], "x_{i}^{2}"),
    ("fraction: numerator, Tab, denominator",
     ["template.frac", "a", "caret.next_slot", "b"], r"\dfrac{a}{b}"),
    ("leave the fraction and keep typing",
     ["template.frac", "1", "caret.next_slot", "2", "caret.next_slot", "+", "3"],
     r"\dfrac{1}{2}+3"),
    ("square root", ["template.sqrt", "2"], r"\sqrt{2}"),
    ("a fraction inside a square root",
     ["template.sqrt", "template.frac", "a", "caret.next_slot", "b"],
     r"\sqrt{\dfrac{a}{b}}"),
    ("parentheses", ["template.paren", "x", "+", "1"], r"\left( x+1 \right)"),
    # LaTeX has nowhere to record "limits are stacked"; display style stacks
    # them anyway.
    ("summation with both limits",
     ["template.sum", "i", "caret.next_slot", "n", "caret.next_slot", "a"],
     r"\sum _{i}^{n} a"),
    ("Greek by name", [r"\alpha", "+", r"\beta"], r"\alpha +\beta "),
    ("backspace deletes the character before the caret",
     ["a", "b", "edit.backspace"], "a"),
    ("backspace in an empty slot unwraps the template, keeping the content",
     ["template.frac", "a", "caret.next_slot", "caret.prev_slot",
      "caret.home", "edit.backspace"], "a"),
    ("undo restores what the last edit changed",
     ["a", "template.frac", "b", "edit.undo", "edit.undo"], "a"),
    ("redo after undo", ["a", "b", "edit.undo", "edit.redo"], "ab"),
]


@pytest.mark.parametrize("name,keys,expect",
                         TYPING, ids=[c[0] for c in TYPING])
def test_typing(name, keys, expect):
    assert type_in(Equation(), keys).latex().strip() == expect.strip()


def test_walking_the_caret_never_changes_the_equation():
    eq = type_in(Equation(), ["template.frac", "a", "caret.next_slot", "b"])
    for _ in range(40):
        eq.move_left()
    for _ in range(40):
        eq.move_right()
    assert eq.latex().strip() == r"\dfrac{a}{b}"


@pytest.mark.parametrize("kind", Equation.templates())
def test_template_survives_a_latex_round_trip(kind):
    """Filled templates round-trip exactly.

    An *empty* hole is a different matter: LaTeX has no way to write "there is
    a row here and it is blank", so a template is only required to be stable in
    the state a person leaves it in.
    """
    eq = Equation()
    assert eq.insert_template(kind)
    eq.insert_text("x")
    for _ in range(8):
        if not eq.next_slot():
            break
        eq.insert_text("x")
    latex = eq.latex()

    again = Equation()
    again.load_latex(latex)
    assert again.latex().strip() == latex.strip()


@pytest.mark.parametrize("chord,command,label", Equation.shortcuts())
def test_shortcut_names_a_real_command(chord, command, label):
    if command.startswith("template."):
        assert command[len("template."):] in Equation.templates()
    else:
        eq = Equation()
        eq.insert_text("xy")
        eq.command(command)          # must not raise; may return False


def test_load_clears_the_history():
    eq = Equation()
    eq.insert_text("a")
    eq.load_latex("b")
    assert not eq.undo()
    assert eq.latex().strip() == "b"


def test_the_editing_model_and_the_office_output_agree():
    eq = type_in(Equation(), ["template.frac", "a", "caret.next_slot", "b"])
    assert eq.omml() == equation.tex_to_omml(eq.latex())
