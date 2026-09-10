"""Independent face -> TeX specification; safe static checks on LAB.

The isolated model runner also calls main() to execute every actual command.
No import of the native module occurs during pytest collection.
"""
import json
import os
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
SPEC = json.loads((ROOT / "docs/palette_intent.json").read_text(encoding="utf-8"))
STR = r'"(?:[^"\\]|\\.)*"'


def symbol_contract():
    result = {}
    for entry in SPEC["symbols"].split():
        face, names = entry.split(":")
        for name in names.split(","):
            assert "\\" + name not in result
            result["\\" + name] = face
    return result


def source_cells(source=None):
    if source is None:
        source = (ROOT / "src/palettes.cpp").read_text(encoding="utf-8")
    pattern = re.compile(
        r'(?P<kind>sym|tpl|raw)\(\s*(?P<cmd>' + STR + r')\s*,\s*u8(?P<face>' + STR + r')'
        r'|I\{\s*(?P<direct>' + STR + r')\s*,\s*u8(?P<direct_face>' + STR + r')')
    cells = []
    for m in pattern.finditer(source):
        if m["kind"]:
            prefix = {"sym": "symbol", "tpl": "template", "raw": "latex"}[m["kind"]]
            cells.append((prefix + "." + json.loads(m["cmd"]), json.loads(m["face"])))
        else:
            cells.append((json.loads(m["direct"]), json.loads(m["direct_face"])))
    return cells


def matrix_tex(rows, environment="matrix"):
    return "\\begin{" + environment + "}" + r"\\".join(
        "&".join(row) for row in rows) + "\\end{" + environment + "}"


def native_contract():
    # command -> (face, keyboard input, independently authored expected TeX)
    contract = {"symbol." + name: (face, "", name)
                for name, face in symbol_contract().items()}
    contract.update({"template." + name: tuple(row)
                     for name, row in SPEC["native_templates"].items()})
    for size in SPEC["native_matrix_sizes"]:
        rows, cols = map(int, size.split("x"))
        # Single characters, all different, so swapped cells cannot pass.
        values = "abcdefghijklmnopqrstuvwxyz0123456789"[:rows * cols]
        table = [list(values[i * cols:(i + 1) * cols]) for i in range(rows)]
        contract["template.matrix" + size] = (size.replace("x", "×"), values, matrix_tex(table))
    for name, (face, tex) in SPEC["native_styles"].items():
        contract["style." + name] = (face, "", tex)
    for face, tex in SPEC["native_raw"].items():
        contract["latex." + tex] = (face, "", tex)
    for name, (face, rows) in SPEC["matrix_actions"].items():
        contract["matrix." + name] = (face, "", matrix_tex(rows))
    return contract


def check_catalogue(cells):
    contract = native_contract()
    commands = [command for command, _ in cells]
    assert len(commands) == len(set(commands)), "duplicate palette command"
    assert set(commands) == set(contract), (
        "missing intent or removed key", sorted(set(commands) - set(contract)),
        sorted(set(contract) - set(commands)))
    for command, face in cells:
        assert face == contract[command][0], (command, face, contract[command][0])


def test_complete_native_intent_catalogue():
    check_catalogue(source_cells())


def test_swapped_commands_and_unreviewed_keys_fail():
    import pytest
    cells = source_cells()
    for left, right in [("symbol.\\leftarrow", "symbol.\\rightarrow"),
                        ("template.overline", "template.underline"),
                        ("template.sup", "template.sub"),
                        ("style.sans", "style.mono")]:
        swapped = [(right if cmd == left else left if cmd == right else cmd, face)
                   for cmd, face in cells]
        with pytest.raises(AssertionError):
            check_catalogue(swapped)
    with pytest.raises(AssertionError):
        check_catalogue(cells + [("template.unreviewed", "?")])


def fill_slots(equation, values):
    for index, value in enumerate(values):
        if index:
            assert equation.next_slot(), ("missing slot", index, values)
        equation.insert_text(value)


def execute(E, command, values):
    eq = E.Equation()
    kind, payload = command.split(".", 1)
    if kind == "template":
        assert eq.insert_template(payload), command
        fill_slots(eq, values)
    elif kind == "symbol":
        assert eq.insert_symbol(payload), command
    elif kind == "latex":
        assert eq.insert_latex(payload), command
    elif kind == "style":
        eq.insert_text("a")
        eq.select_all()
        assert eq.restyle_selection(payload), command
    elif kind == "matrix":
        assert eq.insert_template("matrix2x2")
        fill_slots(eq, "abcd")
        for _ in range(3):
            assert eq.prev_slot(), command
        assert eq.command(command), command
    else:
        raise AssertionError(command)
    return eq.latex()


def main():
    if os.environ.get("EQNEDIT64_ISOLATED_TEST_SESSION") != "1":
        print("FAIL: native palette intent requires an isolated font session")
        return 90
    import eqnedit_core as E
    check_catalogue([(command, face) for _, _, _, items in E.palettes()
                     for command, face, _ in items])
    failures = []
    contract = native_contract()
    for command, (_, values, expected) in contract.items():
        try:
            actual = execute(E, command, values)
            # Only spelling canonicalization is delegated to the parser.
            # The expected operator, direction and operands come from SPEC,
            # never from the command being tested or from its output.
            canonical = E.tex_normalize(expected)
            assert actual == canonical, (actual, canonical, expected)
            assert E.tex_normalize(actual) == actual, "output is not stable"
        except AssertionError as exc:
            failures.append(f"{command}: {exc}")
    # Prove that a correctly round-tripping but wrong insertion is rejected.
    for left, right in [("template.overline", "template.underline"),
                        ("template.sup", "template.sub"),
                        ("symbol.\\leftarrow", "symbol.\\rightarrow")]:
        wrong = execute(E, right, contract[right][1])
        assert wrong != E.tex_normalize(contract[left][2]), (left, right)
    for failure in failures:
        print("FAIL: " + failure)
    print(f"Native palette intent: {len(contract)} keys, {len(failures)} failures")
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
