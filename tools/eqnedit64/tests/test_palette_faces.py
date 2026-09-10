# -*- coding: utf-8 -*-
"""A palette key must show the character it inserts, in the shipped font.

Two failures hid behind each other in 3.0.16, and both were invisible to every
existing test because the suites press the COMMAND behind a key and never look
at the KEY.

1.  The differential-geometry Hodge-star key showed U+2605 BLACK STAR while
    `\\star` renders as U+22C6 STAR OPERATOR, not U+2605 BLACK STAR.  The
    key promised one glyph and inserted another.

2.  U+2605 is not in Latin Modern Math.  `pick_button_font` concatenates
    EVERY palette face into one sample and demands a font that owns the whole
    sample, so a single unavailable character does not blank one key -- it
    rejects Latin Modern Math outright and drops the ENTIRE palette to the
    next candidate.  Cambria Math fails the same sample, so all 245 keys were
    drawn in Segoe UI Symbol.  That is why the prime family read as typewriter
    quotes: a reader could no longer tell x'' from x-double-prime.

The all-or-nothing font gate is the right design -- a half-rendered palette
would be worse -- but it makes one wrong character cost the whole window, so
the sample has to be checked rather than trusted.

Both checks are static: they read the two source files and the font's cmap
table.  Nothing here imports eqnedit_core, so no font is registered and this
file is safe to run in an ordinary desktop session and on a runner with no
compiler.

Run:  python -m pytest tests/test_palette_faces.py -q
"""
from pathlib import Path
import re

from fontTools.ttLib import TTFont

ROOT = Path(__file__).resolve().parents[1]
PALETTES = ROOT / "src/palettes.cpp"
SYMBOLS = ROOT / "src/math_symbols.cpp"
MATH_FONT = ROOT / "assets/latinmodern-math.otf"

# A C++ string literal, captured raw.
STR = r'"((?:[^"\\]|\\.)*)"'
CELL = re.compile(r"(?:tpl|sym|raw)\(\s*(?:u8)?" + STR + r"\s*,\s*u8" + STR)
STYLE_CELL = re.compile(r"I\{\s*" + STR + r"\s*,\s*u8" + STR)
PALETTE_TAB = re.compile(r"p\.push_back\(\{u8" + STR + r"\s*,\s*u8" + STR)
SYMBOL_ENTRY = re.compile(r'\{' + STR + r"\s*,\s*(0x[0-9A-Fa-f]+)\}")

# `\not` inserts U+0338 COMBINING LONG SOLIDUS OVERLAY, which has no standalone
# shape: drawn alone it lands on the key's own edge.  The word is deliberate.
FACE_IS_A_WORD = {"\\not"}


def unescape(literal):
    """The value a C++ compiler would give the captured literal."""
    return literal.replace('\\\\', '\\').replace('\\"', '"')


def cells():
    """Every (command, face) the palettes declare, tab faces included."""
    source = PALETTES.read_text(encoding="utf-8")
    for pattern in (CELL, STYLE_CELL):
        for match in pattern.finditer(source):
            yield unescape(match.group(1)), unescape(match.group(2))
    for match in PALETTE_TAB.finditer(source):
        # A tab has a title, not a command; it still lands in the font sample.
        yield "tab:" + unescape(match.group(1)), unescape(match.group(2))


def symbol_table():
    """command -> the code point the model inserts for it."""
    source = SYMBOLS.read_text(encoding="utf-8")
    return {unescape(name): int(code, 16)
            for name, code in SYMBOL_ENTRY.findall(source)}


def owned_code_points():
    covered = set()
    for table in TTFont(MATH_FONT)["cmap"].tables:
        covered |= set(table.cmap.keys())
    return covered


def test_every_face_glyph_is_in_the_embedded_math_font():
    covered = owned_code_points()
    missing = sorted(
        {(character, command)
         for command, face in cells()
         for character in face
         if ord(character) > 0x20 and ord(character) not in covered})
    assert not missing, (
        "%d palette face character(s) are absent from Latin Modern Math. "
        "pick_button_font tests the whole sample at once, so each one of "
        "these drops EVERY key in EVERY palette to a fallback font: %s"
        % (len(missing),
           ", ".join("U+%04X on %s" % (ord(character), command)
                     for character, command in missing)))


def test_a_symbol_key_shows_the_character_it_inserts():
    table = symbol_table()
    wrong = []
    for command, face in cells():
        payload = command[len("symbol."):] if command.startswith("symbol.") \
            else command
        payload = payload.strip()
        if payload in FACE_IS_A_WORD or payload not in table:
            continue
        if face != chr(table[payload]):
            wrong.append((payload, face, table[payload]))
    assert not wrong, (
        "%d key(s) draw a different character from the one they insert: %s"
        % (len(wrong),
           ", ".join("%s shows %s but inserts U+%04X"
                     % (command,
                        " ".join("U+%04X" % ord(c) for c in face),
                        code)
                     for command, face, code in wrong)))


if __name__ == "__main__":
    test_every_face_glyph_is_in_the_embedded_math_font()
    test_a_symbol_key_shows_the_character_it_inserts()
    print("PASS: palette keys show what they insert, in the shipped font")
