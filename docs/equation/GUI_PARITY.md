# EQNEDT64 vs Equation Editor 3.0: what the interface has to do

Equation Editor 3.0 has no published source — it is a cut-down MathType,
licensed and proprietary — so this inventory comes from the running program.
Its menu tree was walked with `GetMenu` / `GetSubMenu` / `GetMenuStringW`,
which is what a person does by opening every menu, done exhaustively instead of
by memory. Observing an interface and reimplementing it is the same thing we
already did for the key chords and the type sizes; nothing here is taken from
the binary.

The point of writing it down is that the alternative is discovering the gaps
one at a time while debugging, which is expensive and leaves holes.

## Two references, and which one settles a question

There are two models here, and they answer different questions.

**Appearance follows TeX.** Type sizes, the space around an operator, where a
fraction sits about the baseline, how tall a radical grows -- these go to TeX.
Nearly everyone who reads an equation has read far more mathematics set by TeX
than by Equation Editor, so TeX's look is the one that reads as ordinary.

The font is part of that: Latin Modern Math, which is Computer Modern in
OpenType.  Its MATH table does not approximate TeX's parameters, it *is* them
-- TeX reports a display numerator shift of 8.12389 pt at 12 point and the
table says 677 per thousand, which is 8.124.  Setting in it uses the same
numbers TeX lays out from rather than imitating the result.

**Usability follows Equation Editor 3.0.** Key chords, what the caret does,
how the palettes are arranged, that Ctrl+S saves -- these come from the program
this replaces, because that is what the hands using it already know.

**The file format follows neither.**  MTEF is Equation Editor's storage, not
its behaviour, and nothing here needs to be shaped by it.

Where the two references disagree, the split above decides.  Two worked
examples:

* Equation Editor leaves one point of bar past a fraction's contents
  ("Fraction bar overhang"), and the OpenType MATH table has no such notion.
  Appearance -- but TeX has no opinion either, so Equation Editor's is kept.
* A fraction inside a fraction is set smaller by TeX and not by Equation
  Editor, which instead offers two templates and lets the writer choose.  That
  design cannot cross over: the input here is LaTeX, where `rac` is one
  command, so the size must come from the nesting.  TeX's rule.

Differences are settled by measurement, not by eye.  TeX will state the width,
height and depth of any box it sets, and `radia.equation.tex_metrics` returns
the same three numbers, so the comparison is arithmetic and needs no window --
which also means it can be re-run from a script.

## The complete menu, as observed

    File    Save (Ctrl+S) | Exit
    Edit    Undo (Ctrl+Z) | Cut (Ctrl+X) | Copy (Ctrl+C) | Paste (Ctrl+V)
            | Clear (Delete) | Select All (Ctrl+A)
    View    100% | 200% | 400% | Zoom... | Toolbar | Redraw | Show All
    Format  Align Left | Center | Right | at = | at .  | Matrix... | Spacing...
    Style   Math | Text | Function | Variable | Greek | Matrix-Vector
            | Other... | Define...
    Size    Full | Subscript | Sub-Subscript | Symbol | Sub-Symbol
            | Other... | Define...
    Help    Help Topics (F1) | About

Two things are worth noticing before the gap list.

**Size** is exactly the five numbers `SvgStyle` carries — Full 12, Subscript 7,
Sub-Subscript 5, Symbol 18, Sub-Symbol 12. That is where they came from, and it
confirms the ratios are the editor's own rather than a guess.

**Redraw** is the one item deliberately not reproduced. It exists because the
old editor's incremental display could go stale; every paint here recomputes
the layout from the tree, so there is nothing for it to repair.

## Where EQNEDT64 stands

| Equation Editor | EQNEDT64 | |
|---|---|---|
| Undo / Redo | model has both, Ctrl+Z / Ctrl+Y wired | done |
| Copy | Ctrl+C, full Office payload | done |
| Clear (Delete) | Delete / Backspace, backspace unwraps | done |
| templates, symbols | 26 templates, 163 symbols, palette bar | done |
| Show All | slot boxes always drawn | partial — not a toggle |
| Redraw | absent by design | n/a |
| **Select All** | — | **missing** |
| **Cut** | — | **missing** |
| **Paste** | — | **missing** |
| **100/200/400% + Zoom** | fixed 2x | **missing** |
| **Style: Math/Text/Function/Variable/Greek** | tree has the typefaces, no way to set them | **missing** |
| **Style: Matrix-Vector** | same — this is the bold-vector style | **missing** |
| **Size: Full…Sub-Symbol** | `SvgStyle` holds them, no way to change them | **missing** |
| **Format: align left/centre/right/at =/at .** | — | **missing** |
| **Format: Matrix…** | fixed 2x2 and 3x3 templates only | **missing** |
| **Format: Spacing…** | TeX inter-atom spacing, not adjustable | **missing** |
| **Save** | — | **missing** (`.tex`, UTF-8) |
| Toolbar toggle | palette bar always shown | low value |
| Help / About | — | low value |

Nothing in the "done" column is in doubt; the missing rows are the work.

## Order of work, and why

**1. Selection.** `Select All`, `Cut`, and half of `Clear` all rest on it, and
without it a range cannot be deleted or replaced — which is most of editing.
The model is a range within ONE slot: a template is a single item in its
parent's slot, so selecting a whole fraction is the same operation as selecting
a run of characters, and no model spanning tree levels is needed.

Brings: Shift+arrows, Shift+click, drag, Ctrl+A, typing over a selection,
Ctrl+C of a partial selection, Ctrl+X.

**2. Paste.** Ctrl+V, reading the clipboard's LaTeX. Today there is no way to
bring an equation in from anywhere else, which makes the editor a one-way door.

**3. Zoom.** Mouse wheel plus the 100/200/400 steps. Small, and the fixed 2x is
already wrong on a large monitor.

**4. Style.** The typefaces exist in the tree and nothing sets them. This is
what a bold vector needs, and vectors are most of what this lab writes. Chords
matter more than a menu here, since there is no menu bar.

**5. Size.** `SvgStyle` is already the contract; this is a dialog over numbers
that already work, so it is cheap once there is somewhere to put it.

**6. Format.** Alignment at `=` is the one that matters for a derivation across
several lines. Matrix... and Spacing... are further down.

**7. Save.** `.tex`, UTF-8, the LaTeX as it stands. Round-tripping is already
tested, so this is file dialogs.

**Show All** becomes a toggle when there is somewhere to put it; the boxes
themselves are done.

## What is NOT reproduced, deliberately

- **Redraw** — see above.
- **A menu bar** — Edit and Help were both declined; the palette bar is the
  visible interface and the rest is chords.
- **`.eqn` / MTEF output** — the lab has no `.eqn` assets, does not edit old
  equations, and Office displays the ones already embedded in documents.
