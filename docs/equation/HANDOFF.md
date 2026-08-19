# EQNEDT64 — handoff

Written 2026-08-20, at radia 4.95.56 (published). This is what another session
needs to pick the editor up, and what is still open.

## What it is

A replacement for Microsoft Equation Editor 3.0, shipped in the radia wheel as
`eqnedt64.exe`. It has **two references, and they are different documents**:

- **appearance follows TeX** — measured, in numbers, not by eye
- **usability follows Equation Editor 3.0** — its own key table, read out of the
  binary rather than guessed at
- **the file format follows neither** — an equation is stored as **LaTeX**,
  usually inside a Markdown file, which is what lets everything else read it

Two notation rules are decided (CLAUDE.md § Equation Notation): a vector is
`\vec\bm` (bold *italic*); `\mathbf` is upright bold and a different thing;
`\dfrac` is the default spelling of a fraction, and a nested one goes out as
plain `\frac` so the two rules agree level by level.

## Where the code is

| | |
|---|---|
| `src/ext/equation/tex_parser.cpp` | LaTeX → node tree |
| `src/ext/equation/latex_emitter.cpp` | node tree → LaTeX |
| `src/ext/equation/mtef_parser.cpp`, `line_pass.cpp` | MTEF (`.eqn`) → tree, and the repair passes |
| `src/ext/equation/mtef_omml.cpp` / `mtef_rtf.cpp` / `mtef_mathml.cpp` | Office clipboard formats |
| `src/ext/equation/mtef_svg.cpp`, `math_font.cpp` | layout; MATH-table typesetting parameters |
| `src/ext/equation/eq_edit.cpp` | the editing model (caret, templates, undo) |
| `src/ext/equation/eq_chords.cpp` | key chords → command names |
| `src/ext/equation/eq_window.cpp` | the window and the palette bar |
| `src/radia/equation/office.py` | Word / PowerPoint writers, `markdown_to_pptx` |

## Build

```
cmake --build build-msvc --config Release --target _equation eqnedt64 -j
```

`C:\temp\build_eq64.bat` wraps that with `vcvars64.bat`. The post-build step
copies `eqnedt64.exe` and `_equation.pyd` into `src/radia/`.

**A header change needs a full rebuild.** A partial one links stale objects and
produces 333 failures with a `RuntimeError` that looks like a code defect and is
not one.

## How it is checked

- `tests/equation/` — 1388 pass, 1 xfail. `test_tex_metrics.py` is the
  appearance half: it compares against XeLaTeX box dimensions produced by
  `validation_test/equation/tex_reference.tex`. **Numbers, not screenshots.**
- `validation_test/equation/convert_corpus.py` — runs a real `.eqn` corpus
  before and after a change and scores what is left. A change is accepted only
  when **no document gains a defect**. The lab corpus is private; only derived
  LaTeX is written, to scratch paths, and none of it is committed.
- Corpus health went 620/779 → 774/779 clean. What remains: 5 stray style
  markers (five different shapes) and 1 genuinely empty fence in a document
  that was already wrong.

## OPEN: it crashes, and we cannot yet say where

Sugahara hit a crash on 2026-08-20 while trying the editor. **Which operation
was not recorded.**

```
Application Error, eqnedt64.exe
exception 0xC0000005 (access violation)
fault offset 0x000000000000c435   (RVA in eqnedt64.exe itself)
~52 s after launch
```

**It cannot be symbolized today.** `CMakeLists.txt:1112` builds the target with
`/O2 /W3 /utf-8` — no `/Zi`, and no `/DEBUG` on the link — so no PDB exists, and
the offset cannot be turned into a function. There was no minidump either;
WER local dump collection is not enabled for this image.

**Do this first, before trying to reproduce:**

1. Add `/Zi` to the target's compile options and `/DEBUG` to its link options
   (keep `/O2`; RelWithDebInfo-style, so the shipped binary stays optimized).
   Make sure the `.pdb` lands next to the `.exe`.
2. Enable local dumps for the image:
   `HKLM\SOFTWARE\Microsoft\Windows\Windows Error Reporting\LocalDumps\eqnedt64.exe`
   with `DumpType=2` (full).
3. Then ask Sugahara which operation it was, or fuzz the command surface —
   `Equation::commands()` is enumerable, so a harness can drive every command
   from several caret positions and find it without a human.

The editing model is covered by tests at the API level (`test_edit.py`,
`test_key_dispatch.py`, `test_usability.py`), so the fault is most likely in the
window/painting layer, which those tests do not touch — `eq_window.cpp` is 979
lines and has no test of its own. That is the gap to close.

## OPEN: the palette icons are too busy

Sugahara, 2026-08-20: *ボタンのアイコンはもう少し簡素化してごちゃごちゃしないがよい*.

The current design (`eq_window.cpp:75-230`) draws each button's **actual
contents** — a real fraction with two slots, a real matrix — rendered by the
same layout engine that draws the document, so a button cannot drift from what
inserting it produces. That property is worth keeping; the visual density is
not. The knobs:

```cpp
kPaletteScale = 1.5      // everything scales together
kBtnW = 46, kBtnH = 24   // bar button, before scaling
kCellW = 34, kCellH = 30  // popup cell
kBarPt  = 9.0             // type size inside a bar button
kCellPt = 12.0            // type size inside a popup cell
```

A button currently shows *the first few members* of its group. Showing **one**
representative member, larger, would simplify without giving up the
generated-from-the-real-template property.

## Other open usability items

- `PageUp` / `PageDown` (EE3 key-table group 10) and `Ctrl+Tab` (group 7) are
  still unbound. Decoding them needs the EE3 dispatcher at `0x427E37`.
- Sugahara's standing judgement: **usability is the remaining focus**
  (改善比率 見た目:操作感 = 1:2). The appearance half is measured and largely
  settled; the interaction half is where the work is.

## Known limits (documented, not bugs to chase)

- `\sum`'s stacked-limits flag has nowhere to live in LaTeX, so it does not
  survive a save. Display style stacks them anyway.
- An empty row has no LaTeX spelling, so a blank trailing row of `cases` or a
  matrix disappears on save. Filled rows round-trip exactly.
- A slashed fraction is not offered as a template: MTEF writes it as
  `{}^{a}/{}_{b}`, which does not read back as one fraction, and a template that
  changes shape when saved is worse than no template.
- The SVG layout implements TeX's inter-atom spacing and per-glyph boxes, not
  the whole of TeX. Font metrics come from GDI, so the layout path is
  Windows-only; the others are not.
- Text-style fraction numerator sits 0.1437 pt off TeX (xfail, arithmetic
  recorded in the test). `\overleftarrow` is 0.95 pt off — TeX's own two arrows
  differ by 0.96 pt on the same body.

## Rules that must not be broken

- **Never extract EE3's resources** (icons, bitmaps). This ships to a public
  GitHub repo and PyPI. Observing and reimplementing *behaviour* is fine;
  copying its art is not.
- The lab `.eqn` corpus is private. Write only derived LaTeX, to scratch paths.
  Commit none of it.
- Heredocs eat backslashes even when quoted. Write patch scripts with the Write
  tool and build every backslash from `chr(92)`; assert the search text is
  present before replacing. This has silently produced no-op edits and a
  `"\|"`-vs-`"\\|"` defect that swapped two delimiters for weeks.

## Acceptance targets (both met at 4.95.56)

(a) A real Equation Editor document converts to `.tex` and can be edited in
EQNEDT64 to the same result — checked on the 第8回 `.eqn`.
(b) A handwritten note (PDF) becomes correct TeX and from there a native
Microsoft (OMML) equation in PowerPoint — seven notebook pages, 72 equations,
none failing to convert.
