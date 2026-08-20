# EQNEDT64 — handover

Written 2026-08-20, at radia **4.95.56** (published to PyPI). Everything a new
session needs to continue the editor, including the parts that were expensive to
learn and would cost days to re-derive.

---

## 1. What this is, and the rule that decides arguments

A replacement for Microsoft Equation Editor 3.0 (EQNEDT32), shipped in the radia
wheel as `eqnedt64.exe` plus the `radia.equation` Python API.

**It has two references, and they are different documents:**

| | reference | how it is settled |
|---|---|---|
| **appearance** | TeX | measured against XeLaTeX box dimensions — numbers, never screenshots |
| **usability** | Equation Editor 3.0 | its own key table, read out of the binary |
| **file format** | neither | an equation is **LaTeX**, usually inside a Markdown file |

When a question is "should it look like this?", the answer comes from TeX. When
it is "should it behave like this?", the answer comes from EE3. The file format
follows neither, and that is deliberate: storing LaTeX is what lets the notes,
the papers, the slides and the MCP tooling all read the same equation.

Sugahara's standing weighting: **見た目:操作感 = 1:2.** The appearance half is
measured and largely settled. **The interaction half is where the remaining work
is.**

### Two notation rules (CLAUDE.md § Equation Notation)

1. **A vector is `\vec\bm`** — bold *italic* under the arrow. Applying the
   vector style applies both the arrow and the face.
2. **`\mathbf` is a different thing** — upright bold, for a matrix name. These
   were one typeface until the rule was written down.
3. **`\dfrac` is the default spelling of a fraction.** The outermost fraction is
   *drawn* at display size, so a bare `\frac` on the clipboard would come out
   smaller than the picture the author accepted. A **nested** fraction goes out
   as plain `\frac`, because LaTeX steps a nested one down by itself — the two
   rules then agree level by level. `\tfrac` is preserved when asked for.

These govern what the editor writes into files other people read. A silent
disagreement between the picture on screen and the LaTeX on the clipboard is the
same failure class as a silent fallback: the reader gets a result they cannot
audit.

---

## 2. File map

### C++ (`src/ext/equation/`)

| file | role |
|---|---|
| `tex_parser.cpp` | LaTeX → node tree |
| `latex_emitter.cpp` | node tree → LaTeX (owns the `\dfrac`/`\bm` rules) |
| `mtef_node.h/.cpp`, `mtef_common.h` | the node tree; MTEF record and typeface constants |
| `mtef_parser.cpp` | MTEF (`.eqn`) → node tree |
| `line_pass.cpp/.h` | the repair passes that put a torn-apart MTEF document back together |
| `mtef_omml.cpp` | → OMML (Word / PowerPoint / Excel native math) |
| `mtef_rtf.cpp` | → RTF (what Word actually accepts from the clipboard) |
| `mtef_mathml.cpp` | → MathML (PowerPoint, Excel) |
| `mtef_svg.cpp` | layout → SVG; the whole typesetting engine lives here |
| `math_font.cpp/.h` | reads the OpenType **MATH** table (TeX's `fontdimen` values) |
| `mtef_gdi.cpp` | GDI glyph metrics and EMF/PNG rendering |
| `tex2mtef.cpp`, `mtef2tex.cpp` | `.eqn` write / read helpers |
| `md_doc.cpp`, `md_blocks.cpp`, `md_layout.cpp` | which spans of a `.md` are math; block structure; document layout |
| `eq_edit.cpp/.h` | **the editing model** — caret, templates, selection, undo |
| `eq_chords.cpp/.h` | key chords → command names |
| `eq_window.cpp/.h` | **the window and the palette bar** (979 lines, no test) |
| `eqnedt64_main.cpp` | `WinMain` |
| `equation_pybind.cpp` | the pybind11 surface |
| `mtef_dump.cpp` | `dump_tree` / `tex_dump_tree` — print the parsed tree |

### Python (`src/radia/equation/`)

| | |
|---|---|
| `__init__.py` | re-exports the C++ types and the office writers |
| `office.py` | `markdown_to_docx`, `markdown_to_pptx`, `copy_to_clipboard`, `omml_paragraph`, `split_math` |

### Tests and validation

| | |
|---|---|
| `tests/equation/` | 33 files, **1388 pass, 1 xfail** |
| `validation_test/equation/tex_reference.tex` (+ `_matrix`, `_wide`) | the TeX side of the appearance check |
| `validation_test/equation/compare_boxes_with_tex.py`, `score_against_tex.py` | run the comparison |
| `validation_test/equation/convert_corpus.py` | `.eqn` corpus before/after harness |
| `validation_test/equation/read_ee3_key_table.py`, `read_ee3_strings.py` | EE3 reverse-engineering readers |
| `validation_test/equation/probe_window_chords.ps1` | drives the real window with keystrokes |

---

## 3. Build

```
cmake --build build-msvc --config Release --target _equation eqnedt64 -j
```

`C:\temp\build_eq64.bat` wraps that with `vcvars64.bat`. A post-build step copies
`eqnedt64.exe` and `_equation.pyd` into `src/radia/`.

**A header change needs a full rebuild.** A partial build links stale objects and
produces **333 failures with a `RuntimeError`** that looks exactly like a code
defect and is not one. If a change "breaks everything", rebuild before debugging.

---

## 4. How correctness is established

### Appearance — TeX numbers, not pictures

`tests/equation/test_tex_metrics.py` compares layout output against box
dimensions that XeLaTeX itself reports (`\the\wd0 \the\ht0 \the\dp0`) for
`validation_test/equation/tex_reference.tex`.

The reference **must** be XeLaTeX + `unicode-math` + `\setmathfont{Latin Modern
Math}`. `\usepackage{lmodern}` is the *wrong* reference — it selects Type1
lmex10 with different metrics, and a session already burned time chasing a
disagreement that was the reference's fault.

Typesetting parameters are read from the font's **MATH table**, not guessed:
`math_font.cpp` maps MathConstants onto the quantities TeX calls `fontdimen`.
Radicals follow TeX's own `make_radical`. Large operators, limits, fences,
accents and matrices each have their own locked comparison.

### MTEF reading — a real corpus, and a rule about accepting changes

```
python validation_test/equation/convert_corpus.py --health
python validation_test/equation/convert_corpus.py --diff
```

Runs the lab's `.eqn` corpus before and after a change and scores what is left.

**A change is accepted only when zero documents gain a defect.** Corpus health
went **620/779 → 774/779** clean under that rule. What remains: 5 stray style
markers (five *different* shapes) and 1 genuinely empty fence in a document that
was already wrong before this work.

Two warnings learned the hard way:

- **The health check itself lied once.** Its empty-fence regex counted *short*
  fences as *empty*, so 19 of 20 flagged documents were fine. It was fixed in
  its own commit (`b33c3a518`) because it was pointing the next piece of work at
  the wrong place. **Verify the checker before trusting the score.**
- **A defect-marker delta cannot see content.** Two files that recovered a whole
  missing integrand scored "neutral". Read the actual diff, not only the score.

### The debugging tool that settles arguments

```python
eq.dump_tree(data, run_passes=...)   # MTEF bytes → parsed tree
eq.tex_dump_tree(latex)              # LaTeX → parsed tree
```

Every disagreement between the LaTeX, OMML and SVG paths so far has come from
reading a *different field of the same tree*. Printing it settles the question in
seconds. Reach for this before theorising.

---

## 5. MTEF — the format knowledge

MTEF v3 is the Equation Editor 3.x / MathType binary format. It is supported in
**one direction**: to read equations out of documents that already contain them.
`tex_to_mtef` exists so an equation can be handed back as a `.eqn`. Fidelity of
the *write* path is not a goal.

**Header**: 5 bytes, `03 01 01 03 0a`.

**Records**: a tag byte whose low nibble is the type and high nibble the options.

| type | | | |
|---|---|---|---|
| 0 `END` | 1 `LINE` | 2 `CHAR` | 3 `TMPL` |
| 4 `PILE` | 5 `MATRIX` | 6 `EMBELL` | 7 `RULER` |
| 8 `FONT` | 9 `SIZE` | 10 `FULL` | 11 `SUB` |
| 12 `SUB2` | 13 `SYM` | 14 `SUBSYM` | |

Options: `OPT_LINE_NULL = 0x01`, `OPT_CHAR_EMBELL = 0x02`,
`OPT_LINE_LSPACE = 0x04`, `OPT_NUDGE = 0x08`.

### The discovery that made the reader work

**EQNEDT32 writes templates in pieces.** Display fractions, big-operator limits,
integral bodies and matrix rows are written **outside** their template,
separated by SIZE markers (`FULL` / `SUB` / `SYM`). A document is a **run of
top-level object lists**, not one list.

Two consequences that were each a separate bug:

- The parser must **loop over the whole stream**, not parse the first object
  list and stop. Guard the loop on "did `pos_` advance?" so a malformed stream
  cannot spin.
- A `MATRIX` header carries valign, hjust, vjust, rows, cols and then row/column
  partition arrays at **2 bits per line, packed** —
  `((n + 1) * 2 + 7) / 8` bytes each. Reading them as bytes crams every row into
  one cell.

### The repair passes

`line_pass.cpp` reassembles what the writer tore apart. Order matters and is
fixed in `PassPipeline::PassPipeline()`:

```
MatrixCellPass → IntegralSlotPass → DisplayFractionPass
              → FenceMergePass → BigOpDisplayPass → BigOpDisplayAltPass
```

The display fraction goes before the fence merge: it puts the fence back inside
the denominator, where the fence pass can then fill it from its own siblings.
`PassPipeline::process` additionally unwraps single-line PILEs, splices chunk
LINEs, trims dead trailing markers, drops dead SUB blocks, and recurses into
template slots.

### The rule every repair follows

**Refuse rather than guess.** No separator, no terminator, no closing character,
an empty denominator, a cell without line boundaries → **nothing is touched.**

> A fraction that is merely not repaired is visibly wrong.
> One repaired *wrongly* reads as finished.

Two traps found while writing these passes, both of which made output *worse*
than doing nothing:

- **Aliasing**: the slot being split *is* one of the fields written. Move the
  slot out first, or `slot.clear()` at the end erases the limit just set.
- **Choosing the wrong "deepest bare" big operator** made limits *vanish*. An
  operator with a block arriving in its own list is already spoken for — that is
  what `symFollows` checks.

---

## 6. Office output

| target | format | note |
|---|---|---|
| Word | **RTF** | Word's RTF spells OMML as control words: `<m:f>` is `{\mf` |
| Word / PowerPoint / Excel | **OMML** | `<m:oMath>` inside `<w:p>`; in PowerPoint inside `<a14:m>` inside `<a:p>` |
| PowerPoint, Excel | **MathML** | some paths read this instead |
| anywhere else | EMF / PNG | last resort |

`copy_to_clipboard` puts **every** format on the clipboard at once, so one Copy
serves whatever the target prefers.

**Why not a picture**: an equation pasted as OMML *is* an equation — Office's own
tools edit it, it follows the theme font and colour, it scales with the text, and
the reader needs nothing installed. A picture is none of those; an OLE object
needs Equation Editor on the reader's machine.

`count_equations` counts `<m:oMath[ >]` in `word/document.xml` /
`ppt/slides/slideN.xml` and is what the tests assert on.

### Slides: script-first Markdown

`markdown_to_pptx` builds a deck from one Markdown file:

| notation | meaning | goes to |
|---|---|---|
| `# heading` | starts a slide | slide title |
| `> quote` | **the script — what you say** | speaker notes |
| everything else | what the audience sees | slide body |

A slide nobody has scripted has **no notes page at all**, rather than an empty
one, so the coverage checks in `radia_mcp.presentation` can tell the difference.
Two layout rules are in the code because both were wrong in the first real deck:
a paragraph that is one display equation gets `buNone` and is centred, and the
body asks PowerPoint to shrink to fit rather than run off the page.

---

## 7. Python API surface

```python
import radia.equation as eq

eq.tex_to_omml(r"\frac{a+b}{c}")    # Office-native math
eq.tex_to_rtf(...) / tex_to_mathml(...) / tex_to_svg(...) / tex_to_emf(...) / tex_to_png(...)
eq.copy_to_clipboard(...)           # every format at once
eq.markdown_to_docx(md, "note.docx")
eq.markdown_to_pptx(md, "talk.pptx", title="…")
eq.mtef_to_latex(data) / mtef_to_omml(data) / tex_to_mtef(latex)
eq.read_eqn(path) / write_eqn(path, data)
eq.dump_tree(data) / tex_dump_tree(latex)
eq.tex_normalize(latex)             # LaTeX → tree → LaTeX fixed point
```

The editing model:

```python
e = eq.Equation()
e.insert_text("x"); e.command("template.sub"); e.insert_text("i")
e.latex()                            # 'x_{i}'
```

Current counts (verify rather than trust after a change):
**129 shortcuts, 61 templates, 10 symbol palette groups / 163 items,
9 template palette groups / 56 items.**

`Equation` methods worth knowing: `press`, `command`, `chord_steps`,
`caret_geometry`, `move_to_point`, `extend_to_point`, `selection_geometry`,
`selected_latex`, `insert_latex`, `set_style`/`styles`, `undo`/`redo`,
`svg`, `omml`, `shortcuts`, `templates`.

Undo snapshots the whole LaTeX rather than inverting commands: equations are
tens of characters, so a snapshot costs nothing and cannot drift out of step with
the tree. It rests on LaTeX → tree → LaTeX reaching a fixed point, which
`tex_normalize` exposes and `test_edit.py` checks.

---

## 8. CLOSED — the crash, and the interaction-layer tests that found it

Sugahara tried the editor on 2026-08-20 and it died on some operation
(`0xC0000005`, fault offset `0xc435`, ~52 s in). **Which operation was not
recorded**, there was no PDB, and no dump. It is fixed, and the way it was
found is the part worth keeping.

### The standard practice this project now follows

Nothing here is invented; it is the ordinary GUI-testing toolkit, and naming it
that way is deliberate — the lab had never shipped a GUI before, so the wheel
was there to be picked up rather than carved.

| practice | where it lives here |
|---|---|
| **Test pyramid** — most tests at the model level, few at the UI | 1397 model tests under `tests/equation/`; the window has one |
| **Humble object** — keep logic out of the untestable view | `eq_edit` / `eq_chords` decide everything; `eq_window` only reads modifiers and paints |
| **Programmatic driving, not record-and-replay** | `--selftest` sends window messages by handle; no keyboard, no foreground, so it runs beside a working user and on a headless desktop |
| **GUI ripping** — derive the cases from the widget tree | every chord comes from `chords()`, every cell from `Equation::*_palettes()`; nothing is written out twice |
| **Monkey / stress testing** (Android's Monkey is the canonical one) | seeded random walks over keys, chords, palette, mouse, wheel, resize, DPI |
| **Deterministic seeds** | xorshift64\*, seeds 1..N — a failing walk replays exactly |
| **A journal flushed before each step** | the last line names the step that killed the process |
| **Runtime memory verification** | `-DRADIA_EQ_ASAN=ON` builds the editor under AddressSanitizer |
| **Crash dumps** | WER `LocalDumps\eqnedt64.exe`, `DumpType=2`, into `C:\temp\wer_dumps\eqnedt64` |
| **Resource-leak oracle** | `GetGuiResources` sampled every 100 steps; a paint that leaks a GDI handle per frame fails |
| **Push the bug down the pyramid** | what the monkey found got a fast model-level regression test |

Not adopted, and why: **UI Automation / WinAppDriver** is the standard external
driver on Windows, but it needs a UIA provider and this window is custom-drawn
with no controls to expose. Message injection is the right level for it today;
if the editor ever needs to be driven by an outside tool, exposing UIA is the
move. **Visual/approval testing** of the painted output is not here either —
the appearance is already pinned at the model level by `test_tex_metrics.py`
against TeX's own box dimensions, which is stronger than a screenshot diff.

### Running it

```bash
eqnedt64.exe --selftest [--log <path>] [--walks N] [--steps N] [--clipboard]
```

Exit code is the failure count, or the exception code if a step crashed the
process. `--clipboard` opts into the copy/paste chords; without it the run
leaves the user's clipboard alone. `tests/equation/test_window_selftest.py`
wraps it for pytest.

Every chord × 4 caret states, every palette cell, the mouse over and outside
the canvas, blink/resize/minimize/DPI, then the walks. Roughly 100 s for
4 × 2000 steps.

### What it found, first run, 8.5 s in

A random walk died. The journal's last line said `chord template.nthroot`; the
dump symbolized to `take_selection`; and the ASan build named the line:

```
AddressSanitizer: access-violation ... READ
  #3 mtef::Equation::take_selection    eq_edit.cpp:692
  #4 mtef::Equation::insert_template   eq_edit.cpp:940
```

**The selection anchor was a bare index with no record of which slot it
indexed.** `clamp()` only ever clamped `index_`. So Tab, Ctrl+Up, a click or
undo — anything that moved to another slot or rebuilt the tree — left the
anchor pointing into a list it no longer belonged to; a shorter list then made
`take_selection()` read past the end of the vector and erase a range that was
not there. The heap was wrecked at that moment and the process died seconds
later somewhere unrelated, which is exactly the shape of the original report.

**The fix** (`eq_edit`): the anchor carries its slot (`anchor_path_`), and one
private `selection_range()` is the only place `lo`/`hi` are computed. A stale
anchor means **no selection** — not a selection quietly shrunk to fit, which
would delete a different range than the one highlighted.
`tests/equation/test_selection_anchor.py` locks it at model level.

Fixed in the same pass, both found by the sweep: **`WM_SYSKEYDOWN` was not
handled**, so every Alt chord (`Ctrl+Alt+Space`, the quad) was dead — the press
went to `DefWindowProc` and became menu activation. And **`Ctrl+Shift+S` was
published for `style.script` while the window uses it for Save As**, so that
chord could never reach the table; script moved to `Ctrl+Shift+P`.

---

## 9. OPEN — the palette is too busy

Sugahara, 2026-08-20: *ボタンのアイコンはもう少し簡素化してごちゃごちゃしない
がよい*.

Each button currently draws **its actual contents** — a real fraction with two
slots, a real matrix — rendered by the same layout engine that draws the
document. That property is deliberate and worth keeping: a hand-written table of
sample LaTeX would drift from what the templates really are, and a button would
start lying about what it inserts.

The density is a separate matter. Knobs, `eq_window.cpp:90-96`:

```cpp
kPaletteScale = 1.5      // everything scales together
kBtnW = 46, kBtnH = 24   // bar button, before scaling
kCellW = 34, kCellH = 30 // popup cell
kBarPt  = 9.0            // type size inside a bar button
kCellPt = 12.0           // type size inside a popup cell
```

A bar button shows *the first few members* of its group. Showing **one**
representative member, larger, would simplify the bar without giving up the
generated-from-the-real-template guarantee. Note `eq_window.cpp:213` — cells
already shrink-to-fit before centring, because centring alone let contents spill
over neighbours.

---

## 10. Other open usability items

- `PageUp` / `PageDown` (EE3 key-table **group 10**) and `Ctrl+Tab`
  (**group 7**) are still unbound. Decoding them needs the EE3 dispatcher at
  `0x427E37`.
- Only chords that are **unambiguous in EE3's own Help** are bound. The rest of
  `templates()` is reachable from the palette and is deliberately left unbound —
  do not invent chords.

---

## 11. Known limits — documented, not bugs to chase

- `\sum`'s stacked-limits flag has nowhere to live in LaTeX, so it does not
  survive a save. Display style stacks them anyway.
- An empty row has no LaTeX spelling, so a blank **trailing** row of `cases` or
  a matrix disappears on save. Filled rows round-trip exactly.
- A slashed fraction is not offered as a template: MTEF writes it as
  `{}^{a}/{}_{b}`, which does not read back as one fraction, and a template that
  changes shape when saved is worse than no template.
- Spacing commands (`\,` and friends) emit nothing; Office's own spacing is used.
- The SVG layout implements TeX's inter-atom spacing and per-glyph boxes, not the
  whole of TeX. Font metrics come from GDI, so the **layout path is
  Windows-only**; the other paths are not.
- Text-style fraction numerator sits **0.1437 pt** off TeX (xfail, with the
  arithmetic recorded in the test). `\overleftarrow` is **0.95 pt** off — TeX's
  own two arrows differ by 0.96 pt on the same body, so this is at the noise
  floor of the reference itself.

---

## 12. Rules that must not be broken

- **Never extract EE3's resources** (icons, bitmaps, art). This ships to a public
  GitHub repository and to PyPI. Observing and reimplementing *behaviour* is
  fine; copying its artwork is copyright infringement.
- **The lab `.eqn` corpus is private.** Write only derived LaTeX, to scratch
  paths under `C:\temp`. Commit none of it.
- **Do not run the keyboard probe without explicit consent**
  (`probe_window_chords.ps1` takes `-IAmNotUsingThisMachine`), and **never on
  100号機** — kubota and yano have live RDP sessions there.
- **Heredocs eat backslashes even when quoted.** Write patch scripts with the
  Write tool and build every backslash from `chr(92)`; `assert old in s` before
  replacing, and grep for a marker afterwards. This has silently produced
  no-op edits, and a `"\|"`-vs-`"\\|"` defect that swapped `\middle|` and
  `\middle\|` — round-tripping correctly the whole time, so only the *drawing*
  was wrong and nothing disagreed out loud.
- Pre-presentation material (manuscripts, slides, notes) lives in
  `W:\02_学会資料\`, never in this repository.

---

## 13. Acceptance targets — both met at 4.95.56

**(a)** A real Equation Editor document converts to `.tex` and can be edited in
EQNEDT64 to the same result. Checked on the 第8回 `.eqn`: every equation correct,
round-trip stable.

**(b)** A handwritten note (PDF) becomes correct TeX and from there a native
Microsoft (OMML) equation in PowerPoint. Seven pages of a lab notebook came
through as 72 equations, none failing to convert; PowerPoint opens and exports
the file without complaint.

These are acceptance targets for retiring EQNEDT32, not features. **Both are
met, and the crash that blocked retirement is fixed (§8) with the tests that
found it.** What remains before retirement is usability, not correctness.

---

## 14. Suggested order of work

1. ~~`/Zi` + `/DEBUG` + WER local dumps~~ — **done**, plus a `-DRADIA_EQ_ASAN=ON`
   build (§8).
2. ~~A test for `eq_window.cpp`, and the command-fuzz harness~~ — **done**:
   `--selftest` and `test_window_selftest.py`. It found the crash on its first
   run; the fix and two dead chords are in §8.
3. ~~Paste size in PowerPoint~~ — **done** (§15): the GVML clipboard format,
   pasting at 24 pt as a native equation.
4. **Palette simplification** (§9) — one representative member per button.
5. `PageUp`/`PageDown`, `Ctrl+Tab` (§10).
6. The 5 remaining stray-marker corpus documents (§4), lowest priority: they were
   wrong before this work and are five different shapes, so each is its own
   investigation.

---

## 15. CLOSED — a pasted equation is 24 pt

Sugahara, 2026-08-20: *powerpointに貼り付けたときは、24ptにしてほしいよ、18ptでは
小さい*.

**The clipboard could not say the size.** Measured, not assumed — PowerPoint
16.0 through its own object model, pasting into a real slide and reading the
run size back:

| clipboard MathML | size in PowerPoint |
|---|---|
| as we emit it | 28 pt |
| `<mstyle mathsize="24pt">` added | 28 pt |
| `mathsize` on `<math>` | 28 pt |

28 pt is that placeholder's own level-1 size. **PowerPoint ignores MathML
sizing entirely and uses the destination's**, which is why an 18 pt body box
gave 18 pt equations.

**PowerPoint's own copy says it.** Copying a shape out of PowerPoint puts
`Art::GVML ClipFormat` on the clipboard: an OPC package (a ZIP, `PK\x03\x04`)
holding `clipboard/drawings/drawing1.xml`, whose runs carry explicit
`<a:rPr sz="2400"/>`. So the route was the same one the RTF took —
**transcribe what the application itself puts on the clipboard.**

`gvml_clip.cpp` writes that package: a store-only ZIP (no deflate to carry for
a few kilobytes; PowerPoint accepts it — measured) holding `[Content_Types]`,
`_rels/.rels` and a `lockedCanvas` shape whose paragraph holds the `<a14:m>`
OMML. Two things kept it small:

- **The OMML writer did not have to change.** PowerPoint repeats `sz` on every
  `<m:r>`; a paragraph default (`<a:defRPr sz>` plus `<a:endParaRPr sz>`) was
  measured to give the same 24.0 pt, so `tex_to_omml` is used as it is.
- **No theme part.** PowerPoint's own package carries `theme1.xml` and a
  relationship to it; dropping both still pastes.

Offered **before** the MathML, since PowerPoint takes the richest format it
recognises. `mtef::kPasteSizePt` is 24 — 18 is small to read from the back of a
room, and PowerPoint's own default shrinks from 18 down by outline level.
Exposed as `radia._equation.tex_to_gvml(latex, size_pt, display)` and
`PASTE_SIZE_PT`, and used by both clipboard paths (the editor's Ctrl+C and
`radia.equation.office.copy_to_clipboard`, which grew a `size_pt` argument).

`tests/equation/test_paste_size.py` checks the package without Office, and its
last test pastes into a real slide and asserts 24.0 pt — the claim is about
PowerPoint, so something has to ask PowerPoint.
