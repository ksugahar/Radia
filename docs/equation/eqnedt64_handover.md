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
| `eq_window.cpp/.h` | **the window and the palette bar**, and `--selftest`, which drives it (§8) |
| `gvml_clip.cpp/.h` | → `Art::GVML ClipFormat`, the only clipboard format that can state a paste size (§15) |
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
| `tests/equation/` | 36 files, **1409 pass, 1 xfail** |
| `tests/equation/test_window_selftest.py` | runs `eqnedt64.exe --selftest`; deselected from the quick suite because it takes ~100 s |
| `validation_test/equation/tex_reference.tex` (+ `_matrix`, `_wide`) | the TeX side of the appearance check |
| `validation_test/equation/compare_boxes_with_tex.py`, `score_against_tex.py` | run the comparison |
| `validation_test/equation/convert_corpus.py` | `.eqn` corpus before/after harness |
| `validation_test/equation/read_ee3_key_table.py`, `read_ee3_strings.py`, `read_ee3_menus.py` | EE3 reverse-engineering readers |
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
went **620/779 → 774/779 → 779/780** clean under that rule. **No document
carries a stray style marker or an empty fence any more.** One document is left,
and its defect class was there all along — the checker could not see it (below).

**The five stray markers were NOT five different shapes, and the documents were
not "already wrong".** That reading stood for months and was wrong on both
counts. Four of the five came from one source (Harrington ch2) with one
signature, and Equation Editor renders them correctly — so the fault was the
reader's. Two are fixed (§18); the rest are described there.

Two warnings learned the hard way:

- **The health check itself lied once.** Its empty-fence regex counted *short*
  fences as *empty*, so 19 of 20 flagged documents were fine. It was fixed in
  its own commit (`b33c3a518`) because it was pointing the next piece of work at
  the wrong place. **Verify the checker before trusting the score.**
- **A defect-marker delta cannot see content.** Two files that recovered a whole
  missing integrand scored "neutral". Read the actual diff, not only the score.
- **A score says a document carries a defect; it never says whose fault it is.**
  Render the document in Equation Editor and look —
  `validation_test/equation/render_in_ee3.ps1`. Five documents were written off
  as malformed input until EE3 drew the first of them perfectly (§18).
- **The checker missed half a fraction.** It looked for an empty NUMERATOR and
  not an empty denominator, so a document reading `a/()` — with its denominator
  escaped from the fraction entirely — scored clean. The corpus briefly read
  **100 %** with that document still wrong. The check was added; the honest
  score is 778/780. **A run at 100 % is a reason to look harder at the
  checker, not to stop.**

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

## 9. CLOSED — the palette was too busy (see §17)

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

## 10. DECIDED — the three keys that stay unbound

`PageUp`, `PageDown` and `Ctrl+Tab` were carried as an open item ("decoding them
needs the EE3 dispatcher at `0x427E37`"). They are now **closed as decided, not
as done**: this editor has nothing for them to do.

**What the evidence says.** The key table's navigation records (kind 6) group
cleanly, and reading them all at once is what makes the gaps legible:

| group | members | what it is |
|---|---|---|
| 0 | `Shift+Tab`→1, `Tab`/`Insert`→3 | previous / next **slot** |
| 1 | `Left`→1 `Right`→3 `Up`→6 `Down`→7 | move by **character** |
| 2 | the same four with `Shift` | **select** by character |
| 3 | `Ctrl+Shift+Left`→1, `Ctrl+Shift+Right`→3 | **select** by item |
| 6 | `Enter`→13 | break the line |
| **7** | **`Ctrl+Tab`→0** | **only member** |
| 8 | `Ctrl+Left`→1 `Ctrl+Right`→2 `Ctrl+Up`→3 `Ctrl+Down`→4 | move by **item** |
| **10** | **`PageUp`→2, `PageDown`→3** | **only members** |
| 11 | `End`→0, `Home`→1 | ends of the slot |

Everything else in that table is bound here. The two odd groups are not.

**Ctrl+Tab is a TAB STOP, and this editor has no tab stops.** The string table
holds `Tab symbol`, `Tab Stop Changes` (an undo description) and *"Using tab
formatting without left alignment may produce unexpected results"* — EE3 aligns
with tab stops, and since plain `Tab` already moves between slots, the chord for
inserting one is the obvious `Ctrl+Tab`. Alignment here is `format.left` /
`center` / `right` / `at_eq`, applied to a line; there is no tab stop for a key
to insert.

**PageUp / PageDown want something to page through.** There is none: one
equation, always fully visible — which stopped being an aspiration and became
true in §16, and is the same reason the mouse wheel is free to zoom rather than
scroll.

**How this was settled, so nobody re-derives it.**
`validation_test/equation/read_ee3_menus.py` dumps every EE3 menu with the
accelerator text it shows the user: **no menu mentions any of the three**, so
they are not commands a user could find and name — they are internal navigation.
(The Help file `EQNEDT32.HLP` would be the other authority; it is a
phrase-compressed WinHelp file that this session did not decode. The menus were
enough: everything EE3 offers by name is there, and these are not.)

An observation probe was built and then set aside: EE3 consults its key table in
the **message loop**, so an injected `SendMessage` bypasses it entirely, and
`PostMessage` plus `AttachThreadInput` is needed before a modifier registers at
all — `SetKeyboardState` writes the calling thread's input state, which a window
in another process never reads. It got that far; it did not get far enough to
trust what a keypress did. If the question reopens, start there, and hold the
modifier **across** the post: clearing it immediately races, and `Ctrl+F` arrives
as a plain `f`.

**The standing rule is unchanged**: only chords that are unambiguous in EE3's own
documentation are bound. The rest of `templates()` is reachable from the palette
and is deliberately left unbound — do not invent chords. Binding these three to
something invented would give the editor keys that no EE3 muscle memory expects
and no feature here supports.

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
4. ~~The window at awkward sizes~~ — **done** (§16): it fits what it holds and
   grows to what it holds.
5. ~~Palette simplification~~ — **done** (§17), along with two drawing bugs it
   exposed.
6. ~~`PageUp`/`PageDown`, `Ctrl+Tab`~~ — **decided** (§10): they stay unbound,
   because nothing here is a tab stop and there is nothing to page through.
7. ~~The remaining stray-marker corpus documents~~ — **none left** (§18). Two
   documents still lose a fraction's denominator; that class is described there.

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

**But that test placed the package on the clipboard itself, and the editor's
Ctrl+C was broken for a day without it noticing — see §21.**

---

## 16. CLOSED — the window fits what it holds, and grows to what it holds

Sugahara, 2026-08-20: *窓サイズを自動的に内容に調整するのか？小さいときには
み出さないが大事*.

It did neither. The equation was drawn at the asked-for zoom whatever the
window was, from a fixed left margin and centred vertically — so a long
equation ran off the right edge, a tall one painted over the palette bar and
the status strip, and there was no scrollbar and no way to bring either back.
The comment on the mouse wheel even said *"there is nothing to scroll — one
equation, always fully visible"*, which was the intention and not the truth.

**Two changes, in that order of importance.**

**It cannot overflow.** `view_of()` computes the canvas — between the palette
bar and the status strip — and drops the drawing scale until the equation fits
inside it. The margin is given up first, because on a narrow window the white
space matters less than the maths. The status strip shows the scale actually on
screen, with `(fit)` when it is not the one that was asked for: a window that
silently disobeys the zoom setting would be a mode you cannot see, which is the
same trap the `Style:` cell exists to avoid. A clip rectangle is set over the
canvas as well, so no layout can paint over the buttons even if its extents
were wrong.

`view_of()` is also the one place the position and scale are computed, and the
painter and the mouse now both call it. The formula used to be written out
three times, and only one copy needed to drift for a click to land somewhere
other than where the caret appeared.

**It grows.** `fit_window_to_content()` enlarges the window so the equation
fits at the zoom asked for, capped at the monitor work area. It only ever
grows, never shrinks: a window that got smaller as you deleted would move under
the hand that is typing. And it stops the moment the user takes hold of an edge
(`WM_ENTERSIZEMOVE` sets `user_sized`) — a window someone has deliberately
sized is theirs.

**Both are checked by the self-test, every paint.** `check_inside()` runs after
each of the thousands of paints in the sweeps and walks — including at 90×40,
which nobody would choose and which is exactly where a rule that only holds at
the default size would go unnoticed. `sweep_autosize()` checks the window grows
for an equation wider than it (760 → 1855 px) and then holds still once the
user has sized it.

The oracle was checked against a build with the fitting disabled, because a
test that cannot fail proves nothing: it reported 6 overflows and exit 5. That
negative control is not in the tree — run it by adding a `v.fit = 1.0;` after
the two fit lines in `view_of()`.

---

## 17. CLOSED — one icon per button, and the two bugs that hid behind three

Sugahara, 2026-08-20: *ボタンのアイコンはもう少し簡素化してごちゃごちゃしない
ほうがよい*, and then, seeing the result, *ステータスバーがきちんとみえない*.

Each bar button drew the **first three members** of its group. Three real
templates, each with its own empty slots, shrunk into a 69×36 button came out as
a smudge: the fences button was `([{ }])` drawn on top of itself and the
integrals were indistinguishable from the sums.

**The fix is fewer things, not smaller ones.** One representative member, drawn
at the size the popup already uses for its cells. A button still wears a REAL
member rendered by the real layout, so it cannot start advertising something the
template no longer is — that property is the reason there is no hand-written
table of sample LaTeX, and it is kept.

**Which member.** "First in the list" is Equation Editor's ordering, and it is
not always the member that says what the group IS: matrices opened on a 1×2
(two boxes), relations on "approximately", arrows on a double LEFT arrow.
`PaletteGroup::representative` names one, `PaletteGroup::icon()` resolves it,
and a group that names nothing wears its first member. The named ones are
`\neq`, `\rightarrow`, `\forall`, `\subset`, `\infty`, `vec`, `matrix2x2`.

An icon's empty slots are drawn **wider** than the editor draws them
(`kBarSlotEm = 1.6` against the editor's `0.55`): a decoration is only as wide
as what it decorates, so an overbrace over an editing-width slot was six points
across. It changes nothing about what the button inserts.

### The two bugs this uncovered

**A palette icon was painting into the row of buttons above it.**
`layout_brace_deco` placed an overbrace by shifting it up by the brace's
**ascent** where it should have used its **descent**, leaving the glyph a whole
brace-height too high — outside the box the layout then reported. Consequences
beyond the palette: the picture path sizes its bitmap from that box, so
**`\overbrace` was cropped out of every pasted picture entirely** — a blank strip
of paper where the annotation should be. The underbrace branch had the
arithmetic right all along, which is why only one of the two was ever wrong.

Caught now by `tests/equation/test_ink_inside_box.py`, which renders a construct
and counts ink: a decoration must put something in the third of the picture it
belongs in. It fails on the old build and passes on the new one.

**The status strip was sliced in half by its own border.** `draw_layout` sets
`TA_BASELINE` — right for a laid-out glyph, and it never put it back.
`DrawTextW` is documented to require `TA_TOP`, and with a baseline left behind
it silently draws a whole line higher, so "Style: Math" straddled the separator
with its lower half hidden. A leaked DC setting, one caller away from where it
was set. `draw_layout` now restores the text alignment and background mode it
changes, and `paint_status` states what it needs rather than inheriting it.

While there: the **zoom moved to the right edge** of the strip, where Word,
Excel, PowerPoint, Acrobat and every browser put it. It reads the scale actually
on screen and marks it `(fit)` when the window, not the setting, is what shrank
the equation (§16).

### How the icons were checked

`PaletteGroup.icon` is exposed to Python so `tests/equation/test_palette_icons.py`
can assert what the **window** resolves. Membership alone would have proved
nothing: the first version of the icon table reached the compiler as `"\neq"`
with ONE backslash — C++ read it as a newline followed by `eq` — so every symbol
icon silently fell back to the member it was meant to replace. Nothing failed;
the buttons simply did not change.

`eqnedt64 --selftest` also journals what each button wears
(`[bar N] <group> icon=<command> glyphs=… w=… asc=… desc=…`). That dump is what
identified the stray arc as the Braces glyph rather than anything to do with the
Logic button it was drawn in — the model, the SVG and the picture paths all
agreed and all looked innocent, because the drawing was simply somewhere else.

---

## 18. CLOSED (mostly) — the corpus documents that were never malformed

Five documents carried a stray `\scriptstyle`, and the note here said they "were
wrong before this work" and were "five different shapes, so each is its own
investigation". Both halves were wrong.

**Four of the five are one shape**, all from Harrington ch2, all with the same
signature: an integral whose limits arrive as a size-wrapped block AFTER the
line the integral sits in.

**And the documents are fine.** `render_in_ee3.ps1` opens one in Equation Editor
and captures what IT draws:

```
V = ∫₋ₐ^a dx′ ∫₋ₐ^a σ(x′,y′) / (4πε₀ √((x−x′)² + (y−y′)²)) dy′
```

A correct double integral. The reader was at fault, not the input. That is the
lesson worth keeping: **a corpus score says a document carries a defect marker;
it never says whose fault the defect is.** Nothing in the harness could have
told the difference, and nobody asked the one program that knows.

### What was wrong

Equation Editor writes an operator's limits as a block AFTER the thing they
belong to — `SIZE SUB`, the limit lines, `SIZE SYM`, the operator's glyph — and
`BigOpDisplayPass` reunites them. It reaches the owner two ways:

- **as a sibling**: the reverse scan meets an integral and takes it, no questions;
- **inside a LINE**: `deepestBareBigOp` looks for one that has "no limits yet".

The second asked `!hasLower && !hasUpper`. But those record which SLOTS the
template wrote, and an integral written with **variation 2 — the common case —
puts its INTEGRAND in the slot this reader calls `upper`**. So every ordinary
integral looked like it already had limits, and one sitting inside a line could
never be handed the block that followed. The two paths now ask the same
question: has a display block been given to it, and is one coming later in its
own list.

Result, at zero cost to anything else: **778 of 780 documents byte-identical**,
2 recovered, 0 regressions.

```
before: V = \int ^{dx'}\scriptstyle  - aa\int \displaystyle \int \limits_{-a}^{a}4\pi ...
after : V = \int \limits_{ - a}^{a}\int \limits_{ - a}^{a}4\pi ...
```

### Where the regression lock lives, and why not in tests/

`validation_test/equation/test_nested_integral_limits.py`, gated on the corpus
being present (`RADIA_EQN_CORPUS`, or the lab path).

A `tests/` fixture was written first, using `tex_to_mtef` to build the shape and
read it back. **It passed against the broken reader**, because our writer emits
the tidy form where limits are real slots — the torn-apart shape only exists in
documents EQNEDT32 itself wrote. It was deleted rather than kept: a test that
cannot fail claims a coverage it does not have, which is worse than no test.

### What is left (3 documents)

- `harrington_ch2_operator_l_plate.eqn` — same family, and it also loses the
  fraction; the limit block is not the only thing adrift.
- `harrington_ch2_polarizability_xx_matrix.eqn` — carries the empty fence as
  well: `\left( \left. \right\rangle \right)`.
- `perturbation_alpha_native.eqn` — a different shape: a fraction with an empty
  denominator under a pile of four size markers.

Start each by rendering it in EE3.

---

## 19. CLOSED — nudge is two bytes, and the reader was skipping four

The last stray-marker document did not look like the other four. Its tree held a
`CHAR` with **typeface -128 and code 0x8083**, and the equation's actual content
— both inner products of a fraction — was nowhere in it at all.

Those numbers are the tell. `0x80` is a nudge byte and `0x83` is a record tag:
the reader was standing two bytes off and manufacturing characters out of
whatever it landed on. **`skipNudge()` skipped four bytes. A nudge is two** —
`dx` and `dy` as signed bytes, with `(-128, -128)` as an escape meaning two
16-bit values follow instead. Every nudged record therefore ate the first two
bytes of what came next and desynchronised the rest of the stream.

Nudged records are rare, which is why this survived 779 of 780 documents.

```
before: \alpha = 1 - \dfrac{\left. \scriptscriptstyle \scriptstyle
                             \scriptscriptstyle \displaystyle \right\rangle }{}
after : \alpha = 1 - \dfrac{\left\langle f_{0},g \right\rangle }{}
                     \left\langle f_{0},Mf_{0} \right\rangle
```

Measured over the corpus: **778 of 780 byte-identical, 2 recovered, 0
regressions.** The other recovered document is the polarizability matrix, which
had been carrying the corpus's only empty fence.

This is the failure mode worth remembering: **a desynchronised reader does not
crash and does not produce garbage.** It produces well-formed LaTeX for a
different equation. Nothing downstream can tell. The only signals were a
typeface that cannot exist and content that never appeared.

### Where the score went to 100 %, and why that was the warning

With this fixed, the corpus reported **780/780, no defect marker**. It was not
true. `perturbation_alpha_native.eqn` still reads

```
\dfrac{\left\langle f_{0},g \right\rangle }{}\left\langle f_{0},Mf_{0} \right\rangle
```

— an empty denominator, with the denominator's content sitting *after* the
fraction. `DisplayFractionPass` never put it back. The checker had a rule for an
empty **numerator** (`\dfrac{}`) and none for an empty denominator, so the
document scored clean.

The rule was added, and the honest score is **778/780**. Two documents lose a
denominator this way; both did so before any of this work, and neither is a
regression.

**Read a 100 % as "check the checker".** It has now been wrong twice: once
counting short fences as empty (§4), once unable to see half a fraction.

### What was left, and what closed it

`perturbation_alpha_native.eqn` turned out to be a **third shape** of display
fraction: numerator inside the template, denominator **immediately after it**,
with no size markers between. The two shapes `DisplayFractionPass` already knew
both need a separator and a terminator to find where the denominator ENDS, and
refuse without them — right, when the boundary is a guess. Here it is not: a
LINE is one self-contained chunk, and the denominator is exactly that chunk. The
rule is narrow (denominator empty, very next live sibling a real LINE) and moved
one document:

```
before: \dfrac{\left\langle f_{0},g 
ight
angle }{}\left\langle f_{0},Mf_{0} 
ight
angle
after : \dfrac{\left\langle f_{0},g 
ight
angle }{\left\langle f_{0},Mf_{0} 
ight
angle }
```

which is what Equation Editor draws. 779 of 780 byte-identical, 0 regressions.

---

## 20. OPEN — one document, and why the safe rules are exhausted

`harrington_ch2_operator_l_plate.eqn` is the last corpus document that is wrong,
and it is wrong three ways at once. Equation Editor draws:

```
L(f) = ∫₋ₐ^a dx′ ∫₋ₐ^a dy′ · f(x′,y′) / (4πε₀ √((x−x′)² + (y−y′)²))
```

and the reader gives

```
L(f) = \int ^{dx'}\int \limits_{-a}^{a}\dfrac{f(x'}{},y')4\pi \varepsilon _{0}
       \sqrt{(x - x'})^{2}+(y - y')^{2} - aa\int
```

**What makes it hard is that everything crosses a level boundary.** The tree,
after the repair passes:

```
LINE n=10           L ( f ) =  INT₁  INT₂  FRAC  ","  "y"
CHAR ")"            ← still the numerator's, one level up
LINE n=14           ← the DENOMINATOR, one level up
CHAR ")"  SCRIPT ²
LINE[null] LINE(-a) LINE(a) SIZE SYM  ∫     ← INT₁'s limit block
SIZE SUB   LINE(-a) LINE(a) SIZE SYM  ∫     ← INT₂'s, claimed (§18)
```

Three separate things, none reachable by a rule that stays inside one list:

1. **INT₁'s limit block opens with a null LINE**, not a size switch, so
   `BigOpDisplayPass` never starts on it, and `BigOpDisplayAltPass` \u2014 which does
   read that shape \u2014 begins from an operator in its OWN list, and this operator
   is inside a line. Both halves are true at once. **Letting a block open on a
   null LINE was tried and reverted**: it claimed the limit blocks of 32 working
   summations (`\sum\limits_{n}` became `\sum n`).
2. **The numerator continues past the template** (`,y` inside the line, `)` in
   the parent) \u2014 and it is detectable, because the numerator holds an unclosed
   `(`. `DisplayFractionPass` already uses exactly that signal for the
   denominator; the symmetric rule for the numerator would have to reach up a
   level to find the closing character.
3. **The denominator is in the parent list**, past that `)`, so neither the
   separator/terminator path nor the immediate-sibling rule (§19) can see it.

The chunk splice in `PassPipeline::process` is the place where a cross-level
repair belongs \u2014 it already opens PILEs and chunk LINEs so that a template and
its continuation land in one list. `LINE n=10` is not opened because
`continuedBy` does not recognise what follows it as a continuation.

**Why this is where to stop rather than push on.** A wrongly-repaired fraction
reads as finished; the corpus score is marker-based and cannot tell a correct
reassembly from a plausible one; and the one rule tried in this direction cost
32 working documents. The next attempt needs the EE3 rendering of every document
it changes as its acceptance test, not the marker count.

`harrington_ch2_polarizability_xx_matrix.eqn` is clean by every marker and by
its fence, but its `\sum` still leaves the `n` at the end of the line rather
than under the sign \u2014 the same class as (1) above, and invisible to the score.

---

## 21. CLOSED — Ctrl+C did not paste into PowerPoint at all

Sugahara, 2026-08-21: *powerpointに貼り付けた場合に24ptはできている？*

The honest answer at the time was **no**, and finding that out took asking the
question the way the user asks it.

`radia.equation.office.copy_to_clipboard` → paste → **24.0 pt, a text shape**.
The editor's own Ctrl+C → paste → **PowerPoint refuses**: *"the clipboard is
empty or contains data that cannot be pasted here."* Nothing on the slide, no
picture, no fallback.

**One byte.** The editor's `put()` asked `GlobalAlloc` for `size + 1` and wrote
a NUL at the end. That is right for the clipboard's TEXT formats — RTF and
MathML are read as C strings — and fatal for the GVML package, which is an OPC
**ZIP**: PowerPoint will not open an archive with a byte after the end.

Measured directly, same package, two buffers:

| payload | paste |
|---|---|
| exact bytes | pasted, 24.0 pt |
| bytes + one NUL | **failed** |

`put()` takes a `terminate` flag now; GVML, PNG and CF_DIB pass `false`.

### Why the existing test could not see it

It builds the package with `tex_to_gvml` and puts **that** on the clipboard
itself. It proves the package is right. It says nothing about whether the
editor puts that package on the clipboard — and the editor did not.

`validation_test/equation/test_editor_clipboard_bytes.py` closes the gap by
driving the real gesture: launch the window, post Ctrl+C, read the bytes back,
compare with `tex_to_gvml` **byte for byte**. It needs no Office and takes three
seconds. On the broken build it says:

```
the editor wrote 2599 bytes where the package is 2598; a byte after the end of
an OPC archive makes PowerPoint refuse the paste outright
```

Driving the window from a test needs two things that are easy to get wrong and
are written down in that file: `SetKeyboardState` writes the CALLING thread's
input state, so the queues must be attached first (`AttachThreadInput`), and the
modifier has to stay down **across** the posted message — clearing it straight
away races and `Ctrl+C` arrives as a plain `c`.

### The general lesson

**A test that stands in for the caller does not test the caller.** This is the
third time in two days the same shape has appeared: a `tests/` fixture built
with our own writer passed against a broken reader (§18), a corpus score could
not see whose fault a defect was (§18) or half a fraction (§19), and here a
clipboard test supplied its own clipboard. When the claim is about what a user
does, something has to do what the user does.
