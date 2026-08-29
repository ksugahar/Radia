# Eqnedit64 changelog

## Unreleased

### Native/Web geometry and Office MathML parity

- The native Geometry tab now exposes the browser editor's same 11
  differential-geometry entries in the same order.  The native over/underline
  palette remains independently available under Basic.
- Both products publish inline 24 pt fallback MathML with the same large-operator
  semantics: sums use above/below limits and ordinary integrals use TeX-style
  side limits.  Native copy also supplies CF_HTML containing the exact same
  MathML as its registered `MathML` formats.  The browser HTML adds conditional
  OMML so PowerPoint preserves fractions, radicals, n-ary operators, overlines,
  and underlines as editable Office Math; its Office HTML import remains 18 pt
  because browsers cannot publish Windows' registered MathML format.
- Palette parity and real PowerPoint OOXML/render coverage now include the
  sum/integral and overline/underline constructs that exposed the mismatch.

### Always-visible math alphabets

- Added a fixed `R x` / `I x` / `B x` group to the native and Web palettes.
  It applies upright Roman (`\mathrm`), explicit math italic (`\mathit`), or
  vector bold (`\mathbf`) without changing subject tabs.
- Native canvas selections are restyled structurally and unselected input is
  persistent; native/Web TeX selections are wrapped and empty selections put
  the caret inside a new command. Explicit Roman and italic now survive
  parse/save/reopen instead of collapsing to `\text` or an implicit variable.
- Completed the extended alphabet family in both products: `\mathsf`,
  `\mathtt`, `\mathcal`, `\mathbb`, `\mathfrak`, `\bm` (`\boldsymbol` input
  alias), and the `\mathnormal` reset. Native rendering, TeX round trips, and
  Office MathML preserve these semantics.

### Session-safe embedded font

- Replaced `AddFontMemResourceEx` after reproducing two Server 2022
  `fontdrvhost.exe` access-violation crashes while the old hidden UI test still
  returned success. The standalone EXE now extracts its embedded Latin Modern
  Math bytes to a content-addressed per-user cache and uses file-backed,
  process-private registration without installing a font or writing registry
  entries.
- CI now runs 32 isolated application lifecycles and fails on a changed
  session font-host PID or a new Application Error event for `fontdrvhost.exe`.
  Category-label pixel checks also reject tiny-fragment output instead of
  accepting four changed pixels as visible text.

### Visible native toolbar text

- Palette category tabs and palette buttons retain native keyboard, focus,
  radio, and accessibility behavior while Eqnedit64 now draws their pixels
  itself. This avoids a Windows GDI/theme session failure where the controls'
  text and fonts were correct but every toolbar label appeared blank.
- The hidden interaction test renders each real button once with and once
  without its label and requires changed text pixels, closing the gap left by
  the old memory-DC-only font probe.

### Japanese text size

- `\text{...}` Japanese now uses an explicitly measured CJK face. It no longer
  passes through Latin Modern Math font linking, which could turn five 12 pt
  characters into a 306 pt-wide run on LAB.

### Radia-owned Web edition

- Imported the current laboratory-homepage JavaScript editor into `web/` and
  made Radia the source of truth for both native and browser editions.
- Added a checked HTML mount fragment and repository contract tests for palette
  grouping, Tab slot traversal, inline MathML Office copy without competing
  PNG, and the separate PNG action.
- Added `radia_mcp.presentation` policy/tool coverage so MCP workflows know the
  canonical sources, homepage publication role, TeX-only contract, and the
  split between batch OMML deck generation and interactive Eqnedit64 CLI use.

### Session-safe automated GUI tests

- The menu-response and real-window fuzz suites now run all items/seeds in one
  hidden process instead of launching hundreds of short-lived processes.  This
  retains real `WM_COMMAND`/window-procedure coverage while limiting embedded
  private-font registration to one lifecycle per suite.
- Both suites record progress for exact hang/seed replay and fail if the active
  Windows session's `fontdrvhost.exe` is replaced while they run.  They never
  take foreground focus or open Explorer during their background checks.
- The executable and Python module now require their embedded Latin Modern Math
  resource.  The old fallback that copied an OTF to `%TEMP%`, registered it by
  filename, and left the copy behind has been removed.

### Cross-version product contract

- Eqnedit64 and the laboratory Web/JS editor now share an explicit product
  priority: visible editable Office Math through MathML, dual GUI/palette and
  TeX-source learning, then the native editor's Microsoft Equation 3.0-derived
  shortcut compatibility and extensions.  Common changes must inspect and test
  both implementations; clipboard parity is judged by the Office result rather
  than identical platform formats.
- `Tab` / `Shift+Tab` in the native TeX source now move through the next /
  previous literal empty `{}` without inserting a tab character.  The same hand
  motion already works in the Web source and in the native structural canvas.

### Reliable PowerPoint copy

- `Ctrl+C` on the structural canvas now copies the complete equation when
  nothing is selected.  A structural selection still copies only that range;
  source-pane copy keeps the normal text-selection rule.
- The external PowerPoint gate no longer calls the clipboard publisher helper
  directly.  It sends the real no-selection GUI copy command, then requires
  editable Office Math with fraction and radical structure, so menu/command
  routing regressions cannot pass behind a working helper.
- Restored the registered `MathML` / `MathML Presentation` clipboard route used
  by the last version confirmed visible in real PowerPoint.  The newer inline
  CF_HTML route could create a native math zone whose text box nevertheless
  looked empty.  Visible native equations are again the compatibility gate;
  PowerPoint may report 18--24 pt depending on its conversion path.
- The external gate now exports the pasted shape with PowerPoint's own renderer
  and requires the fixed fraction test equation's visible silhouette, including
  a long fraction rule.  XML, `MathZones`, font metadata, or a tofu/replacement
  glyph can no longer let an empty-looking text box pass.  The test refuses to
  attach when a user's PowerPoint process is already running.

### Arbitrary-size matrices

- Matrices are no longer limited to 2x2 and 3x3 entry points.  The Basic
  palette now offers common row, column, rectangular, and 4x4--6x6 shapes;
  row/column add and remove actions grow any matrix up to 99x99 while keeping
  one row and one column.  `Ctrl+Alt+Arrow` performs the same structural edits.
- Ordinary matrices now support vertical caret movement.  Every resize is a
  named Undo/Redo step, preserves unaffected cells, and is exercised through
  both the headless model and real hidden `WM_COMMAND`/keyboard paths.
- Canonical TeX preserves trailing `&` cells and emits `{}` only for the final
  empty row of a one-column matrix, preventing 1xN and Nx1 matrices from
  shrinking on save/reopen without making ordinary empty cells verbose.  The
  deterministic 25,000-operation fuzzer found and now guards this bug.
- A repository-local running Eqnedit64 is forcibly stopped at build start so
  the signed portable executable can be replaced; the resolved-path guard does
  not terminate unrelated or separately deployed copies.

### Canonical Eqnedit64 operation

- Eqnedit64 and UTF-8 `.tex` are now the laboratory's canonical equation
  workflow. `build\accept_release.ps1 -Deploy` runs the normal suite, external
  paste APIs, ASan, 500,000 hidden real-window operations, release identity,
  and preserved Eqnedit32 checksum/signature gates before updating registered
  operational copies.
- A new offscreen visual-scale gate renders a representative research equation
  at 96, 120, 144, and 192 dpi equivalents and checks real pixels for equation
  ink, the structural caret, selection feedback, and palette/source fonts.
  The hidden interaction suite also proves that the shortcut coach returns
  focus and does not interrupt a 64-character typing burst.
- The final gate records the real offscreen GDI paint benchmark and rejects a
  cached representative-equation repaint at 5 ms or slower, keeping the
  keystroke path well inside a 60 Hz frame budget.
- Eqnedit32 remains an unchanged, Microsoft-signed comparison oracle under
  `reference`; it is no longer an operational editor and can neither be a
  deployment target nor a fallback document format. `.tex` is the asset of
  record; new `.eqn`/MTEF assets are outside the product contract.

### Live TeX source

- **The native editor now has the web editor's two live editing surfaces.**
  The structural canvas and a permanently visible TeX source pane update the
  same equation in both directions.  Complete, incomplete, multiline, and
  Japanese source edits are parsed on every change without replacing the
  half-typed text or moving its caret; leaving the pane, editing the canvas,
  Undo/Redo, and saving show canonical TeX again.
- A continuous source-edit burst takes one structural Undo checkpoint.  The
  source pane is permanently visible, and either pane is edited directly by
  clicking it.  There is no GUI/TeX mode or focus-switch shortcut; the TeX
  menu and its numbering toggle have now been removed altogether.  New files
  use `equation`; an opened `equation*` file retains its envelope on save.  TeX
  commands belong in the source pane; the canvas no longer has a separate
  `\command`-then-Space mode.
- The hidden interaction test covers complete and half-typed source, both
  synchronization directions, burst Undo, permanent layout, focus, absence of
  canvas command interpretation, and menu state.  Source editing is also mixed into the deterministic hidden GUI
  fuzzer, closing the crash-history gap that caused the earlier pane to be
  removed.
- The palette now teaches TeX without another prose layer: highlighting any
  cell shows canonical `TeX: ...` generated by executing that cell's real
  insertion command, and clicking selects only the newly inserted spelling in
  the permanent source pane.  Structural shortcuts provide the same feedback
  without taking focus or the caret from the canvas; entering the source pane
  collapses the transient selection before typing, so it cannot accidentally
  replace the command.  When available, the shortcut follows on the same short
  status line.  Basic now follows the web editor's structures -> brackets ->
  decoration order; arrows and sets lead the Sets/Symbols group.

### Drawn like TeX

- **The canvas now draws with Latin Modern Math** -- the OpenType Computer
  Modern, the typeface TeX itself sets with -- embedded in the executable and
  loaded per process, so the single portable file needs no font install and no
  TeX installation.  Geometry had already been matched to 0.03 em, but geometry
  is not appearance: the same layout in Times/Cambria still looked nothing like
  a pdfLaTeX run.  Variables are drawn as their Unicode math-italic code
  points, the way TeX sets them.  GUST Font License; see
  assets/THIRD_PARTY_NOTICES.md.
- **Big operators, fences, and the radical use the font's designed sizes**
  from the OpenType MATH table instead of the base glyph scaled up.  Scaling
  widens a sign exactly as much as it heightens it: the display integral was
  24% too wide for its height, which is what pushed its curls into correctly
  placed limits; a tall parenthesis was a thin rubbery 0.38 w/h against TeX's
  0.76.  Integral aspect is now within 0.3% of pdfLaTeX, and display operators
  are centred on the math axis (61.2% above the baseline vs TeX's 61.5%).
- **TeX's style chain is followed**: a fraction sets its parts one style down,
  so nested fractions come out 12, 8 and 6 pt exactly as pdfLaTeX sets them
  (they all stayed 12 pt before), and a non-display fraction uses the tighter
  num2/denom2 shifts of tex.web 704.
- **Integral limits sit where TeX puts them**: lower +0.554 em and upper
  +0.996 em from the operator's origin.  They had been measured from the
  glyph's widest ink, 0.65 em too far right.
- **Every horizontal rule is default_rule_thickness** (0.040 em).  All four
  rule sites floored the thickness at 0.6 pt -- 25% too thick at 12 pt -- and
  the radical additionally scaled its bar with the content's height, which TeX
  never does.
- **Accents sit on the letter**: a combining accent is drawn entirely above
  its own baseline, and placing it by that baseline left the whole rise as
  empty space -- \tilde{a} floated 1.1 base-heights up where TeX puts 0.30.
  Accents are also drawn at the base's size, not script size, and the hat is
  U+02C6 rather than the ASCII circumflex.
- **A superscript clears its base's italic correction** (tex.web 756), read
  from the font's MATH table: the 2 of f^2 sat on the f's hook, 0.097 em short
  for f and 0.202 em for V, while x and a looked fine -- which is why it read
  as a problem with one letter rather than a missing rule.  Subscripts
  deliberately do not get it.
- **A function name is spaced as an operator**: sin ωt was set as sinωt
  because \sin arrives as letters whose first atom is ordinary.
- The radical's clearance follows tex.web 737 (rule + x-height/4, surplus
  split) instead of a flat 0.344 em that was nearly twice TeX's.
- 19 geometry measurements against pdfLaTeX, all within 0.03 em, no known
  differences; tools/tex_sweep.py renders 44 constructs both ways and compares
  scale-free ink numbers, so a construct nobody thought to probe still gets
  looked at.

### Everything the editor emits now compiles

- **\nsubset does not exist in LaTeX or any package.**  It had been emitted
  since the program began, so every equation containing one failed to typeset.
  U+2284 is now spelled \not\subset and the parser folds the pair back into
  one character.
- \hbar was missing outright; \angle and \circ could be written but not
  read back.  All added, with \mho.  The 其の他 palette leads with working
  physics symbols instead of card suits.
- \oiint and \oiiint are no longer offered as templates: one needs esint,
  the other exists only in packages that replace the whole math font.  Pasted
  ones still display.
- tests/test_tex_compiles.py runs every symbol and template through pdflatex.
  The package list (amssymb, cancel -- 5 commands) is checked both ways, so a
  new dependency fails and a stale entry fails too.

### Clipboard

- **PowerPoint paste uses the visible, editable registered-MathML route.**  An
  attempted inline-MathML CF_HTML path could leave a native Office Math zone in
  a text box that looked empty in real PowerPoint.  Both registered MathML
  names are therefore restored.  PowerPoint can choose an 18--24 pt conversion
  and a math paragraph; visible editable content takes priority over forced
  left alignment and a forced 24 pt size.
- **The metafile carries glyph outlines, not text**, so pasting into
  PowerPoint keeps TeX's letterforms.  It used to name its font, which is
  loaded into Eqnedit64's process alone -- everything else substituted (13.8%
  of pixels differed).  Recorded as filled polygons via the font's cmap, with
  a corrected frame mapping; build/test_emf_outlines.ps1 reads the metafile's
  records and fails on any text record.
- The SVG export writes size-variant glyphs as outline paths, so it no longer
  depends on the reader having the font.

### Simpler insertion surface

- **The native editor and web editor now share five palette concepts**:
  Basic, Analysis, Sets/Symbols, Geometry, and Greek.  The native toolbar shows
  five category tabs on its first row and at most five matching popup palettes
  on its second instead of exposing all 18 palettes at once.  All 215 cells and
  every existing shortcut remain available.
- The two top-level Template and Symbol menus are one Insert menu with the
  same five groups.  `tests/test_palettes.py` proves that every one of the 18
  palettes belongs to exactly one nonempty category; the hidden Win32 test
  switches every real tab and checks the corresponding child-window styles.

### Signed portable release

- Every normal build now carries an Authenticode developer signature from
  `ksugahar`; signing or local trust validation
  failure stops the build before `dist/Eqnedit64.exe` is replaced.  The
  distribution remains one statically linked executable.  The setup script
  creates a non-exportable, self-signed CurrentUser code-signing key for lab
  development; it does not claim public-CA or SmartScreen reputation.

### Release-gate reliability

- The one-command gate now includes 24 seeded runs and 72,000 operations
  through the real hidden window procedures before exercising PowerPoint and
  IrfanView.  The same fuzzer accepts an explicit executable, and
  `test_asan.ps1` runs the hidden checks and full fuzz load against the ASan
  build.
- The ASan builder no longer reports failure after a successful link by
  invoking a second source-less compiler command.  PowerPoint verification now
  accepts equivalent namespace placement on the inline `a14:m` element,
  reports each structural predicate, and retains a failing PPTX as evidence.

### Editing

- **The structural canvas now owns an explicit IME lifecycle.**  IMM32 places
  composition/candidate UI at the structural caret, committed `WM_IME_CHAR`
  input passes through the normal Unicode edit path exactly once, cancellation
  changes nothing, and focus loss clears composition state.  Hidden messages
  cover start, Japanese commit, end, and cancel without activating a desktop
  IME; visual placement with a user's IME remains a real-session check.
- **Editing marks make invisible structure optional, not permanent ink.**
  Empty slots use dotted boxes, explicit TeX spaces use small blue-grey dots,
  and `aligned` ampersand positions use blue-grey guides.  View -> Editing
  Marks (`Ctrl+Shift+8`) toggles them without changing TeX; clipboard PNG/EMF,
  Office paste, and exports continue to render with marks disabled.
- **Unavailable edit commands explain themselves.**  Cut and Copy require a
  selection on the focused canvas/source surface; Paste requires a supported
  TeX or text clipboard format; empty Select All and unavailable Undo/Redo are
  disabled.  Hovering a disabled item or invoking its accelerator shows the
  short reason in the status bar without changing the clipboard.  The
  whole-equation Google Slides export remains a separate command.
- **Undo and Redo say what they will do.**  History entries now retain names
  such as `Typing`, `Template`, `Style Change`, and `TeX Edit`; the Edit menu
  shows that name, disables an unavailable direction, and the operation/flight
  log records the same name when the action runs.
- **Ctrl+click selects one safe structural atom.**  Clicking inside a
  fraction, radical, fence, script, or matrix selects the innermost containing
  template together with all its slots.  Cut, delete, typing, and template
  replacement therefore cannot tear a visual template away from its TeX tree.
  Model and hidden real-Canvas tests cover selecting and deleting a fraction.
- **Drag selection follows long equations beyond the canvas edge.**  A
  captured 50 ms timer continues panning and hit-testing while the pointer is
  held still at an edge, clamps the equation to its content range, and releases
  the timer and capture together.  The hidden Win32 test selects through a
  128-character equation in a narrow canvas without moving the real cursor.
- **A deletion leaves the caret where it happened**: Backspace on x_1^2
  removed the 2 but left the caret outside, so typing gave x_{1}^{}9 instead
  of x_{1}^{9}.  Fixed for Backspace and Delete together -- the last two
  defects here were each fixed for one key and left standing for the other --
  and repeated presses now peel the equation to empty instead of stopping
  partway.
- Equations start at the top of the canvas rather than its vertical centre.
- The operation log's start/stop popups are gone; the status bar and the
  title flag say everything they said.

### Finding what a person would otherwise have to find

- **A flight recorder for freezes**: every interaction goes into an always-on
  ring buffer, every command is bracketed, and a watchdog dumps the ring to
  %LOCALAPPDATA%\Eqnedit64\hang-<time>.log when a handler runs 8 seconds
  without returning -- while the program is still stuck.  The last line names
  the freezing press.  Proven with a planted Sleep(60000).
- **Menu coverage in three layers**: --menu-audit dispatches all menu items
  and reports any with no effect (now 0, with every skip named and reasoned);
  test_menu_responds.ps1 presses every item in one hidden process with a
  suite deadline and per-item progress logging;
  test_menu_navigation.ps1 walks the real menu loop, and skips loudly when it
  cannot take the foreground instead of reporting phantom freezes.
- **Fonts are chosen by drawing, not by name.**  This machine reached a state
  where CreateFontW("Cambria Math") returned a font GDI rasterised to nothing
  -- every toolbar button blank on a build that had displayed them the evening
  before.  The math font retries until a probe glyph actually measures; the
  button font takes the first candidate that provably puts ink on a bitmap,
  and logs every candidate that failed.
- --debug-colors paints regions in flat unmistakable colours (canvas pure
  green, equation box pure yellow) so a screenshot can be judged by pixel
  count.
- The build stamp is in the title bar, About shows the running executable's
  path, and --version answers both without opening the editor: a day-old
  build left in a scratch directory was mistaken for a regression, and two
  copies of this program looked identical from across the desk.
- build/deploy.ps1 refreshes every copy listed in deploy_copies.txt from
  dist, refuses to ship a dist that differs from the tested build, and names
  a destination held open by a running instance instead of reporting a copy
  that did not happen.

### Kept from the Eqnedt32 review

(The counts below are as first shipped; the sections above have since grown
them to 165 symbols and 222 cells, and removed the two esint templates.)
Reviewed Eqnedt32's own implementation -- its menus, accelerators, dialogs,
string table and toolbar art, read statically out of the binary -- and adopted
the one design decision that mattered most.

- **The toolbar is now the catalogue, not a shortcut duplicate.** Eqnedt32 has
  no Insert menu at all: nineteen wide toolbar buttons, symbols on the top row
  and templates on the bottom, each dropping a small grid, and every one of its
  ~300 items is two clicks away.  Eqnedit64 had flat menus of 27 templates and
  24 hand-picked symbols, and the toolbar had been trimmed to three buttons on
  the reasoning that a command with a chord does not need a button.  Measured,
  that left **117 of 160 symbols reachable only by knowing and typing the
  command** -- neither a shortcut nor a mouse route.  There are now 18 palettes
  (10 symbol, 8 template) covering all 160 symbols and all 53 templates, and
  the Template and Symbol menus present the same catalogue as submenus so the
  keyboard route is complete too.  `tests/test_palettes.py` holds the rule:
  every symbol and every template must have exactly one palette home, and
  `--ui-interaction-test` inserts all 220 cells through the real WM_COMMAND
  path.
- Added 24 template kinds whose nodes the engine already rendered and emitted
  but which no menu, button, or shortcut could reach: `norm`, `floor`, `ceil`,
  `dirac`, `iiint`, `oiint`, `oiiint`, `coprod`, `bigcup`, `bigcap`,
  `overbrace`, `underbrace`, `overrightarrow`, `overleftarrow`,
  `overleftrightarrow`, `over`, `ddot`, `dddot`, `prime`, `dprime`, `tprime`,
  `strike`, `frown`, `smile`.

Five round-trip defects that the palettes would otherwise have exposed to a
user's mouse, all found by round-tripping every palette command before shipping
the buttons:

- `\overbrace` and `\underbrace` were written by the emitter but never read by
  the parser, so a horizontal brace saved and reopened came back as the literal
  text "overbrace".  Same for `\cancel` (the strike-through embellishment).
- The prime family emitted nothing at all: the embellishment writer required
  both a prefix and a suffix, and `x'` has only a suffix, so the mark was
  silently dropped and a prime template produced a bare `x`.
- `\|x\|` outside `\left`/`\right` read back as `|x|`, quietly turning a norm
  into an absolute value.  `\Vert` (U+2016) is now in the symbol table and
  `\|` canonicalises to it.
- `\bar{x}` was emitted as `\overline{x}`, collapsing TeX's distinction between
  a one-character accent and a stretchy decoration -- so `bar` and `overline`
  became the same thing on the next save.  `Ctrl+-` is `\bar`; `Ctrl+T` `-` is
  `\overline`.
- **Explicit spacing was discarded on load**: `a \quad b` came back as `ab`, so
  spacing an author had deliberately written was lost the first time the file
  was opened and saved.  `\!`, `\,`, `\:`, `\;`, `\ `, `\quad` and `\qquad` are
  now characters with TeX's own widths (3mu, 4mu, 5mu, 18mu), round-trip, and
  render as `<mspace>` in MathML.

Not yet matched from Eqnedt32, and why: left-hand pre-scripts and labelled
arrows (`\xrightarrow`) have no node type; long division has no renderer; the
Spacing dialog's 19 named parameters and the Style/Size definition dialogs are
a settings model this build does not have.

## 3.0.0 — 2026-08-22

- A new document starts empty.  The sample equation meant the first thing
  anyone did was delete it.
- Fixed Backspace destroying a whole template in one press: `E = mc^{2}` with
  the caret after it lost the `c` together with the `2`, and `\sqrt{x^{2}}`
  lost its `x` and its `2`.  One press now removes one thing the reader can
  see, descending to the innermost slot that still holds something; the
  template itself goes only once every slot is empty.
  This is stated as a rule in `docs/GUI_SPEC.md` and checked as a rule:
  `tests/test_edit.py` counts glyphs in the rendered SVG across scripts,
  fractions, fences, radicals, big operators and `cases`, so no template can
  fall outside it.  The earlier fix here addressed the one reported state and
  its regression was written to that same state, which is why the neighbouring
  state stayed broken and undetected; recorded as UXP-0010.
- Derived the layout constants from TeX instead of by eye, and checked them
  against pdfLaTeX.  The fraction rule was 0.39 em wider than TeX's and the
  numerator 0.24 em too close to it; scripts were set at 7/12 and 5/12 where
  TeX uses 8/12 and 6/12; superscripts rode 0.04 em high and subscripts sat
  0.07 em low.  Fractions now follow tex.web 704 (`num1`, `denom1`,
  `axis_height`, a 3-rule-thickness clearance, and a rule exactly as wide as
  the fraction).
- Added `tools/tex_geometry.py` and `tests/test_tex_geometry.py`: the same
  equation is typeset by pdflatex and by the canvas, glyph origins and rules
  are compared in em, and `tests/tex_reference.json` carries the pdfLaTeX
  side so the check runs without a TeX installation.  Remaining differences
  are listed with their cause and the test fails if one silently starts
  matching.  Nothing before this compared the canvas with anything but
  itself.
- Joined the radical's vinculum to the glyph rather than to its advance
  width.  Cambria Math's flag is a 0.065 em plateau ending 0.06 em past the
  advance, so a bar drawn at the nominal rule thickness from the advance came
  out thinner than the flag and started past its end.  Thickness and start
  are now measured from the glyph.
- Made the canvas repaint the way Eqnedit32 does.  The draw loop created and
  destroyed an HFONT for every glyph on every keystroke, and `WM_PAINT`
  painted straight onto the window after filling the whole client area.
  Fonts are cached, `WM_PAINT` composes offscreen and blits, and the thirteen
  whole-window invalidations no longer ask for an erase.  `--paint-bench`
  times both paths in one binary: 1.4x to 1.8x.

- Fixed `\rangle` serializing without its trailing space, which turned
  `\langle a \rangle x` into `\langle a\ranglex` and brought it back from the
  file as the word "ranglex".  Undo restores a snapshot by re-parsing its own
  output, so the same hole corrupted Undo as well as save and reopen.
- Fixed 45 symbol commands being written to `.tex` as bare Unicode glyphs
  (`⟹`, `⌊`, `†`, `…`) that pdfLaTeX cannot set in math mode.  A character
  parsed from a named command now keeps that command and is written back
  unchanged, so `\epsilon`, `\ldots`, and `\to` no longer come back as
  `\varepsilon`, `\cdots`, and `\rightarrow`.
- Fixed `\parallel` degrading to `\|` and then to `|`, and `\ast` to `*`, on
  one save-and-reopen cycle.  Separated `\cdots` (U+22EF) from `\ldots`
  (U+2026) and moved `\langle`/`\rangle` from the full-width CJK brackets
  U+2329/U+232A to the mathematical U+27E8/U+27E9 the canvas already drew.
- Fixed `\div`, `\Re`, and `\Im` being read as operator words and set upright
  as "div", "Re", and "Im"; named glyphs now win over operator names.
- Fixed `cases` losing its second column: `x & x > 0` was flattened to
  `x x > 0` and the environment was rewritten as `gathered`.  The `cases`
  template and the parser now share one two-column shape that round-trips as
  `cases`.
- Fixed a mismatched closing delimiter being replaced by a mirror of the
  opening one, so the half-open interval `\left( a,b \right]` survives, and
  fixed `\left. x \right.` growing the parentheses of the tmPAREN fallback.
- Fixed opening a file whose read failed loading an empty equation, keeping
  the path, and destroying the original on the next save.  `read_file` now
  reports failure separately from an empty file.
- Made outer whitespace normalization idempotent so a first serialization
  matches the next parse.
- Added `tests/test_symbols.py`, which sweeps every command in the symbol
  table for round-trip identity, a one-pass fixed point, and ASCII TeX
  output, plus `cases` columns and mixed and invisible fences.  The raw-TeX
  fuzzer's corpus named 12 of the 160 commands, which is why the defects
  above shipped with a green suite.
- Added `build\test_tex_document.exe` to `build_tests.bat` and to the
  required checks; it was built by every release and never run.
- Made saving atomic: the equation is written to a sibling temporary file and
  renamed over the target, so a write that fails part-way no longer leaves a
  truncated `.tex` where the saved equation used to be.
- Fixed editing the raw TeX pane clearing the undo stack.  It called
  `load_latex()`, a document load, so one keystroke there discarded the whole
  canvas history.  A pane edit now takes one checkpoint per burst, and Ctrl+Z
  returns to the equation as it stood before the pane was touched.
- Fixed `\lim` and other operator names being wrapped in braces when they
  carried a subscript.  `{\lim}_{x \to 0}` is an ordinary atom, so in display
  style the limits sat beside the operator instead of under it.
- Fixed `\middle` printing the letters of its own command name: the bra-ket
  `\left\langle a \middle| b \right\rangle` came back with the word "middle"
  in it.  The delimiter itself is kept instead.
- Added `LegalCopyright` to the version resource, which was empty.
- Turned on `/WX` for the executable, the C++ test, and the headless module,
  so "builds without warnings" is enforced rather than asserted.
- Extended `--ui-interaction-test` to cover the chords that live only in
  `canvas_keydown` -- zoom, alignment, input style, and the `Ctrl+T`/`Ctrl+K`/
  `Ctrl+G`/`Ctrl+B` two-stroke families -- through the real key path, using a
  per-process test modifier flag rather than desktop input.
- Marked `CMakeLists.txt` as an IDE convenience build that cannot produce a
  release, so it is not mistaken for the shipping path.
- Fixed six variant Greek commands sharing a code point with their plain
  form.  LaTeX's `\phi` is the straight symbol and `\varphi` the ordinary
  letter, and the pair was inverted; `\epsilon`, `\vartheta`, `\varpi`,
  `\varkappa`, and `\varrho` had no distinct glyph at all, so the canvas and
  the MathML sent to Office both showed the wrong letter.  Saved TeX is
  unchanged, because a parsed command is written back verbatim.
- Made the variant Greek letters italic like every other lowercase Greek
  letter; they sit outside the Greek block, so a range test had left them
  upright.  `typeface_for_code` now has one definition instead of a copy in
  the parser and another in the editor.
- Turned `\left\langle a \middle| b \right\rangle` into the bra-ket node the
  tree already had, so Tab moves between bra and ket and the bar and angle
  brackets stretch to the taller side.  Other uses of `\middle` keep their
  delimiter inline.
- Documented that `-` and `=` are one physical key on a Japanese keyboard, so
  `Ctrl+Shift+=` is zoom-out there and the math input style is on the `+`
  key.  Zoom-out has no other chord; the input style does.
- Extended `--ui-interaction-test` to the layout-resolved chords as well:
  math style, zoom out, overline, and the `Ctrl+T` members that need Shift.
- Removed a dead embellishment branch in the function-name emitter, gave
  Ctrl+A its status update, and stopped `wWinMain` reporting an uninitialised
  exit code when the message queue errors.
- Fixed a root index containing `]` being written unbraced.  TeX ends an
  optional argument at the first unbraced `]`, so `\sqrt[{]}]{x}` came back
  from the file as `\sqrt{]}x`: index gone, wrong radicand, and the operand
  outside the root.
- Fixed a script base that serialises to several tokens being braced
  inconsistently.  A text or vector character came out as `\text{a}^{2}` when
  typed and `{\text{a}}^{2}` when read back from a file, so the document
  changed every other time it was opened.  Bracing is now decided from the
  emitted text rather than from the node shape.
- Gave `cases` its own left-aligned layout.  It shared the centred matrix
  layout, so the canvas disagreed with the PDF it was about to produce, and a
  literal `\left\{ \begin{matrix} ... \right.` was silently rewritten as
  `cases`.
- Kept the save temporary file short and distinctive.  The first spelling
  pushed paths within a dozen characters of MAX_PATH over the limit, turning
  a working save into a refusal, and a bare `~` would have overwritten the
  editor backup files some tools keep beside the document.  When no sibling
  temporary is possible at all, the write falls back to writing in place
  rather than refusing to save.
- Widened the raw-TeX fuzzer's corpus to the shapes with special parse rules
  -- bra-kets, mismatched and invisible fences, a braced matrix -- and to
  named glyphs, an operator name, and a variant Greek letter.  The three
  defects above were all found by that corpus within one run.

- Fixed transparent raster paste by forcing every 32-bit DIBV5 alpha byte to
  255 after GDI rendering; regression tests now require opaque white pixels,
  dark ink pixels, and zero transparent pixels in normal copy and `--texclip`.
- Added 24 pt Presentation MathML clipboard formats so PowerPoint's normal
  paste creates an editable native Office Math object at 24 pt instead of its
  18 pt Unicode-LaTeX default; all 3,072 deterministic TeX fuzz cases now also
  require well-formed 24 pt MathML.
- Fixed Backspace and Delete on an empty outer sub/superscript so only the
  empty script shell is removed and its visible base and other scripts remain.
- Added an exact Win32 regression for Backspace from a populated superscript:
  `H^{2}|` becomes `H^{}|`, then `H|`, without desktop input automation.
- Kept an automatically recognised function upright when its argument is typed
  immediately after it: character-by-character `sinx` now saves as `\sin x`,
  while a longer recognised name such as `sinh` still wins.
- Moved integral upper/lower limits to the measured TeX display-integral
  positions. At 12 pt the Eqnedit64 offsets are -13.200/+10.800 pt with a
  +5.280 pt upper italic correction versus TeX -13.117/+10.793/+5.313 pt.

- Rebuilt the editor as a native 64-bit Windows application.
- Made UTF-8 TeX the source of truth for opening, editing, copying, pasting,
  and saving equations.
- Added structural caret movement, drag selection, selection replacement,
  undo/redo, and familiar template shortcuts.
- Added multiline equations with `aligned`, Enter, alignment tabs, and
  vertical row movement.
- Made the editing canvas left-aligned by default so the equation does not
  shift while typing; saved `equation` TeX remains normally centred.
- Replaced the space-padded static status text with a native partitioned
  Windows status bar, including hidden layout/font/simple-mode regression tests.
- Prevented Cut from deleting a selection when clipboard publication fails.
- Added structural paste for Unicode text and common TeX math wrappers.
- Added multi-format clipboard copy: Office-recognisable delimited LaTeX for
  native PowerPoint/Word/Excel equations, scalable EMF, 32-bit DIBV5, and a
  raw registered LaTeX payload from the same copy command.
- Replaced the installer release with one statically linked portable
  `Eqnedit64.exe`; it performs no registry or Start-menu registration.
- Added GDI canvas rendering, SVG export, hit testing, and scale-aware layout.
- Restored Eqnedit32's structural double-click rule: select the current input
  slot, and select the whole equation only from the outermost slot.
- Added a licensed application icon and Windows version resources.
- Added background-only interaction tests that do not control the desktop.
- Added API-driven external paste regression: hidden PowerPoint COM verifies
  editable Office Math XML, IrfanView CLI verifies a nonblank pasted image,
  and the test snapshots and restores the user's clipboard formats.
- Added a separate Google Slides copy command (`Ctrl+Alt+C`) that publishes a
  300 dpi PNG rendered from a 24 pt base style plus point-sized clipboard HTML,
  without changing PowerPoint's normal native-equation paste selection.
- Added a hidden `--texclip` command that replaces clipboard TeX with a
  standalone 300 dpi PNG rendered from the same 24 pt base style, preferring
  the registered raw `LaTeX` format and falling back to Unicode text; the same
  raster is also published as DIBV5 for Windows image-editor paste support.
- Added deterministic raw-TeX fuzzing (3,072 generated and damaged inputs)
  covering normalization fixed points, finite layout metrics, and valid SVG.
- Fixed UTF-8 byte splitting, supplementary-plane character loss, nested-script
  reassociation, malformed escaped-literal instability, and empty-root-index
  two-pass normalization found by that fuzzer.
- Added hidden Win32 message-path tests for character/surrogate input, caret
  keys, undo, templates, source editing, view commands, and status updates.
- Preserved text, function, and vector input styles inside structural slots
  such as fraction numerators, with save/reopen regression coverage.
- Added automatic Function styling in Math input for legacy names such as
  `sin`, `log`, and `exp`.
- Restored the documented prime and double-prime chords to the shortcut guide.
- Added opt-in operation logging, F12 problem markers, a guided manual test,
  and a dedicated operation-debug shortcut.
- Upgraded operation logs to semantic v2 traces with monotonic event timing,
  focus, style, alignment, zoom, equation mode, and shortcut
  prefix state while retaining v1 analysis compatibility.
- Added a privacy-aware LLM usability-review pipeline that detects evidence
  windows around explicit markers, invalid shortcuts, no-ops, immediate undo,
  navigation reversals, and correction bursts; human decisions are retained in
  a machine-readable preference ledger and require background regressions.
- Restored Eqnedit32's documented Ctrl+T, Ctrl+K, Ctrl+G, and Ctrl+B key
  families, then added non-conflicting Eqnedit64 shortcuts.
- Added an optional shortcut coach that teaches the chord after a menu or
  toolbar action.
- Removed MTEF/.eqn compatibility readers, writers, datasets, and all installer
  paths.  The signed Eqnedit32 original and help remain as a staged-retirement
  compatibility reference until Eqnedit64 is accepted.
