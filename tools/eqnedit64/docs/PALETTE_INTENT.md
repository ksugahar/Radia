# Palette intent acceptance contract

`palette_intent.json` is the reviewed, independent expectation for a key's
visible mathematical meaning and emitted TeX. It is test input, never production
input and never regenerated from the implementation or recorded output. New
keys require a reviewed expectation in the same change. Unknown keys fail.

## Reading the specification

- `symbols`: `visible character:TeX name[,permitted alias]`. For example a left
  arrow promises `\leftarrow`, not just an arbitrary stable arrow. `not` is an
  explicitly word-labelled combining-overlay command.
- `native_templates`: command identifier -> `[face, slot input, expected TeX]`.
  Type the listed distinct characters into successive empty slots using Tab.
  `frac` promises `\frac{a}{b}`; `nthroot` promises `\sqrt[b]{a}`; over/under
  annotations keep the base `a` and annotation `b`. Different operands make
  swapped numerator/denominator, limits, scripts and matrix cells observable.
- `native_matrix_sizes`: all offered dimensions. Fill row-major with distinct
  characters from `abcdefghijklmnopqrstuvwxyz0123456789`.
- `native_styles`: select `a`, apply the style, compare the promised wrapper.
  The persistent roman/italic/vector buttons additionally promise `\mathrm`,
  `\mathit`, and `\mathbf`; check their UI dispatch and actual model output.
- `native_raw`: face -> literal insertion. These are complete native fragments;
  a script without a preceding atom needs an explicit empty base.
- `matrix_actions`: start with rows `[a,b]`, `[c,d]`, caret in the first cell;
  verify all surviving cells and the placement of the inserted/deleted row or
  column, not merely the new dimensions.
- `web_overrides`: category/face -> insertion snippet. Other Web keys use the
  shared symbol expectations. Check both caret insertion and selection
  replacement/wrapping through the actual Web insertion function. The three
  persistent typeface buttons are checked separately.

## Independent evidence, not circular approval

The static gate compares every native command's face with this specification.
The native runtime gate executes the catalogue's actual command, fills its
slots, and compares emitted TeX with the independently specified expression.
The parser canonicalizes the expected spelling (aliases, spacing, grouping);
it does not choose the expected operator or operands. This gate is not an
independent proof that the parser itself is correct. Parser round trips and
rendering checks remain separate obligations.

In particular, `E.tex_normalize(expected)` is a shared implementation
dependency: a normalizer defect can move both sides of the comparison together.
The mutation controls cover selected swaps, not every such correlated failure.
A green run alone therefore does not establish oracle correctness; review the
authored expectations against the mathematical meaning independently.

Mutation controls deliberately interchange left/right arrows, over/under
lines, sub/superscripts and typefaces. Each must fail despite producing valid,
distinct, round-tripping mathematics. The existing distinct-rendering sweep
still detects renderer fall-through; the font-cmap gate still detects missing
button glyphs. Neither replaces semantic expectations.

Run safe static and Web checks in the Ubuntu Eqnedit64 lane. Run native insertion
in the existing single-process model suite on disposable Windows CI, never on
the interactive LAB font session. These tests do not claim that Win32 pointer
hit-testing or PowerPoint paste has been exercised.

Font cmap coverage requires every Windows Unicode subtable (platform 3,
encoding 1 or 10), never the union of Mac Roman and Unicode tables. Tab faces
are included: use U+2044 with the radical and `x` with U+2191, not the Mac-only
U+00BD or U+00B2 entries of the shipped font. This predicts which font
`pick_button_font` will accept; it does
not observe which font a running process received. `--self-test` now makes that
observation: it calls the production chooser, reads the physical face GDI
resolved at 96, 120, 144, 192, 288, and 384 dpi, and exits 243 naming the
substitute and DPI if the palette is not drawing in Latin Modern Math.
Diagnostics separate physical-face substitution, glyph coverage, and missing
ink; they do not infer a root cause from the selected fallback alone. The ink
probe allocates its surface from the selected font's measured line height and
text extent, so a large Windows ascent cannot clip the probe at high DPI.

The candidate before this correction passed English Windows CI but was reported
to return 243 on Japanese Windows on both 100 and mdx1. The Windows cmap gap
and fixed 64px ink surface are established defects; code page 1252 versus 932
as the reason for the differing machines remains a hypothesis. CI success is
not proof of success on LAB. Require a Japanese Windows isolated-session check
of the signed candidate before publication; never run this check in LAB's
interactive session.

Do not route this through `flight_note`. Those notes live in a ring buffer that
reaches a file only on a crash or an eight-second watchdog freeze, so a healthy
run records nothing and an instruction to read `font.buttons` cannot be carried
out.

## Visual acceptance, including AI review

Do not delegate basic clipping discovery to the user. AI reviews the actual
owner-draw proofs first: all 19 palettes, at 96/144/192 DPI, in normal and
selected states. `test_palette_visual` creates 120 PNGs and a command index;
CI retains them as `eqnedit64-palette-proofs`. The images use the production
34-by-28 logical pixel cells and font chooser. Added gray borders are proof
guides, not a screenshot of Windows menu chrome. Run the proof executable
only in session 0 or the explicitly isolated CI environment. The additional
six sheets exercise the 19 palette selectors and three persistent style faces
through the same ink-fitting helper used by the toolbar, including the actual
italic/bold font variants. Selector proofs use 52-by-30 logical pixels; the
58-pixel persistent buttons are conservatively tested at that narrower width.

The cell renderer crops measured ink, preserves aspect ratio and a clear
margin, and never enlarges a glyph beyond its nominal size. The gate checks
visible ink AND an empty outer pixel border in both states, not just a font
name or any nonzero pixel. AI additionally checks recognizability and semantic
distinctions; automation alone cannot certify those.

Thin/medium/wide spacing commands use one/two/three explicit spacing marks.
The font's U+2423 hairline disappeared during downsampling at 144 DPI in the
first implementation; the new visible-ink gate caught it. Spacing marks must
remain visible in normal and highlighted states without relying on that glyph.

Both editions take their grouping, order, labels, insertion adapter and preview
TeX from `shared/palettes.json`. Run `build/generate_palettes.py` after editing
it; CI's `--check` rejects stale C++ or embedded JS output. The generator never
reads the independent `docs/palette_intent.json` oracle. Every entry must declare
`preview_tex`; an explicitly empty string means a literal face (matrix dimensions
and contextual controls, for example). A failed nonempty native preview must
fail the visible-ink gate, not silently switch to a literal label.

The common catalogue has 245 keys in 19 palettes and five categories. Existing
Web-only snippets remain in an explicit Web-additions tab; four contextual
matrix row/column operations are disabled with an explanation on Web, where
editing remains source-based. This is common catalogue ownership, not a claim
that both editing engines or platform capabilities are identical. No new
Python-facing API is required. The browser uses the declared body slot when
wrapping selected text (root and overset indices precede their bodies in TeX).

Web uses MathJax for the shared examples and retains literal accessible names.
Preview failures are visible, not raw TeX. Native uses its existing model and
GDI fitter. Over/underbrace examples use a two-letter base and omit annotation
ink so the brace can remain larger; native previews disable editing placeholders.
Inserted expressions still include those editable annotation slots.
The image audit exposed a renderer fallback that changed wide braces to rules.
Braces must retain their glyph at every width; horizontal scaling is shared by
GDI and SVG, and the decoration regression checks one-, two- and eight-letter
bases in both directions. Passing an ink/clipping gate alone missed this defect.

Pixel gates detect clipping and empty output, not minimum semantic feature
size. Small primes and harpoon distinctions still require image review. No
hover-performance improvement or cache claim is made without measurements;
the current native raster path is not cached.

A human may still resolve ambiguous legibility or monitor-specific perception.
AI images are not evidence of popup positioning, pointer behavior, or the
health of LAB's interactive font session. Record untested dimensions explicitly.

## Prime attachment

A prime decoration applies to its entire selected/template body. Native output
must protect a compound body (including existing superscripts) as
`{a^{2}}^{\prime }`, not `a^{2}'`. Apply the same rule for two and three primes.
The Unicode prime symbol itself emits `\prime`, not an apostrophe: the latter
creates a superscript and cannot substitute for a glyph inside a superscript.
Require the first save/reparse to preserve the emitted TeX. Web selection
decoration uses the same grouped-base form; unselected Web caret insertion
retains its existing empty-atom protection after an exponent. No automatic
rewrite of arbitrary manually typed TeX is implied.
