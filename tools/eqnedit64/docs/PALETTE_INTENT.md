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
selected states. `test_palette_visual` creates 114 PNGs and a command index;
CI retains them as `eqnedit64-palette-proofs`. The images use the production
34-by-28 logical pixel cells and font chooser. Added gray borders are proof
guides, not a screenshot of Windows menu chrome. Run the proof executable
only in session 0 or the explicitly isolated CI environment.

The cell renderer crops measured ink, preserves aspect ratio and a clear
margin, and never enlarges a glyph beyond its nominal size. The gate checks
visible ink AND an empty outer pixel border in both states, not just a font
name or any nonzero pixel. AI additionally checks recognizability and semantic
distinctions; automation alone cannot certify those.

Fraction/radical, accent/line and typeface commands show native-model previews
instead of ambiguous literal abbreviations. Fractions use distinct `a`/`b`
slots, decorated variables use `x`, typefaces use `A`; style selection must be
cleared before rendering. The preview is derived from the insertion command,
not from a separate hand-written TeX implementation. Literal catalogue faces
remain the fallback/accessibility labels and the Windows glyph coverage input;
the intent oracle still tests the actual insertion independently. Web rendering
is unchanged in this native fix; equivalent visual goals apply to Web review,
but GDI cropping must not be copied into the browser renderer.

A human may still resolve ambiguous legibility or monitor-specific perception.
AI images are not evidence of popup positioning, pointer behavior, or the
health of LAB's interactive font session. Record untested dimensions explicitly.
