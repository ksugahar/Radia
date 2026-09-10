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

Font cmap coverage predicts which font `pick_button_font` will accept; it does
not observe which font a running process received. `--self-test` now makes that
observation: it calls the production chooser, reads the physical face GDI
resolved, and exits 243 naming the substitute if the palette is not drawing in
Latin Modern Math. The existing hidden-executable CI step therefore fails the
moment a face character leaves the embedded cmap again.

Do not route this through `flight_note`. Those notes live in a ring buffer that
reaches a file only on a crash or an eight-second watchdog freeze, so a healthy
run records nothing and an instruction to read `font.buttons` cannot be carried
out.

One judgement is still a person's. The automated gate proves the typeface;
whether each key is legible at its drawn size — the primes distinguishable from
quotation marks, the two-letter style faces readable — is decided by looking at
the staged EXE before publication.
