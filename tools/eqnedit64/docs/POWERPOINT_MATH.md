# Native PowerPoint math export

Version 3.1.5 is a review candidate. Version 3.1.4 is the published release.

```text
Eqnedit64.exe equation.tex equation.pptxmath
```

The output is an XML `m:oMathPara` in the 2006 Office math namespace, with
DrawingML `a:rPr` formatting and canonical Office math `m:t` text runs. It is intended for PowerPoint's `a14:m`
extension, not Word import and not HTML. The serializer shares equation
structure with the normal-copy RichEdit HTML route. Both routes emit matched
visible delimiter pairs as stretchable structures.
ISO `baseJc` is vertical and uses `center` here. Horizontal paragraph
justification remains `left`.

The strict export rejects unsupported commands/environments, invalid colours,
unclosed groups/environments and excessive depth. The known non-covered
boxed, binom and cancel constructs are rejected rather than flattened.
Other unsupported inputs must not be assumed valid simply because the editor
can show them. Function argument collection retains the shared one-atom or
balanced-visible-fence contract; arbitrary TeX parity is not established.

The CLI writes atomically, emits a non-zero exit and a diagnostic including
the input TeX on rejection, and preserves an existing output file. The
`eqnedit64.export_powerpoint_math(tex, output)` helper is a checked subprocess
client with no renderer, clipboard or image fallback. This path does not load
the private rendering font. Office math names Cambria Math, with 18 pt native
runs; no font is embedded into a PPTX.

## Markdown decks

`radia_mcp.presentation.presentation_markdown_to_pptx` requires
`eqnedit64>=3.1.4` and `radia-mcp[document,equation]`. It supports plain slide
headings, paragraphs, flat bullets/numbered lists, fenced/inline code, `>`
speaker notes, `$TeX$` inline and `$$TeX$$` display equations in body and notes.
Unsupported Markdown raises. Escape literal currency dollars as `\$`.

The adapter only assembles documents; native Eqnedit64 owns all mathematical
serialization. Math-bearing shapes use `mc:AlternateContent` with an `a14`
choice and no fallback representation. Modern PowerPoint is required.
`python-pptx` preserves this XML on reopening/saving, but its ordinary shapes
iterator does not expose shapes nested in the choice. XML-aware consumers
must inspect that choice; they must not silently treat those shapes as absent.
Prose follows `font_pt`; native math stays at 18 pt.
Display equations are left-aligned with 8 pt before/after paragraph spacing.
Matched visible pairs in one row, including plain parentheses, explicit
`\left...\right`, use a stretchable Office delimiter structure. The `\big`, `\Big`, `\bigg`, `\Bigg` commands (plain and `l/r/m` variants)
are fixed-size tokens in both clipboard and strict OOXML output. At the 18 pt
Office base, the four classes use upright symbol runs of 21.6, 32.4, 43.2 and
54.0 pt. They never grow with the operand or become a stretchable `m:d`.
Each fixed symbol is scoped inside its own `m:box` with an 18 pt control
property, so its enlarged font does not change an adjacent fraction's layout.
Native TeX round trips preserve the exact command. MathJax's font-specific
big-family heights are mapped to these shared Office classes; internal
binomial layout is not classified as an explicit size command.
Strict export rejects missing/non-delimiter arguments. Interactive native
parsing ignores an unfinished size command and retains an ordinary token,
so a half-typed command does not make the whole equation invalid.
General `\middle` separators and one-sided fences currently retain visible
characters without stretching; they are outside the stretchable-fence gate.

## Validation and limits

The declared corpus has 66 fixtures: 63 clipboard byte comparisons plus 3 clipboard exclusions;
strict OOXML has 63 structural/text/style/colour comparisons and 3 required
rejections. The big-family fixture independently checks all 16 commands and
all four fixed run sizes in both clipboard HTML and DrawingML.
Run grouping, control formatting and dialect-specific vertical layout are
normalized; byte identity is not claimed. Independent golden annotations
check functions, upright/bold styles and enclosing nary fences.

The isolated Eqnedit64 CI lane runs CLI atomic-error checks, real native-backed
Markdown-to-PPTX packaging and `python-pptx` reopen checks. The Open XML SDK
checks shipped math fragments without renaming text or formatting, outer
package structure and a strict DrawingML run profile. DrawingML run-property
SDK limitations are counted separately; text is never rewritten to make validation pass. It does not completely
validate the `a14` payload, and these checks do not certify PowerPoint import.
The result JSON explicitly records that limitation.

Actual PowerPoint open/edit/save/reopen in an explicitly approved Office
session is a separate mandatory release gate. A missing session is **NOT RUN**,
not a pass. An approved desktop test uses a separate test presentation,
protects the user's documents, and changes no Office settings. Publication
and retirement-PR merge require their own approvals.

Primary format reference: [Microsoft DrawingML math extension](https://learn.microsoft.com/en-us/openspecs/office_standards/ms-odrawxml/853b19c7-68a9-4f9a-a2ae-5e6cb0d02e62).

Clipboard HTML explicitly declares English for mathematical/control runs and
Japanese for runs containing Japanese kana, CJK ideographs or full-width text.
This prevents inherited unrelated proofing languages without changing Office settings.

Native array column specifications (`\begin{array}{...}`) are unsupported;
use the supported `cases` or `matrix` environments. A single stretchable
delimiter encloses only one balanced group; adjacent groups remain separate.

## Desktop evidence and renewed candidate gate (2026-10-15)

Prior approved desktop batches opened, edited, saved and reopened the generated
math decks. Their saved math contained canonical `m:t` text (12 and 11 text
nodes in the two reviewed decks), replacing the producer's former `a:t`.
That evidence motivated this correction; it does not qualify the changed
candidate. The new candidate requires a fresh approved PowerPoint batch and
clipboard recheck, exact-SHA CI, signed O: hand testing and Claude review
before release. Those renewed desktop and signing gates are pending.
