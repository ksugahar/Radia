# Native PowerPoint math export

Version 3.1.4 is a review candidate, not a published release.

```text
Eqnedit64.exe equation.tex equation.pptxmath
```

The output is an XML `m:oMathPara` in the 2006 Office math namespace, with
DrawingML `a:rPr` and `a:t` runs. It is intended for PowerPoint's `a14:m`
extension, not Word import and not HTML. The serializer shares equation
structure with the normal-copy RichEdit HTML route; the latter is unchanged.
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

## Validation and limits

The declared corpus has 52 fixtures: 49 structural/text/style/colour semantic
comparisons against the clipboard emitter, and 3 required rejections.
Run grouping, control formatting and dialect-specific vertical layout are
normalized; byte identity is not claimed. Independent golden annotations
check functions, upright/bold styles and enclosing nary fences.

The isolated Eqnedit64 CI lane runs CLI atomic-error checks, real native-backed
Markdown-to-PPTX packaging and `python-pptx` reopen checks. The Open XML SDK
checks ISO math structure after validation-only run normalization, outer
package structure and a strict DrawingML run profile. It does not completely
validate the `a14` payload, and these checks do not certify PowerPoint import.
The result JSON explicitly records that limitation.

Hidden PowerPoint open/save/render in an approved isolated Office session is
a separate mandatory release gate. A missing session is **NOT RUN**, not a
pass. Do not run on an interactive desktop or reuse a user's active Office
process. Publication and retirement-PR merge require their own approvals.

Primary format reference: [Microsoft DrawingML math extension](https://learn.microsoft.com/en-us/openspecs/office_standards/ms-odrawxml/853b19c7-68a9-4f9a-a2ae-5e6cb0d02e62).
