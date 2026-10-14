# Eqnedit64 PowerPoint normal-paste contract

## Accepted decision (2026-08-30)

Eqnedit64.exe and the Web/JS editor use the same normal-copy result in
Microsoft PowerPoint:

- editable native Office Math;
- inline `m:oMath`, never centred `m:oMathPara`;
- left-aligned paragraph;
- 18 pt for the equation, trailing sentinel, and next insertion point in a
  standard blank PowerPoint presentation.

Left alignment has higher product priority than forcing 24 pt.  Google Slides
image copy remains a separate 300 dpi / 24 pt contract.

## Why the former 24 pt result was real but unsuitable

Eqnedit64 3.0.6 published the registered Windows clipboard formats `MathML`
and `MathML Presentation`. PowerPoint gave those formats priority, retained
`mathsize="24pt"`, and stored a standalone paste as `m:oMathPara`. The visible
result was a centred display equation. This native Windows route cannot be
made the common EXE/Web contract because the Web Clipboard API cannot publish
the same registered native format.

From 3.0.7 onward both editions used inline MathML in CF_HTML. This produces
the desired editable, left-aligned `m:oMath`, but PowerPoint's ordinary UI
Paste applies the destination formatting and produces 18 pt. `mathsize`, CSS,
`mstyle`, a styled NBSP, and Unicode-text MathML did not make ordinary UI Paste
retain 24 pt. `PasteSourceFormatting` did retain 24 pt, but requires a
different paste command and is not the normal Ctrl+V workflow.

OOXML/OMML is PowerPoint's saved representation after import, not a documented
public clipboard interchange format. Office.js `Ooxml` coercion is supported
for Word, not for PowerPoint.

Primary references:

- Microsoft MathML clipboard/import rules:
  <https://learn.microsoft.com/en-us/office/math/mathml>
- W3C Clipboard API mandatory and Web custom formats:
  <https://www.w3.org/TR/clipboard-apis/>
- Microsoft Paste Special and destination/source formatting:
  <https://support.microsoft.com/en-us/word/paste-special>
- Office.js coercion support (`Ooxml` is Word-only):
  <https://learn.microsoft.com/en-us/javascript/api/office/office.coerciontype?view=powerpoint-js-preview>

## Regression-test rule

Never use `slide.Shapes.Paste()` as evidence for the user's Ctrl+V result.
That PowerPoint object-model path retained 24 pt while the built-in UI Paste
command on the exact same clipboard produced 18 pt. It caused the 3.0.8/3.0.9
false-positive release evidence.

The acceptance path is:

1. refuse to attach if a user PowerPoint process already exists;
2. create a temporary presentation window and move it offscreen;
3. select its blank slide;
4. invoke `Application.CommandBars.ExecuteMso("Paste")`;
5. verify 18/18/18 pt (equation/final character/zero-length insertion point),
   left alignment, one MathZone, inline `m:oMath`, no `m:oMathPara`, and a
   nonblank PowerPoint-rendered image;
6. close only the test presentation and restore every clipboard format.

`Shapes.Paste()` may be tested as an API of its own, but must never substitute
for the normal-paste product gate. The Web publication QA must use the same
built-in Paste command and compare its saved OMML and rendered result with the
native edition.
