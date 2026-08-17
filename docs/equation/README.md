# `radia.equation` — equations as LaTeX

An equation is stored as **LaTeX**, usually inside a Markdown file. Everything
else is produced from that: a native Office equation for Word and PowerPoint,
an SVG picture, and a layout an editor can draw with and position a caret from.

```python
import radia.equation as eq

eq.tex_to_omml(r"\frac{a+b}{c}")     # Office-native math, for a .docx / .pptx
eq.tex_to_rtf(r"\frac{a+b}{c}")      # the same equation, for the clipboard
eq.tex_to_svg(r"\frac{a+b}{c}")      # a picture
eq.markdown_to_docx(open("note.md", encoding="utf-8").read(), "note.docx")
```

## Why OMML rather than a picture

An equation pasted as OMML **is an equation**: Word's and PowerPoint's own
tools edit it, it follows the theme font and colour, it scales with the text,
and the reader needs nothing installed. A picture is none of those things, and
an OLE object needs Equation Editor on the reader's machine.

## Three spellings, one walk

Office takes an equation in three spellings, and which one is needed depends on
where it is going. All three were measured by copying an equation out of Word
and reading the clipboard, rather than inferred:

| | | |
|---|---|---|
| OMML | inside a `.docx` / `.pptx` | `<m:f><m:num/><m:den/></m:f>` |
| RTF | clipboard, **Word only** | `{\mf{\mfPr{\mctrlPr}}{\mnum ...}{\mden ...}}` |
| MathML | clipboard, **PowerPoint and Excel** | `<mfrac><mi>a</mi><mi>b</mi></mfrac>` |

The split is not a preference, it is measured. Offering Word's formats to
PowerPoint one at a time: with only `MathML` present PowerPoint produces a
native equation; with only `Rich Text Format` it produces a **text box**. Word
is the other way round -- it reads our RTF and refuses our MathML alone.

That is why `copy_to_clipboard` puts both on at once. The Windows clipboard
holds many formats and each application takes the richest thing it understands,
so one Copy serves every target with no mode switch.

Since the structure is the same in every case and only the spelling differs, the
interesting part is written once in `math_writer.cpp` and each output supplies a
`MathSyntax`. The interface names *meanings* (fraction, radical, script) rather
than element names: OMML and RTF happen to share names, but MathML does not --
its root takes its arguments the other way round, its delimiters are ordinary
operators inside a row, and its runs split into identifier / number / operator.

Duplicating the walk would guarantee the outputs drift apart, and the walk is
where the judgement lives: MTEF stores a script's base as the *preceding*
sibling and a big operator's operand as the *following* run, so both have to be
absorbed or Office draws a placeholder box.

`validation_test/equation/test_paste_into_word.py` is the acceptance test:
clipboard -> paste -> save -> read back `<m:oMath>`, for Word and PowerPoint.
An equation that arrives as text runs is a picture at best, and no amount of
correct-looking markup would show that.

### Pictures, and where they are needed

| target | what it takes |
|---|---|
| Word | RTF, a native equation. No picture. |
| PowerPoint | MathML, a native equation. No picture. |
| Excel | equations live in a floating object; a pasted Word equation lands as a picture, and Excel keeps the **metafile** beside the bitmap |
| Google Slides | no equation object at all. SVG upload is rejected; raster only. EMF reaches it as vector through Google Drawings. |

`tex_to_emf` and `tex_to_png` come from one GDI drawing routine over one
layout -- the same layout an editor will draw on screen with and position its
caret from. In GDI a metafile, a bitmap and a window are all device contexts,
so there is no second renderer and no separate PNG encoder to keep in step.

The pictures are on the clipboard *after* the equation formats deliberately.
"The richest format wins" is not a rule -- each application picks by its own
priority list, so an equation format can lose to a metafile. Word and
PowerPoint were re-measured with the pictures present and still produce native
equations (5/5 each), but `copy_to_clipboard(..., pictures=False)` exists for
an application that would rather have the picture.

**Excel is not verified.** Its equations live in shapes, and pasting into a
shape's text is a UI operation its object model does not expose -- Shape.Select
plus Paste, TextRange.Select plus Paste, and PasteSpecial are all refused. The
payload is the same one PowerPoint accepts and Excel uses the same equation
object, so it is likely to work; that is not the same as measured.

Going "via Word" does not help, and that *is* measured: copying an equation in
Word and pasting into Excel produces a picture (`image1.png` + `image2.emf`),
not an equation, even though Word's clipboard carries strictly more formats
than ours. The missing piece is the paste target, not the format.

## Layout

| | |
|---|---|
| `src/ext/equation/tex_parser.cpp` | LaTeX → node tree |
| `src/ext/equation/math_writer.cpp` | the tree walk every Office output shares |
| `src/ext/equation/mtef_omml.cpp` | OMML spelling (inside a .docx / .pptx) |
| `src/ext/equation/mtef_rtf.cpp` | RTF spelling (clipboard, Word) |
| `src/ext/equation/mtef_mathml.cpp` | MathML spelling (clipboard, PowerPoint / Excel) |
| `src/ext/equation/math_layout.h` | the display list every way of drawing shares |
| `src/ext/equation/mtef_svg.cpp` | node tree → layout, and layout → SVG |
| `src/ext/equation/mtef_gdi.cpp` | layout → EMF and PNG, through GDI |
| `src/ext/equation/eq_edit.cpp` | the editing model |
| `src/ext/equation/md_doc.cpp` | which spans of a `.md` are math |
| `src/radia/equation/office.py` | Word and PowerPoint writers |

## The editing model

`Equation` is a model, not a widget. It carries an insertion point that lives
*inside* the structure, so an arrow key walks into a fraction rather than
skipping over it; templates whose empty slots are reached with Tab; and a
backspace that at the start of an empty slot unwraps the template and splices
its contents into the parent instead of swallowing them.

```python
e = eq.Equation()
e.insert_text("x")
e.command("template.sub")      # Ctrl+L: the x just typed becomes the base
e.insert_text("i")
e.latex()                      # 'x_{i}'
```

A front end binds keys to command names; the model knows nothing about
keyboards. `Equation.shortcuts()` returns Equation Editor 3.0's chords as
`(chord, command, label)` so the familiar ones survive into a new front end
rather than being reinvented. Only chords that are unambiguous in that
editor's own Help are listed; the rest of `Equation.templates()` is reachable
from a palette and is deliberately left unbound.

Undo snapshots the whole LaTeX rather than inverting each command. Equations
are tens of characters, so a snapshot costs nothing and cannot drift out of
step with the tree the way command inversion can. It rests on
LaTeX → tree → LaTeX reaching a fixed point, which `tex_normalize` exposes and
`tests/equation/test_edit.py` checks.

## Markdown

`MarkdownDoc` answers one question: which spans of the file are math.
Headings, emphasis and lists are a text editor's usual business and are
deliberately not modelled.

Loading and saving an untouched file reproduces it **byte for byte**,
delimiters included, so a file written with `\(` `\)` does not come back as
`$` `$`. `set_math_latex(i, ...)` replaces one equation and leaves every other
byte alone.

The scan is conservative in one direction on purpose. Missing an equation
costs the user one manual edit; inventing one rewrites their prose or their
code when the file is saved, silently. So a `$` inside a fenced block or an
inline code span is never math, `\$` is a dollar, the inline rule follows
Pandoc's — no space just inside the delimiters, no digit after the closing one
— and inline math may not span a backtick. Without that last rule, a sentence
like *prices $5 and $6, where `$HOME` is set* has its second dollar open a span
that closes on the dollar inside the code span.

## MTEF

MTEF is the Equation Editor 3.x / MathType binary format. It is supported in
one direction for one reason: to read equations out of documents that already
contain them. `mtef_to_tex` and `mtef_to_omml` are that import path, and
`tex_to_mtef` exists so an equation can be handed back to Equation Editor as a
`.eqn` file. Nothing else here goes through MTEF, and its fidelity is not a
goal.

`dump_tree(data, run_passes=...)` and `tex_dump_tree(latex)` print the parsed
tree. Every disagreement between the LaTeX, OMML and SVG paths so far has come
from reading a different field of that tree, and printing it settles the
question in seconds.

## Known limits

- `\sum`'s stacked-limits flag has nowhere to live in LaTeX, so it does not
  survive a save. Display style stacks them anyway.
- An empty row has no LaTeX spelling, so a blank trailing row of a `cases` or
  a matrix disappears on save. Filled rows round-trip exactly.
- A slashed fraction is not offered as a template: it is written as MTEF's
  `{}^{a}/{}_{b}`, which does not read back as one fraction, and a template
  that changes shape when saved is worse than no template.
- Spacing commands (`\,` and friends) emit nothing; Office's own spacing is
  used instead.
- The SVG layout implements TeX's inter-atom spacing and per-glyph boxes, not
  the whole of TeX. Font metrics come from GDI, so it is Windows-only; the
  other paths are not.
