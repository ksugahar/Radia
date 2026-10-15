# Eqnedit64 for Python

Source version: **3.1.4** (review candidate; released version remains 3.1.3).

Candidate native OOXML export: `Eqnedit64.exe equation.tex equation.pptxmath`.
The checked Python helper is `eqnedit64.export_powerpoint_math(tex, output)`.
The optional radia-mcp adapter assembles Markdown decks with native editable math;
missing/old backends and unsupported input raise. See the native product's
`docs/POWERPOINT_MATH.md` for the format, supported subset and validation limits.
This candidate has not passed the isolated PowerPoint release gate and is not published.

`eqnedit64` packages the same TeX-first equation engine and standalone Windows
application maintained in the Radia monorepo. The signed one-file
`Eqnedit64.exe` is the canonical editor and command-line converter and remains
available separately from GitHub Releases. It does not require Python.

```powershell
python -m pip install eqnedit64
eqnedit64
```

The direct native command-line converter always takes one input and one output:
`Eqnedit64.exe INPUT OUTPUT`. The package's `eqnedit64` console command is an
optional launcher for the bundled copy of that same EXE, not another CLI and
not a pybind11 surface. The input is
a UTF-8 TeX file or `clipboard`; the output is `office`, `slides`,
`clipboard-png`, or a `.png` / `.emf` file path. The old CLI output word `png` and the compatibility switches are rejected
with exit code 94. Use `clipboard-png` for clipboard PNG or a `.png` file path.
The descriptive destinations `powerpoint` / `word` and `google-slides` remain supported.

```powershell
Eqnedit64.exe equation.tex office
Eqnedit64.exe equation.tex equation.png
Eqnedit64.exe clipboard clipboard-png
```

After a wheel install, replacing `Eqnedit64.exe` above with `eqnedit64` runs
the packaged copy with the same arguments.

The wheel provides:

- `eqnedit64.Equation` for structural editing, Tab-slot traversal, undo/redo,
  Backspace, multiline equations, and arbitrary matrix growth;
- `tex_to_svg`, `tex_to_mathml`, `tex_normalize`, palettes, and shortcuts;
- `copy_equation()` and `render_equation()` as thin subprocess users of the
  bundled native CLI;
- the Radia-owned browser editor assets through `web_asset()`;
- the bundled standalone application through the `eqnedit64` console command.

Only UTF-8 TeX is a source format. Current wheels target 64-bit Windows and Python 3.10 or newer.

The wheel embeds the signed executable produced for the matching
`eqnedit64-v3.1.3` GitHub Release; it does not build a second application
implementation. The package is BSD 2-Clause licensed. The embedded Eqnedit
Math font is a TrueType-outline derivative of Latin Modern Math under the
GUST Font License. Its license and derivation manifest are included in the wheel.


## 3.1.2 PowerPoint colour limitation (UXP-0030/0031)

PowerPoint import drops colour on the shared conditional OMML CF_HTML route in
both editions. The user measured all r10 colour variants A–G as not red on
2026-10-08 and approved releasing 3.1.2 with this documented limitation. Colour
support remains enabled in the native/Web canvas, SVG, saved TeX, and native
EMF/DIB output. This decision does not strip colour from the source or renderer.
PowerPoint colour retention is not an acceptance requirement for this release.

User hand tests: H1/H5/H6 EXE OK on ea787664d; H7 EXE canvas red, PowerPoint
not red. Web H5/H6 OK on 8361ec916. Historical CI 37692902458 compared 22/22 primary OMML branches at ea787664d;
its Web side used the headless browser clipboard, not the OS clipboard.
This supports the UXP-0031 EXE-test plus bit-gate qualification for the covered
corpus; it does not prove every TeX input or every Office/browser environment.
Changed payload behaviour requires renewed qualification. The completed 3.1.3
publication record is in the [canonical handover](../../docs/equation/eqnedit64_handover.md).


## 3.1.3 release and coverage limits

Version **3.1.3** is released. The 52-fixture Office corpus has 49 covered
native/Web primary OMML branches that are byte-identical. Headed full Chromium
run [37721927723](https://github.com/ksugahar/Radia/actions/runs/37721927723)
verified Web OS-clipboard writes, styles and functions on `25fca594e`.
The user accepted six signed-candidate retests, including function arguments
and upright/bold styles. Qualification applies to the covered corpus.

Three constructs remain non-covered: native `\boxed` and `\binom` drop the
unsupported command, and native `\cancel` emits a bottom bar instead of Web's
diagonal strike. They are not passes and do not qualify for Web hand-test
exemption. Bare functions may have empty OMML arguments; `sin 2x` collects only
`2`. PowerPoint colour loss on the shared OMML route remains a known limitation.
Full arbitrary-TeX or Office/browser parity is not established.

## PowerPoint normal paste

Ctrl+V must produce editable inline m:oMath on the slide (no m:oMathPara there), left aligned, with 18 pt equation, trailing character and insertion point in a blank presentation. EXE/Web put RichEdit OMML in the CF_HTML msEquation block as the primary payload, with MathML only as fallback. Test Application.CommandBars.ExecuteMso("Paste"), never substitute Shapes.Paste(). Refuse attachment to user processes, inspect saved OMML/rendering, restore clipboard and close owned resources only. Slides image copy remains 300 dpi / 24 pt.
