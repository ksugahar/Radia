# Eqnedit64 command line

Eqnedit64 is both a standalone editor and a one-step TeX converter. The public
converter syntax always has one input and one output:

```text
Eqnedit64.exe INPUT OUTPUT
```

No Python runtime, installer, ActiveX registration, or external TeX command is
required by the standalone executable.

The native executable is the canonical command-line boundary. The optional
Python package console launcher and helper functions invoke this same EXE; they
do not define a second converter. The CLI dispatcher is intentionally not
duplicated through pybind11. The existing pybind11 structural-editing API is a
separate optional interface and is not required for conversion.

## Input

| `INPUT` | Meaning |
|---|---|
| `equation.tex` | UTF-8 TeX file |
| `clipboard` | Registered `LaTeX` clipboard data, then Unicode text as fallback |

Outer `$...$`, `$$...$$`, `\(...\)`, `\[...\]`, `equation`, and `equation*`
wrappers are accepted. Paper metadata such as `\label`, `\nonumber`, `\notag`,
`\tag`, and unescaped `%` comments does not become visible equation content.

## Output

| `OUTPUT` | Result |
|---|---|
| `office` | Editable PowerPoint/Word equation plus TeX, EMF, and opaque DIBV5 clipboard fallbacks |
| `slides` | Google Slides 300 dpi / 24 pt PNG and HTML clipboard data |
| `clipboard-png` | PNG and opaque DIBV5 clipboard data; no file is created |
| `equation.png` | PNG file at the named path |
| `equation.emf` | EMF file at the named path |

`powerpoint` and `word` are aliases for `office`; `google-slides` is an alias
for `slides`. The old output word `png` remains accepted as an alias for
`clipboard-png`, but new scripts and teaching material use the explicit name so
it cannot be mistaken for a file path.

```powershell
Eqnedit64.exe equation.tex office
Eqnedit64.exe equation.tex slides
Eqnedit64.exe equation.tex clipboard-png
Eqnedit64.exe equation.tex equation.png
Eqnedit64.exe equation.tex equation.emf
Eqnedit64.exe clipboard clipboard-png
```

## Editor launch and diagnostics

- No arguments start the editor.
- One `.tex` path opens that document.
- If Explorer supplies two `.tex` paths, Eqnedit64 opens the first document and
  reports the choice in the status bar. It never treats the second `.tex` as an
  output image or overwrites it.
- An invalid conversion output writes a diagnostic to standard error when a
  console or redirected pipe is available; otherwise it displays an error
  dialog. The exit code remains 94.
- `--help` and `--version` write UTF-8 text to standard output when invoked by a
  console or script. A graphical launch without a usable standard stream keeps
  the dialog behavior.

Exit codes are 0 for success, 82 for input read failure, 83 for empty input, 84
for clipboard publication failure, and 94 for an invalid output destination.
The pre-3.0.13 `--copy-*`, `--render-*`, and `--texclip` spellings remain
compatibility inputs only and are not part of the public command-line API.
