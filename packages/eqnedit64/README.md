# Eqnedit64 for Python

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
`clipboard-png`, or a `.png` / `.emf` file path. The legacy word `png` remains
accepted as a clipboard alias, but a path ending in `.png` always creates a
file.

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

Only UTF-8 TeX is a source format. MTEF and `.eqn` are intentionally not
supported. Current wheels target 64-bit Windows and Python 3.10 or newer.

The wheel embeds the signed executable produced for the matching
`eqnedit64-v3.0.14` GitHub Release; it does not build a second application
implementation. The package is BSD 2-Clause licensed. The embedded Latin
Modern Math font remains under the GUST Font License documented in the Radia
source tree.
