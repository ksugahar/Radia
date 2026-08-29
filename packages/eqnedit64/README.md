# Eqnedit64 for Python

`eqnedit64` packages the same TeX-first equation engine and standalone Windows
application maintained in the Radia monorepo. It is intended for Python, MCP,
and command-line workflows; the signed one-file `Eqnedit64.exe` remains
available separately from GitHub Releases for machines without Python.

```powershell
python -m pip install eqnedit64
eqnedit64
```

The wheel provides:

- `eqnedit64.Equation` for structural editing, Tab-slot traversal, undo/redo,
  Backspace, multiline equations, and arbitrary matrix growth;
- `tex_to_svg`, `tex_to_mathml`, `tex_normalize`, palettes, and shortcuts;
- `copy_equation()` and `render_equation()` over the bundled native CLI;
- the Radia-owned browser editor assets through `web_asset()`;
- the bundled standalone application through the `eqnedit64` console command.

Only UTF-8 TeX is a source format. MTEF and `.eqn` are intentionally not
supported. Current wheels target 64-bit Windows and Python 3.10 or newer.

The wheel embeds the signed executable produced for the matching
`eqnedit64-v3.0.10` GitHub Release; it does not build a second application
implementation. The package is BSD 2-Clause licensed. The embedded Latin
Modern Math font remains under the GUST Font License documented in the Radia
source tree.
