# Changelog

## 3.0.4 - 2026-08-29

- Adds an always-visible math-style palette to both the native and browser
  editors for upright, italic, bold, sans-serif, monospace, calligraphic,
  blackboard-bold, Fraktur, and bold-math input.
- Supports `\mathnormal`, `\mathrm`, `\mathit`, `\mathbf`, `\mathsf`,
  `\mathtt`, `\mathcal`, `\mathbb`, `\mathfrak`, `\bm`, and `\boldsymbol`
  across native parsing, rendering, TeX/MathML export, and palette insertion.
- Uses the embedded math font's Unicode mathematical alphabets, including bold
  Greek glyphs, and keeps browser `\bm` input compatible with MathJax.
- Adds a C++ regression gate for stable math-alphabet normalization.

## 3.0.3 - 2026-08-29

- Draws `\mathbf` and vector-style characters with the embedded math font's
  designed Unicode bold glyphs instead of an indistinguishable regular glyph.
- Replaces misleading numeric, unsupported-slot, combining-mark, and CJK
  palette faces with compact glyphs owned by the embedded font.
- Extends the hidden visual gate to compare plain/bold GDI pixels and render
  every owner-drawn palette cell at 96, 120, 144, and 192 dpi.
- Rewrites the product handover as the current TeX-only, standalone,
  separately released Eqnedit64 operating contract.

## 3.0.2 - 2026-08-29

- Replaces crash-prone in-memory font registration with a verified,
  content-addressed per-user font cache and file-backed private registration.
- Adds an isolated 32-lifecycle CI gate that fails if Windows restarts
  `fontdrvhost.exe` or records a new font-host application crash.
- Requires native category-tab labels to occupy a readable pixel height; a few
  stray pixels no longer qualify as visible toolbar text.

## 3.0.1 - 2026-08-29

- Owner-draws native palette tabs and buttons so correct window text cannot
  become invisible when a Windows GDI/theme session stops painting BUTTON text.
- Adds a hidden pixel-difference regression test over every real toolbar button.
- Measures and draws Japanese text with an explicit CJK system face, preventing
  math-font linking from expanding nominal 12 pt text several times over.

## 3.0.0 - 2026-08-29

- First PyPI packaging of the Radia-owned Eqnedit64 product.
- Ships the structural `eqnedit_core` Python extension, the matching signed
  standalone Windows executable, and the browser editor assets.
- Provides checked Python wrappers for Office/PowerPoint and Google Slides
  clipboard paths plus PNG/EMF rendering.
- Uses UTF-8 TeX exclusively; MTEF and `.eqn` are not supported.
