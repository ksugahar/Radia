# Changelog

## 3.0.10 - 2026-08-30

- Make ordinary PowerPoint Ctrl+V create editable, left-aligned inline Office
  Math at the accepted 18 pt destination size in both native and Web editions.
- Replace the false `Shapes.Paste()` acceptance surrogate with PowerPoint's
  built-in UI Paste command, and check 18 pt equation/tail/insertion ranges,
  inline OMML structure, and nonblank rendered output.

## 3.0.9 - 2026-08-30

- Show the semantic product version beside the source build stamp in the
  native title bar and About/version dialogs.
- Extend the real-PowerPoint acceptance test through the zero-length insertion
  point after the equation. Simple, powered, fractional, and full structural
  equations must all remain left-aligned editable math at 24 pt.

## 3.0.8 - 2026-08-29

- Keep the invisible inline PowerPoint sentinel at 24 pt in both the native
  and Web clipboard paths. The equation, trailing caret, saved inline OMML,
  and rendered result now share the same 24 pt contract.
- Extend the real-PowerPoint release gates to inspect the last character and
  the saved OOXML sentinel run, preventing an 18 pt caret from passing as a
  24 pt equation.

## 3.0.7 - 2026-08-29

- Uses the same inline 24 pt MathML plus trailing NBSP clipboard route in the
  native and browser editions; normal native copy no longer publishes the
  registered MathML formats that PowerPoint prioritizes as a centered display
  equation.
- Normalizes browser overlines, underlines, sums, and integrals to the native
  MathML contract. Real PowerPoint now saves byte-identical inline `m:oMath`
  and renders byte-identical PNG output for the native and browser fixtures.
- Strengthens the Office gate to reject `m:oMathPara`, direct browser OMML,
  centered ink, missing structures, and any size other than exactly 24 pt.

## 3.0.6 - 2026-08-29

- Preserves the browser editor's requested left alignment while making its
  editable PowerPoint equation exactly 24 pt through a synchronous CF_HTML
  fragment with conditional OMML and the shared fallback MathML.
- Tightens the homepage release gate to reject the former 18 pt browser
  conversion and any nested HTML document that causes Office style loss.

## 3.0.5 - 2026-08-29

- Replaces the native Geometry tab's accidental over/underline mapping with
  the same 11-item differential-geometry palette used by the browser editor,
  while retaining the dedicated native over/underline palette under Basic.
- Normalizes browser and native fallback MathML to inline 24 pt output, uses
  TeX/MathJax side limits for integrals, and adds a browser-only conditional
  OMML transport so PowerPoint preserves sums, integrals, fractions, radicals,
  overlines, and underlines as editable Office Math.
- Adds cross-product palette parity and real PowerPoint OOXML/render regression
  coverage so a merely non-empty equation no longer qualifies as parity.

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
