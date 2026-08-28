# Changelog

## 3.0.1 - 2026-08-29

- Owner-draws native palette tabs and buttons so correct window text cannot
  become invisible when a Windows GDI/theme session stops painting BUTTON text.
- Adds a hidden pixel-difference regression test over every real toolbar button.

## 3.0.0 - 2026-08-29

- First PyPI packaging of the Radia-owned Eqnedit64 product.
- Ships the structural `eqnedit_core` Python extension, the matching signed
  standalone Windows executable, and the browser editor assets.
- Provides checked Python wrappers for Office/PowerPoint and Google Slides
  clipboard paths plus PNG/EMF rendering.
- Uses UTF-8 TeX exclusively; MTEF and `.eqn` are not supported.
