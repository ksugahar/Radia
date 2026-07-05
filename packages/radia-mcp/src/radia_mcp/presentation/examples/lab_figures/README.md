# Lab PowerPoint figure-creation scripts

5 historical scripts from `public-safe curated corpus`
(promoted into the presentation subpackage on 2026-05-26) that
illustrate the lab's PowerPoint-as-figure-editor workflow.

## Why PowerPoint as the figure editor?

The lab's drawing workflow:

1. **Generate** a `.pptx` programmatically (these scripts, or via
   `presentation_pptx_*` MCP tools that wrap python-pptx).
2. **Open** in PowerPoint and tweak interactively (move shapes,
   align, recolor, rotate, group, etc.).
3. **Export** the finished slide(s) to PNG / PDF for paper or
   poster embedding.

The advantage over matplotlib / Inkscape for diagrammatic figures
(circuit schematics, vector arrows, free-form shapes, annotated
mechanical layouts) is that PowerPoint:
  - Has a familiar UI everyone in the lab already knows
  - Snap-to-grid + smart guides + alignment tools work well
  - Saves as editable vector (.pptx) until the final export
  - Round-trips cleanly to embedded OOXML inside .docx papers

## Scripts (all use `win32com.client` — Windows + PowerPoint required)

| Script | Purpose |
|---|---|
| `lab_export_slide_to_png_com.py` | Single-slide → PNG/JPG export. The original recipe for "save my finished slide as a figure". |
| `lab_drawing_freeform_com.py` | Wide canvas (2000×500 pt) + AutoShape catalog (shape codes 1-99) + BuildFreeform sine wave. The "kitchen sink" demo. |
| `lab_marker_grid_com.py` | Plot-marker grid: circles, X markers, filled dots. Useful for ad-hoc legend / icon sheets. |
| `lab_sine_wave_com.py` | BuildFreeform path tracing a sine wave — the freeform pattern in minimal form. |
| `lab_color_grid_com.py` | 9×9 RGB grid + Export to PNG. Demonstrates the RGB-encoding gotcha (PowerPoint COM uses `r + g*256 + b*65536` — opposite of usual). |

## Cross-platform equivalent

The same patterns are available cross-platform (no PowerPoint COM)
via `presentation/pptx_drawing.py`, surfaced as MCP tools:

  - `presentation_pptx_create_canvas(width_mm, height_mm)`
  - `presentation_pptx_add_shape(...)` / `_add_line` / `_add_freeform`
  - `presentation_pptx_add_text_box(...)`
  - `presentation_pptx_save(pptx_path)`

Use the **python-pptx** MCP tools when:
  - Reproducible in CI / headless server
  - Cross-platform (macOS / Linux contributors)
  - You will hand the .pptx to a human for final touch-up

Use the **win32com COM** scripts here when:
  - Interactive — you want PowerPoint visible as the script runs
  - You need PNG/PDF export via PowerPoint's renderer (best fidelity)
  - You are on Windows AND have PowerPoint installed
