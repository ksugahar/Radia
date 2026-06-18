---
name: panel-preview
description: Render any Radia Layer 3 PySide6 panel (radia_ih, radia_em, radia_pcb, radia_heat, radia_accel) on the real desktop Qt platform and save a screenshot PNG so a human or an AI assistant can VISUALLY review the current UI state. Use this skill when (a) you just edited a `radia_*.py` panel and want to confirm a new widget actually appears / the layout is not clipped, (b) the user asks "what does the X panel look like?", (c) you are about to ship a panel mode and need an end-to-end screenshot for the publish-panel checklist, (d) you are debugging which method exposes which widget. Default mode uses the real desktop Qt platform (font-metric accurate, readable text); pass `--offscreen` for CI / no-display use (NOTE: offscreen Qt on Windows renders all labels as ☐ box characters because the offscreen plugin does not load system fonts -- the companion .txt is the source of truth in that case). Auto-captures every method tab if no method is given. Output goes to `C:/temp/panel_preview/<panel>_<method>.png` by default.
---

# panel-preview

Render any Radia Layer 3 panel and save a screenshot.  Use this when
you need to **see** the panel — for verifying new widgets appear,
catching layout regressions, picking the right method tab to attach
to a doc, or just quickly answering "what does panel X look like?"
without launching Cubit.

## When to use

- After editing a `radia_*.py` panel and you want to verify your new
  widget actually appears in the right method (offscreen render is
  ~2s, much faster than launching Cubit + clicking through).
- The user asks "show me the IH panel", "what does PEEC+BEM look
  like?", "what widgets are visible in EM Kelvin mode?"
- Before saying "done" on a panel-mode addition (also runs as part
  of `publish-panel`).
- Debugging which widget belongs to which method (visibility wiring).

## When NOT to use

- For *behavioural* regressions (build_command output, save/restore,
  hidden widget feeding subprocess args) — that is `panel-qt-test`'s
  job and it runs assertions, not screenshots.
- For real-desktop font-metric clipping detection — that is the
  publish-panel skill's `--real-qt` block.  panel-preview can also do
  it via `--real-qt` but the publish-panel checklist also bundles
  vertical-clip measurement which this skill does not duplicate.

## Quick reference

Default mode = real desktop Qt (readable text in PNGs).  Captures
every method of every registered panel:

```bash
cd s:/Radia/01_GitHub
python .Codex/skills/panel-preview/preview.py
```

Single panel, single method:

```bash
python .Codex/skills/panel-preview/preview.py --panel radia_ih \
                                                --method "PEEC + BEM weak coupling (workpiece)"
```

CI / no-display (offscreen Qt; PNG fonts will be ☐ boxes on
Windows -- the companion .txt is the source of truth in that case):

```bash
python .Codex/skills/panel-preview/preview.py --panel radia_ih --offscreen
```

Specific output dir (default: `C:/temp/panel_preview/`):

```bash
python .Codex/skills/panel-preview/preview.py --panel radia_ih \
                                                --output-dir C:/temp/my_review/
```

## What gets captured

| Panel | Methods captured (default = all) |
|-------|----------------------------------|
| radia_ih | PEEC ind, BEM-A ind, PEEC+BEM, BEM-A+BEM, PEEC+FEM+Kelvin, Full FEM |
| radia_em | (panel-specific method enum) |
| radia_pcb | (panel-specific method enum) |
| radia_heat | (single mode, single screenshot) |
| radia_accel | (panel-specific method enum) |

The script auto-discovers `METHOD_*` constants in each panel module
to build the method list.

## Output

Each capture writes:

- `C:/temp/panel_preview/<panel>_<method-tag>.png` — screenshot
- `C:/temp/panel_preview/<panel>_<method-tag>.txt` — visible widget
   list (label + class name + visibility) for later string-grep
   review without re-launching Qt.

## Reading the output (AI assistant)

After running, list the PNGs and `Read` the ones relevant to the
user's question — Read renders PNG inline as multimodal content.
The companion `.txt` file gives a structured list of which widgets
the panel exposed for that method, which is faster than counting
pixels in the PNG.

```python
# example AI workflow:
# 1. user: "verify --wp-mesh-order widget appears in PEEC+BEM mode"
# 2. assistant runs preview.py --panel radia_ih
#                              --method "PEEC + BEM weak coupling (workpiece)"
# 3. assistant Read()s C:/temp/panel_preview/radia_ih_peec_bem.txt
#    and grep'd for "fes_order" / "Basis order"
# 4. assistant confirms visibility matches the build_command flag
```

## Implementation: preview.py

The actual runner is `preview.py` in this skill directory.  See its
docstring for full CLI surface; the SKILL.md only documents the
common cases.
