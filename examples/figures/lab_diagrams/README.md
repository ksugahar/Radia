# Lab diagram templates (worked examples)

Ready-to-build, version-controlled source for the lab's canonical **flowcharts and
conceptual diagrams** — worked examples of the `radia_mcp.figure` diagram-drawing skill
(`figure_diagram_recipes`). Text source (diffable), labels inherit the paper's Times/newtx
font, vector output — the lab default for a paper diagram (TikZ).

| File | Kind | Recipe topic | Renders |
|------|------|--------------|---------|
| `kelvin_dtn_solve_flowchart.tex` | flowchart (ISO-5807 shapes) | `figure_diagram_recipes("tikz_flowchart")` | the Kelvin/DtN open-boundary solve loop |
| `panel_4layer_architecture.tex` | concept / block diagram | `figure_diagram_recipes("concept_diagram")` | the Cubit panel 4-Layer architecture (`.vol` interface, subprocess boundaries) |

## Build

```bash
pdflatex kelvin_dtn_solve_flowchart.tex     # -> .pdf (standalone, crop-tight)
pdflatex panel_4layer_architecture.tex
```

Embed in a paper at the column width with `\includegraphics[width=<W>cm]{...}` (the title
goes in the LaTeX `\caption{}`, never inside the figure — lab rule).

## How these were authored

Via the `radia_mcp.figure` MCP server:

1. `figure_diagram_recipes("tool_selection")` — TikZ for a paper diagram (vs Graphviz for
   large auto-layout graphs, schemdraw for Python, Mermaid for READMEs).
2. `figure_diagram_recipes("tikz_flowchart")` / `("concept_diagram")` — the templates.
3. `figure_diagram_recipes("design")` — ISO 5807:1985 shape semantics, one flow direction,
   minimise edge crossings, label every decision branch.

For an auto-layout flowchart of a large graph, use `figure_diagram_recipes("graphviz")`
(DOT + `dot -Tpdf`).

> The PDFs are derived artifacts (regenerate from the `.tex`); only the source is tracked.
