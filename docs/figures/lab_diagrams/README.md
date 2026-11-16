# Lab diagram templates (worked examples)

Ready-to-build, version-controlled source for the lab's canonical **flowcharts and
conceptual diagrams**. Text source (diffable), labels inherit the paper's Times/newtx
font, vector output — the lab default for a paper diagram (TikZ).

| File | Kind | Renders |
|------|------|---------|
| `kelvin_dtn_solve_flowchart.tex` | flowchart (ISO-5807 shapes) | the Kelvin/DtN open-boundary solve loop |

## Build

```bash
pdflatex kelvin_dtn_solve_flowchart.tex     # -> .pdf (standalone, crop-tight)
```

Embed in a paper at the column width with `\includegraphics[width=<W>cm]{...}` (the title
goes in the LaTeX `\caption{}`, never inside the figure — lab rule).

Diagram design follows ISO 5807:1985 shape semantics, one flow direction, minimal edge
crossings and a label on every decision branch. Use Graphviz (DOT + `dot -Tpdf`) instead
for an auto-layout flowchart of a large graph.

> The PDFs are derived artifacts (regenerate from the `.tex`); only the source is tracked.
