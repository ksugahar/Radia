---
name: ppt2pdf
description: Convert pptx to PDF for inclusion in LaTeX PDF output. Uses PowerPoint COM automation.
---

# ppt2pdf — PowerPoint to PDF Converter

Convert `.pptx` / `.ppt` files to `.pdf` using PowerPoint COM automation (Windows only).
Output PDF can be directly included in LaTeX via `\includegraphics`.

## Usage

```bash
python s:/Radia/01_GitHub/.Codex/skills/ppt2pdf/ppt2pdf.py <file.pptx> [file2.pptx ...]
python s:/Radia/01_GitHub/.Codex/skills/ppt2pdf/ppt2pdf.py <directory>
```

- Specific file(s): converts only the given file(s)
- Directory: converts all pptx files in that directory
- Output PDF is saved alongside the source file with the same base name

## Dependencies

- Windows with PowerPoint installed
- `pywin32` (`pip install pywin32`)

## Example

```bash
# Convert a specific pptx file (for LaTeX figure)
python s:/Radia/01_GitHub/.Codex/skills/ppt2pdf/ppt2pdf.py "Figures/agent_crosssection.pptx"

# Convert multiple files
python s:/Radia/01_GitHub/.Codex/skills/ppt2pdf/ppt2pdf.py fig1.pptx fig2.pptx

# Convert all pptx in a directory
python s:/Radia/01_GitHub/.Codex/skills/ppt2pdf/ppt2pdf.py "Figures/"
```
