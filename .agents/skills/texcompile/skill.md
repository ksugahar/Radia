---
name: texcompile
description: Compile LaTeX to PDF with automatic file lock release. Kills PDF viewer processes locking the output, cleans corrupted aux files, runs pdflatex + bibtex.
---

# texcompile — LaTeX to PDF with Lock Release

Compile `.tex` files to PDF, automatically releasing file locks held by PDF viewers.

## Usage

```bash
python s:/Radia/01_GitHub/.Codex/skills/texcompile/texcompile.py <file.tex>
```

## What it does

1. Detects if the output PDF is locked by another process
2. Identifies and closes PDF viewers (Acrobat, SumatraPDF, Foxit, etc.) holding the lock
3. Cleans corrupted aux/out files
4. Runs full pdflatex + bibtex pipeline (pdflatex -> bibtex -> pdflatex -> pdflatex)
5. Reports page count and file size

## Dependencies

- TeXLive (auto-detected from `C:\texlive\`)
- Windows (uses PowerShell for process management)

## Example

```bash
python s:/Radia/01_GitHub/.Codex/skills/texcompile/texcompile.py paper.tex
```
