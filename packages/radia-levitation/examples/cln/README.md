# CLN — Cauer Ladder Network for eddy-current modeling

Working folder for the **Cauer Ladder Network (CLN)** research line: extracting
RL-ladder equivalent circuits of eddy-current decay in 3D conductors via
double-double (DD) verified arithmetic, BEM-Foster decomposition, and the
Hankel-Padé / Stieltjes / Wheeler family of moment-matching methods.

This example folder is the canonical working location (since 2026-05-12) for:

- **IGTE 2026 Symposium digest** (`igte_symposium_2026.tex` + `igtesymp.cls`) —
  sphere-only paper showing how DD GPU pipeline pushes the Cauer-extraction
  precision wall from FP64 stage 4-5 to DD stage 12+.
- **Q&A development log** (`MEMORY.md`, 600+ lines) — running record of design
  decisions during digest preparation, kept in sync with the figures and
  numerical tables.
- **Cross-references** to the radia-mcp tool `cln_sphere_dd_pipeline`
  (under `packages/radia-mcp/.../radia_ngsolve/knowledge/cln_sphere_dd.py`)
  which exposes the same pipeline as an LLM-callable tool.

## Files

| File | Purpose |
|---|---|
| `igte_symposium_2026.tex` | IGTE 2026 digest LaTeX source (2 pages, igtesymp class) |
| `igte_symposium_2026.pdf` | Compiled digest (current state) |
| `igtesymp.cls` | IGTE Symposium 2026 LaTeX class file (provided by organisers) |
| `MEMORY.md` | Q&A development log + design decisions |

## How it relates to the rest of the repo

- **`packages/radia-mcp/src/radia_mcp/radia_ngsolve/knowledge/cln_sphere_dd.py`**
  — markdown doc describing the DD pipeline; surfaced via the `cln_sphere_dd_pipeline` MCP tool.
- **`src/ext/axifemm/`** — C++ Henrotte axisymmetric FE (Q1/Q2/P1/P2) used for
  the FEM cross-check (axifemm_p2_triangle Phase B2; commit `81f6415f`).
- **Other CLN literature**: Cauer 1958, Henrici 1958 QD-Padé, Sugahara TEAM 28
  axisymmetric matlab, Stoll Bessel ground truth, Hiruma 3-term FEM-CLN.

## Status

- IGTE digest: sphere-only pivot complete (2026-05-12 long discussion with
  Prof. Nagamine). Cuboid Outlook section dropped in favour of pure
  sphere + DD-pipeline demonstration.
- Pending: reflect DD 540-cell results (`stage 0-5 reliable`, +1 stage over 270-cell baseline) into the digest tables.

## Build the digest

```bash
cd examples/CLN
pdflatex igte_symposium_2026.tex   # or via texcompile skill
```
