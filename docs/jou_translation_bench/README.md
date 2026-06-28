# .jou → build123d translation benchmark

Purpose: measure how well the current mcp-server knowledge bases
(`src/radia/mcp_server/cubit/` and `.../build123d/`) let an LLM
translate Cubit `.jou` scripts to build123d Python.

This is a docs-owned benchmark record, not an `examples/` runtime sample.  The
2026-04-19 report's observed gaps have been folded into the public
`radia_mcp.build123d` and `radia_mcp.cubit` knowledge tables.

**This is not a runtime compatibility layer** — we are NOT executing
`.jou` on build123d. We are checking whether a fresh LLM session,
given only our lab knowledge, can produce idiomatic build123d code
equivalent to the CAD intent of a `.jou`, and whether it can refuse
translation when the `.jou` uses Cubit-only concepts (imprint/merge,
hex sweep, block/sideset).

## Fixtures (`01` → `03`, increasing difficulty)

| File | Difficulty | Expected verdict |
|---|---|---|
| `01_pure_cad.jou` | easy | TRANSLATE. primitives + boolean + move + rename only. |
| `02_sweep_revolve.jou` | medium | TRANSLATE. revolve + fillet; selector logic for fillet is the tricky bit. |
| `03_imprint_hex.jou` | hard | **SKIP**. `imprint all` + `merge all` + hex sweep + block — no build123d equivalents. Keep `.jou` for hex → Radia/ELF. |

## Scoring rubric

For each fixture:

1. **Verdict correctness** — did the LLM translate vs skip in line with
   the expected column above?
2. **Python validity** — does the output import and parse without
   syntax errors?
3. **Geometric correctness** — does running the output under
   `build123d` produce a valid Shape with a sensible bounding box and
   volume?
4. **Idiomaticity** — does the output use build123d algebra or builder
   patterns rather than one-off hacks? Does it cite the source `.jou`?

## Output

The benchmark writer produces:

- `out/NN_<fixture>.py` — generated build123d, or a stub with the skip
  reason if translation was refused. These generated outputs are not canonical
  examples; the 2026-06-27 generated set was pruned from `examples/`.
- `out/report.md` — per-fixture rubric scoring + list of mcp-server
  knowledge gaps observed.

## Running

This benchmark is normally driven by the parent Claude Code session
spawning a sub-agent with Read/Write access and the knowledge base
files available. The sub-agent is told explicitly to use the lab
knowledge, not generic training-data knowledge, and to flag the gaps.

The `.jou` fixtures stay here. Generated Python outputs should be regenerated
outside the tracked examples tree unless they are promoted into a tested
translator fixture.
