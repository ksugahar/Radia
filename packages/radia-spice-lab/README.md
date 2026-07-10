# radia-spice-lab

[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](https://github.com/ksugahar/Radia/tree/main/packages/radia-spice-lab)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`radia-spice-lab` is the Radia monorepo package for circuit-aware
SPICE/LTspice work. It absorbs the former standalone `spice-circuit-lab`
repository while keeping the legacy `spice-circuit-lab`,
`ltspice-convert`, and `mcp-ltspice` entry points available.

It converts `.asc`, `.cir`, and schemdraw Python, checks whether round-trips
preserve topology, and exposes MCP tools so AI agents can author, inspect,
and refine circuits without losing the electrical intent.

## What this library is

`radia-spice-lab` is a public LTspice interoperability layer for Radia, not a
dump of private simulation corpora.  It provides deterministic, testable tools
for moving a circuit between the formats engineers actually use:

- LTspice `.asc` schematics.
- SPICE `.cir` / `.net` netlists.
- Runnable schemdraw Python scripts for publication figures and review.
- Topology signatures that catch silent rewiring after conversion.
- Public-safe LTspice `.measure` log summaries for scalar simulation evidence.
- MCP tools that let an AI agent convert, inspect, lint, and reason about a
  circuit without copying private RAW files.

It does **not** redistribute LTspice, LTspice's bundled examples, textbook
circuits, third-party circuit corpora, private RAW files, or the LAB-private
invention-learning pipeline.  Those materials may be used locally to improve
the converter, but only the authored converter, small fixtures, and public-safe
summary schemas belong in this repository.

The public import names remain:

- `import spice_circuit_lab`
- `import ltspice_converter`

The package is also exposed through the Radia monorepo as
`packages/radia-spice-lab`.

The lineage is `ltspice-converter` -> `spice-circuit-lab` ->
`radia-spice-lab`; the older package and command names remain available for
backward compatibility.

Core capabilities:

- Convert between LTspice `.asc`, SPICE `.cir`, and runnable schemdraw scripts.
- Prefer LTspice's own `-netlist` backend when available, with deterministic
  pure-Python extraction when requested.
- Detect topology drift, not just component-count drift, so silent rewiring
  does not pass as a clean conversion.
- Preserve difficult SPICE forms such as controlled sources, behavioral
  expressions, switch models, and inline subcircuit parameters.
- Summarize LTspice `.measure` logs into a public-safe JSON schema without
  exposing RAW waveforms or private simulation directories.
- Provide public circuit-knowledge helpers and MCP tools for agentic circuit
  design workflows.

Conversion graph:

```
   LTspice .asc  <---->  SPICE .cir  <---->  schemdraw Python script
```

What's new in [v0.4.0](CHANGELOG.md): when **LTspice is installed it is
now the default `.asc -> netlist` backend** (its own `-netlist` is the
ground truth; pure-Python remains the deterministic opt-out via
`--no-ltspice`). Auditing pure-Python against LTspice as an oracle fixed
real extraction bugs around jumper net-ties, inline subcircuit parameters,
special-function `A` devices, `.asy` `Prefix` classes, and BOM-less UTF-16
input. v0.3.14 added a **node-rename-invariant topology check** that catches
silent rewiring even when component counts still match. See
[docs/BENCHMARKS.md](docs/BENCHMARKS.md) for corpus results.

Works without LTspice (pure-Python), but **uses LTspice's own
`-netlist` automatically when LTspice.exe is installed** — that is the
ground truth for its own `.asc` format, so `.asc → netlist` extraction
is canonical (correct vendor-symbol topology, jumpers and special
functions handled exactly as LTspice does). Set `use_ltspice=False`
(or `--no-ltspice`) to force the deterministic pure-Python path. Built
for AI agents to round-trip circuits between schematic and netlist
forms.

## Install

The package lives in the Radia monorepo. Install directly from the
monorepo subdirectory:

```bash
pip install "git+https://github.com/ksugahar/Radia.git#subdirectory=packages/radia-spice-lab"
```

For MCP server support (Claude Code, Cursor, etc.):

```bash
pip install "radia-spice-lab[mcp] @ git+https://github.com/ksugahar/Radia.git#subdirectory=packages/radia-spice-lab"
```

Pinning a specific Radia tag or commit:

```bash
pip install "git+https://github.com/ksugahar/Radia.git@v4.95.9#subdirectory=packages/radia-spice-lab"
```

For development:

```bash
git clone https://github.com/ksugahar/Radia
cd Radia/packages/radia-spice-lab
pip install -e .[test,mcp]
```

## Python API

```python
import spice_circuit_lab as scl

netlist = """* RC Lowpass Filter
V1 in 0 AC 1
R1 in in1 1k
C1 in1 0 1u
.ac dec 20 1 100k
.end"""

# SPICE netlist -> runnable schemdraw script
script = scl.netlist_to_schemdraw(netlist, name="rc")

# schemdraw script -> SPICE netlist
recovered = scl.schemdraw_to_netlist(script, title="rc")

# SPICE netlist -> LTspice .asc text
asc_text = scl.netlist_to_asc(netlist)

# LTspice .asc text -> SPICE netlist. use_ltspice=None (default) auto-uses
# LTspice.exe when installed (canonical, ground-truth topology) and falls
# back to pure-Python otherwise. Force with use_ltspice=True / False.
recovered_netlist = scl.asc_to_netlist(asc_text)
deterministic   = scl.asc_to_netlist(asc_text, use_ltspice=False)
```

Legacy imports still work:

```python
import ltspice_converter as lc
```

## Circuit-knowledge helpers

The package includes small public design helpers for simulation seeds.  These
are engineering rules of thumb, not sign-off designs.

```python
import spice_circuit_lab as scl

seed = scl.buck_seed(24, 5, 1, fsw_hz=100_000)
print(seed.to_dict())
print(seed.to_netlist())

rules = scl.circuit_knowledge("buck converter")
for rule in rules["rules"]:
    print("-", rule)

plan = scl.patentability_search_plan(
    title="snubber-assisted boost converter",
    features=["boost converter", "switch-node RC snubber", "soft start"],
    effects=["reduced ringing", "lower overshoot"],
    domains=["power electronics", "circuit"],
)
print(plan["query_sets"]["google_patents"])
print(plan["jplatpat_keywords_ja"])
```

## Command-line tool

Installing the package wires up the `spice-circuit-lab` console script.
The old `ltspice-convert` command remains available.
No Python knowledge needed.

### Conversion

```bash
# Single file (target inferred from -o or the "opposite" of input)
spice-circuit-lab input.asc -o output.cir
spice-circuit-lab input.cir -o output.asc
spice-circuit-lab input.cir -o output.py        # netlist -> schemdraw script

# Auto output path (same dir, sensible default extension)
spice-circuit-lab input.asc                     # writes input.cir alongside
spice-circuit-lab input.asc --to py             # writes input.py alongside

# Batch (output is a directory; --to picks target format)
spice-circuit-lab *.asc -o build/ --to cir

# Backend for .asc -> netlist: default auto-uses LTspice when installed.
spice-circuit-lab input.asc -o output.cir --no-ltspice   # force pure-Python
spice-circuit-lab input.asc -o output.cir --use-ltspice  # force LTspice
```

### Round-trip check (lint mode)

`--check` reads a file, runs it through the full round-trip, and
reports drift / `.asy` resolution gaps. Use `--strict` to make any
warning exit non-zero -- handy in CI.

```bash
spice-circuit-lab --check input.asc
# -> PASS / PASS (with warnings) / FAIL  on stdout, with [ok]/[warn] details

spice-circuit-lab --check --strict *.asc        # exit 1 if any warning
```

The round-trip arm reports three drift signals: **component count**,
**GND-pin position**, and — as of v0.3.14 — **topology**
(node-rename-invariant connectivity). The topology line is the one that
catches a multi-pin vendor symbol whose `.asy` is missing: the count
stays right but the wiring changes.

Static netlist checks `--check` runs (in addition to round-trip):

- Duplicate instance names (`R1` appearing twice at top level)
- Floating nodes (only one component touches it -- usually a wire
  the user forgot to finish)
- Orphan `.model` declarations (model defined but never referenced)
- Undefined model references (device names a model that is not
  defined inline; standard library names like `1N4148`, `2N3904`,
  `LT1001`, ... are exempt)
- `{PARAM}` references without a matching `.param NAME=...`
- Lines the parser could not classify (with a "did you mean ...?"
  hint for common typos like `Resistor` / `Capacitor`)

### Info / stats

```bash
spice-circuit-lab --info input.asc              # human-readable
spice-circuit-lab --info --json input.asc       # machine-readable JSON
```

### Third-party `.asy` libraries

`--asy-dir` (repeatable) is equivalent to setting the
`LTSPICE_ASY_SEARCH_PATH` env var:

```bash
spice-circuit-lab --asy-dir /path/to/MyLib/sym input.asc -o out.cir
spice-circuit-lab --asy-dir A --asy-dir B input.asc -o out.cir
```

CLI flags take priority over the env var; both can be combined.

### CI / GitHub Actions

Run `--check --strict` on every PR that touches `.asc` files. See
[`docs/example-workflows/asc-check.yml`](docs/example-workflows/asc-check.yml)
for a reusable workflow template you can copy into a circuit-design
repository.

## MCP server

For AI agents (Claude Code, Cursor) that author or refactor SPICE:
the MCP server lets the agent convert between formats AND lint its
own generated netlists in the same conversation, without shelling out
to a CLI.

Install with the `[mcp]` extra and add to your MCP client config:

```json
{
  "mcpServers": {
    "ltspice": {
      "command": "mcp-spice-circuit-lab"
    }
  }
}
```

Exposes eleven tools:

| Tool | Purpose |
|---|---|
| `netlist_to_schemdraw(netlist, name)` | SPICE → schemdraw Python script |
| `schemdraw_to_netlist(script, title)` | schemdraw script → SPICE |
| `netlist_to_asc(netlist, asy_search_dirs?)` | SPICE → LTspice `.asc` |
| `asc_to_netlist(asc_text, use_ltspice?, asy_search_dirs?)` | LTspice `.asc` → SPICE (`use_ltspice` defaults to auto: LTspice if installed, else pure-Python) |
| `check_circuit(text, fmt, asy_search_dirs?)` | Lint: round-trip drift (count, GND-pin, **topology**) + static netlist checks. Returns `{ok, info, warnings}`. |
| `info_circuit(text, fmt, asy_search_dirs?)` | Summary: component counts, symbol kinds, `.subckt` blocks. |
| `compare_topology(netlist_a, netlist_b)` | Node-rename-invariant connectivity diff of two netlists. Returns `{equivalent, ...}`. |
| `parse_measure_log(log_text)` | Public-safe LTspice `.measure` log parser. Returns schema `radia-spice-lab.measure-log.v1`, measure rows, `.step` rows, and duplicate-name warnings. |
| `circuit_knowledge(topic)` | Compact public circuit-design and conversion rules by topic. |
| `buck_seed(vin_v, vout_v, iout_a, fsw_hz?, ripple_fraction?)` | First-pass asynchronous buck sizing plus an LTspice-ready open-loop netlist. |
| `patentability_search_plan(title, features, effects?, domains?, include_japanese?)` | Non-legal prior-art search plan for Google Scholar, Google Patents, J-PlatPat, and web searches. |

Typical agent loop: generate netlist → `check_circuit(..., 'cir')` →
if `warnings` non-empty, fix and re-check → only ship when clean.

`patentability_search_plan` is only a search-query and report-planning aid.
It does not decide patentability and is not a legal opinion.

`compare_topology` answers a different question — *"did my edit change
the wiring?"*  It is invariant to node renaming and benign R/C/L pin
swaps, so after changing a component value
`compare_topology(before, after)` returns `equivalent: true`; if you
accidentally moved a wire it returns `false`.  Use it to confirm an
edit touched only what you intended.

## End-to-end workflow: AI-assisted circuit editing

A typical loop that exercises every layer of the converter — from
the LTspice schematic, through an AI agent, and back to LTspice —
looks like this:

```
                    +----------+    +-----+    +-----+
   LTspice GUI ---> | foo.asc  |--->|.cir |    |.py  | --> PDF/SVG
   (human draws)    +----------+    +-----+    +-----+    (publishable)
                         ^             |
                         |             v
                         |       Claude / Cursor
                         |       (edits .cir)
                         |             |
                         |             v
                         |       check_circuit(...)
                         |             |
                         +-------------+
                              (ships only when warnings = [])
```

### Step-by-step

1. **Start in LTspice**: open or draw a schematic, save as `foo.asc`.

2. **Hand to an AI agent**.  In Claude Code or Cursor with
   `mcp-spice-circuit-lab` configured (`mcp-ltspice` also works), ask the agent to modify the circuit.
   The agent calls:

   ```
   asc_to_netlist(asc_text=<contents of foo.asc>)
   → returns the SPICE netlist text
   ```

3. **Agent edits the netlist** in conversation (add a resistor,
   change a model, swap an op-amp).  Before claiming the job done,
   it validates:

   ```
   check_circuit(text=<edited netlist>, fmt='cir')
   → {"ok": false, "warnings": ["floating node N003"]}
   ```

   If `warnings` is non-empty the agent fixes and re-checks.  It
   only "ships" the result when `ok == true`.

4. **Back to `.asc`** so the user can verify in LTspice:

   ```
   netlist_to_asc(netlist=<clean netlist>)
   → returns .asc text; user saves as foo_v2.asc
   ```

5. **Reopen in LTspice** to visually inspect and run the simulation.
   The regenerated schematic looks like one a human would have drawn
   — `.asc → .cir → .asc` count match is 100 % on real-world corpora
   (see [docs/BENCHMARKS.md](docs/BENCHMARKS.md)).

### Same loop without MCP (just the CLI)

```bash
# 1-2. Extract netlist from LTspice schematic
spice-circuit-lab foo.asc -o foo.cir

# 3. Edit foo.cir by hand (or with any tool), then lint
spice-circuit-lab --check --strict foo.cir

# 4-5. Regenerate the schematic for LTspice
spice-circuit-lab foo.cir -o foo_v2.asc

# Or render to PDF/SVG via schemdraw — no LTspice install needed,
# useful for papers, slides, and web pages:
spice-circuit-lab foo.cir -o foo.py && python foo.py
```

### Why this matters

LTspice's `.asc` is a custom format that does not diff cleanly and is
not portable outside Windows/Mac.  By going through the SPICE netlist
(plain text, diff-able, standard-format) and optionally a schemdraw
Python script (runnable, publishable, AI-readable), this converter
lets you:

- **review circuit changes in git** like any other source file,
- **let an AI agent author or refactor circuits** with a verification
  loop that catches drift,
- **render publication-quality figures** without launching LTspice,
- **work on Linux** where LTspice is harder to install.

The whole point of the v0.3.8 - v0.3.13 work was making this loop
trustworthy enough that the agent's "ship" decision can be taken at
face value.

## Supported elements

| SPICE | schemdraw | LTspice symbol |
|-------|-----------|----------------|
| R     | Resistor  | res, res2      |
| C     | Capacitor | cap, polcap    |
| L     | Inductor  | ind, ind2      |
| V     | SourceV   | voltage        |
| I     | SourceI   | current        |
| D     | Diode     | diode, zener, schottky, varactor, tvs |
| Q (NPN/PNP, 3- or 4-pin substrate) | BjtNpn / BjtPnp | npn, pnp, npn3, pnp3, npn4, pnp4 |
| M (NMOS/PMOS, 3- or 4-pin substrate) | NFet / PFet | nmos, pmos, nmos4, pmos4 |
| J (NJF/PJF, 3- or 4-pin substrate) | JFetN / JFetP | njf, pjf, njf4, pjf4 |
| B (behavioral source) | — | bv, bi, bi2 |
| E, G (VCVS, VCCS) | — | e, e2, g, g2 |
| F, H (CCCS, CCVS) | — | f, h |
| S (voltage-controlled switch) | — | sw |
| T (transmission line) | — | tline |
| K (mutual inductance, directive) | — | — |
| X (subcircuit / opamp / IC) | Opamp | opamp, opamp2, lt1018, ... and arbitrary multi-pin vendor symbols |
| U (digital flop / gate) | — | Digital\\\\srflop, Comparators\\\\..., ... |

`.ac`, `.tran`, `.op`, `.dc`, `.param`, `.model`, `.subckt`/`.ends`
directives are preserved through the round-trip.

## `.subckt` round-trip

The converter preserves an entire `.subckt` block --- header, body
components, models, `.param`s, comments, and `.ends` line --- byte
for byte. This means a netlist file like:

```spice
* myckt with a diac model
V1 in 0 SINE(0 230 50)
X1 in 0 mydiac
.tran 20m

.subckt mydiac T1 T2
* simplified DIAC: two opposing zeners
.model BD D Bv=30
D1 T1 T2 BD
D2 T2 T1 BD
.ends mydiac
.end
```

round-trips through `spice-circuit-lab` cleanly:

```bash
spice-circuit-lab myckt.cir -o myckt.asc      # write LTspice schematic
spice-circuit-lab myckt.asc -o back.cir       # extract back
diff myckt.cir back.cir                     # node names may rename;
                                            # the .subckt block is byte-equal
```

The same holds for `spice-circuit-lab --check myckt.cir`: any drift
inside the subckt body would show up as `component count drift` or a
parser warning.

## Performance

Current `.asc → netlist → .asc` round-trip pass rate on real-world
corpora (component count preserved):

| Source | Pass rate |
|---|---|
| LTspice "Applications" examples (4099 files) | 100% |
| LTspice "Examples" examples (100 samples)    | 100% |
| Textbook circuits (training-adjacent, 28 files) | 100% |
| GitHub repos (unseen, 720 files)              | 100% |

The schemdraw round-trip `netlist → schemdraw script → netlist`
(component count match) runs at **100%** on the same corpora after
the v0.3.11 K-directive / multi-line `.subckt` fixes and the v0.3.13
transmission-line (`T`) fix.

### Third-party symbol libraries

For round-tripping schematics that use third-party LTspice libraries
(not bundled in `lib.zip`), point the
`LTSPICE_ASY_SEARCH_PATH` environment variable at the library's
`sym/` root directory. Multiple paths are separated by the OS path
separator (`;` on Windows, `:` on Linux/macOS):

```bash
# Linux / macOS
export LTSPICE_ASY_SEARCH_PATH="/path/to/LTspiceControlLibrary/lib/sym:/path/to/MyLib/sym"

# Windows (cmd.exe)
set LTSPICE_ASY_SEARCH_PATH=C:\Libs\LTspiceControlLibrary\lib\sym;C:\Libs\MyLib\sym
```

When the env var is set, both `asc → netlist` and `netlist → asc`
use the same `.asy` files, so node-to-pin topology survives the
round-trip for those vendor symbols.

See [docs/BENCHMARKS.md](docs/BENCHMARKS.md) for methodology, the
schemdraw arm, and known failure modes.  Pass rate is not perfect —
file a [GitHub issue](https://github.com/ksugahar/Radia/issues)
with a failing `.asc` and we'll fix it.

## Project history & scope

This package is the successor to the (now-private)
`ksugahar/circuit-converter` repository.  The previous repo also contained
training/research material — textbook circuits, LTspice's own bundled
examples (Educational/, Applications/), and aggregated third-party GitHub
content — none of which can be redistributed.

This re-release is intentionally **converter-only**: just the trained
converter, a small set of author-authored test circuits, and API
documentation.  Higher-level circuit analysis features that lived in
the old repo are not in scope here.

### Distribution policy: trained converter, private corpus

Development uses **private training corpora** (textbook circuits,
LTspice bundled examples, third-party GitHub aggregates) to measure
round-trip fidelity and to surface converter bugs.  **Only the
converter is published** — every corpus stays on the author's
LAB-private storage.  Concretely:

- Public: source code, CLI, MCP server, API docs, the LAB-authored
  test fixtures under `tests/fixtures/` (~58 minimal RC/RLC/filter
  circuits the author wrote from scratch).
- LAB-private: textbook DB JSONs, the ~5,000-file `bench/` corpus,
  the harvest pipeline scripts, the textbook PDFs themselves.

The benchmark numbers in [docs/BENCHMARKS.md](docs/BENCHMARKS.md)
report measurements made on LAB-private corpora; you cannot reproduce
them bit-for-bit without your own corpus, but the converter behaviour
they describe is exactly what you get from installing
`packages/radia-spice-lab` from the Radia monorepo.

## Test fixtures

Everything under `tests/fixtures/` is author-authored: small RC / RLC /
filter circuits in either of three forms (`.asc`, `.cir`, `.gen.py`)
used to exercise the round-trip.  No textbook content, no LTspice
bundled examples, no third-party repo dumps.

To test against larger sets, run the converter on **your own local copy**
of LTspice's `Educational/` and `Applications/` directories — do not
redistribute the results.

## Copyright notice

| Asset | Owner | Status here |
|-------|-------|-------------|
| Converter source code | © 2026 Mitsutoshi Sugahara | **MIT** |
| Test fixtures in `tests/fixtures/` | © 2026 Mitsutoshi Sugahara | **MIT** |
| API reference in `docs/pyltspice_api.md` | reformat of public PyLTSpice/spicelib API | fair-use technical citation |
| LTspice itself | © Analog Devices, Inc. | not redistributed (install separately) |
| LTspice bundled example circuits | © Analog Devices, Inc. | **excluded** |
| Textbook circuits / problem sets used during development | © respective publishers | **excluded** |
| Third-party GitHub circuits used during development | © respective authors | **excluded** |

## Author

[Mitsutoshi Sugahara (菅原光俊)](https://github.com/ksugahar) —
Department of Electric and Electronic Engineering, Kindai University

## License

MIT — see [LICENSE](LICENSE).

