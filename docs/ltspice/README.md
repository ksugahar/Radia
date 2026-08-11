# SPICE and LTspice integration

`radia.ltspice` is Radia's built-in circuit interoperability layer. It converts
`.asc`, `.cir`, and schemdraw Python, verifies that round-trips preserve
electrical topology, parses public-safe simulation summaries, and exposes the
same operations through Python, MATLAB, CLI, and MCP.

This belongs in Radia because circuit extraction is part of the electromagnetic
engineering loop: PEEC and reduced electromagnetic models become circuit
elements, circuits provide excitation and loading, and Simulink coordinates the
coupled system. There is one Python implementation under `src/radia/ltspice`.

## Capabilities

The package provides deterministic, testable movement between the formats
engineers use:

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

The canonical import is `radia.ltspice`. The former `ltspice_converter` and
`spice_circuit_lab` namespaces remain thin compatibility aliases and contain no
second implementation.

Existing `radia-spice-lab.*.v1` JSON schema identifiers are intentionally
retained so saved circuit evidence remains readable across the integration.

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

Release history now follows Radia's unified [changelog](../../CHANGELOG.md).
See [BENCHMARKS.md](BENCHMARKS.md) for the public evaluation methodology.

Works without LTspice (pure-Python), but **uses LTspice's own
`-netlist` automatically when LTspice.exe is installed** — that is the
ground truth for its own `.asc` format, so `.asc → netlist` extraction
is canonical (correct vendor-symbol topology, jumpers and special
functions handled exactly as LTspice does). Set `use_ltspice=False`
(or `--no-ltspice`) to force the deterministic pure-Python path. Built
for AI agents to round-trip circuits between schematic and netlist
forms.

## Install

Install Radia with schemdraw support:

```bash
python -m pip install "radia[ltspice]"
```

For development:

```bash
git clone https://github.com/ksugahar/Radia
cd Radia
python -m pip install -e ".[ltspice]"
```

## Python API

```python
from radia import ltspice

netlist = """* RC Lowpass Filter
V1 in 0 AC 1
R1 in in1 1k
C1 in1 0 1u
.ac dec 20 1 100k
.end"""

# SPICE netlist -> runnable schemdraw script
script = ltspice.netlist_to_schemdraw(netlist, name="rc")

# schemdraw script -> SPICE netlist
recovered = ltspice.schemdraw_to_netlist(script, title="rc")

# SPICE netlist -> LTspice .asc text
asc_text = ltspice.netlist_to_asc(netlist)

# LTspice .asc text -> SPICE netlist. use_ltspice=None (default) auto-uses
# LTspice.exe when installed (canonical, ground-truth topology) and falls
# back to pure-Python otherwise. Force with use_ltspice=True / False.
recovered_netlist = ltspice.asc_to_netlist(asc_text)
deterministic = ltspice.asc_to_netlist(asc_text, use_ltspice=False)
```

Legacy imports remain available during the compatibility window:

```python
import ltspice_converter  # Prefer radia.ltspice in new code.
```

## Circuit-knowledge helpers

The package includes small public design helpers for simulation seeds.  These
are engineering rules of thumb, not sign-off designs.

```python
from radia import ltspice

seed = ltspice.buck_seed(24, 5, 1, fsw_hz=100_000)
print(seed.to_dict())
print(seed.to_netlist())

rules = ltspice.circuit_knowledge("buck converter")
for rule in rules["rules"]:
    print("-", rule)

plan = ltspice.patentability_search_plan(
    title="snubber-assisted boost converter",
    features=["boost converter", "switch-node RC snubber", "soft start"],
    effects=["reduced ringing", "lower overshoot"],
    domains=["power electronics", "circuit"],
)
print(plan["query_sets"]["google_patents"])
print(plan["jplatpat_keywords_ja"])
```

## Command-line tool

Installing the package wires up the `radia-ltspice` console script.
The old `ltspice-convert` command remains available.
No Python knowledge needed.

### Conversion

```bash
# Single file (target inferred from -o or the "opposite" of input)
radia-ltspice input.asc -o output.cir
radia-ltspice input.cir -o output.asc
radia-ltspice input.cir -o output.py        # netlist -> schemdraw script

# Auto output path (same dir, sensible default extension)
radia-ltspice input.asc                     # writes input.cir alongside
radia-ltspice input.asc --to py             # writes input.py alongside

# Batch (output is a directory; --to picks target format)
radia-ltspice *.asc -o build/ --to cir

# Backend for .asc -> netlist: default auto-uses LTspice when installed.
radia-ltspice input.asc -o output.cir --no-ltspice   # force pure-Python
radia-ltspice input.asc -o output.cir --use-ltspice  # force LTspice
```

### Round-trip check (lint mode)

`--check` reads a file, runs it through the full round-trip, and
reports drift / `.asy` resolution gaps. Use `--strict` to make any
warning exit non-zero -- handy in CI.

```bash
radia-ltspice --check input.asc
# -> PASS / PASS (with warnings) / FAIL  on stdout, with [ok]/[warn] details

radia-ltspice --check --strict *.asc        # exit 1 if any warning
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
radia-ltspice --info input.asc              # human-readable
radia-ltspice --info --json input.asc       # machine-readable JSON
```

### Third-party `.asy` libraries

`--asy-dir` (repeatable) is equivalent to setting the
`LTSPICE_ASY_SEARCH_PATH` env var:

```bash
radia-ltspice --asy-dir /path/to/MyLib/sym input.asc -o out.cir
radia-ltspice --asy-dir A --asy-dir B input.asc -o out.cir
```

CLI flags take priority over the env var; both can be combined.

### CI / GitHub Actions

Run `--check --strict` on every PR that touches `.asc` files. See
[`example-workflows/asc-check.yml`](example-workflows/asc-check.yml)
for a reusable workflow template you can copy into a circuit-design
repository.

## MCP server

For AI agents (Claude Code, Cursor) that author or refactor SPICE:
the MCP server lets the agent convert between formats AND lint its
own generated netlists in the same conversation, without shelling out
to a CLI.

The MCP SDK ships with Radia. Add the server to your MCP client config:

```json
{
  "mcpServers": {
    "ltspice": {
      "command": "mcp-server-radia-ltspice"
    }
  }
}
```

Exposes thirty-two tools:

| Tool | Purpose |
|---|---|
| `netlist_to_schemdraw(netlist, name)` | SPICE → schemdraw Python script |
| `schemdraw_to_netlist(script, title)` | schemdraw script → SPICE |
| `netlist_to_asc(netlist, asy_search_dirs?)` | SPICE → LTspice `.asc` |
| `asc_to_netlist(asc_text, use_ltspice?, asy_search_dirs?)` | LTspice `.asc` → SPICE (`use_ltspice` defaults to auto: LTspice if installed, else pure-Python) |
| `check_circuit(text, fmt, asy_search_dirs?)` | Lint: round-trip drift (count, GND-pin, **topology**) + static netlist checks. Returns `{ok, info, warnings}`. |
| `info_circuit(text, fmt, asy_search_dirs?)` | Summary: component counts, symbol kinds, `.subckt` blocks. |
| `compare_topology(netlist_a, netlist_b)` | Node-rename-invariant connectivity diff of two netlists. Returns `{equivalent, ...}`. |
| `balanced_learning_profile()` | Equal-capability public/source MCP learning contract and self-check. |
| `parse_measure_log(log_text)` | Public-safe LTspice `.measure` log parser. It recovers linear `mag()` values from the AC dB wrapper and rejects scalar `ph()` evidence whose sign was lost. |
| `parse_stepped_measure_log(log_text)` | Pairs stepped `Measurement` rows with concrete `.step` assignments and rejects incomplete tables. |
| `sallen_key_filter_family_gate(rows)` | Gates multi-Q low-pass sweeps against the ideal two-pole response. |
| `hysteretic_inductor_cycle_gate(cycle_rows, expected_current_peak_a, expected_copper_energy_j, voltage_thd)` | Gates settled hysteresis cycles by terminal/copper/flux-loop energy, closure, repeatability, and harmonic evidence. |
| `half_wave_rectifier_gate(...)` | Gates capacitor-input half-wave DC, ripple, and load/diode current balance. |
| `bridge_rectifier_gate(summary)` | Gates single-phase full-wave frequency doubling, diagonal-pair conduction, four-diode/load current balance, and output-node KCL. |
| `cockcroft_walton_stage_gate(...)` | Gates loaded two-stage multiplier scaling, adjacent-window settling, ripple, load current, and real-power bounds. |
| `boost_converter_steady_state_gate(...)` | Gates periodic boost ratio, passive power, inductor volt-second balance, and capacitor charge balance. |
| `transient_psrr_gate(...)` | Gates transient input/output ripple attenuation against independent RAW replay and the PSRR definition. |
| `measure_bandwidth_crossing_gate(summary)` | Replays two-sided rise/fall `-3 dB` bandwidth measures from sampled AC magnitude, using dB interpolation on each linear-frequency interval. |
| `bipolar_supply_startup_gate(summary)` | Gates signed dual-rail startup by final regulation, rail balance, tail ripple, overshoot, 10-90% timing skew, and power-good ordering. |
| `bipolar_converter_efficiency_gate(summary)` | Gates signed source power, positive dual-output power, passivity, output balance, and recomputation of reported late-window efficiency. |
| `bipolar_rail_power_quality_gate(summary)` | Gates signed bipolar rails, target regulation, ripple, balance, passivity, and late-window efficiency closure. |
| `monte_carlo_tolerance_family_gate(summary)` | Gates independent uniform component-tolerance statistics against the expected standard deviation and `1/sqrt(N)` averaging of equivalent resistors and symmetric dividers. |
| `series_rlc_complex_impedance_gate(summary)` | Gates a current-driven series RLC using source orientation, three full-complex voltage traces, analytic resonance, minimum impedance, generated-netlist semantics, and converted/reference RAW equivalence. |
| `rc_thermal_noise_psd_gate(summary)` | Gates RC thermal-noise density, finite-band RMS from squared-PSD integration, capacitance scaling, units, and noise-analysis-specific `INTEG` semantics. |
| `distributed_line_delay_loss_gate(summary)` | Separates distributed-line `length*sqrt(LC)` first-arrival delay from series-resistance attenuation, catches the SPICE `m` = milli suffix, and requires observable replay. |
| `second_order_allpass_phase_group_delay_gate(summary)` | Gates flat all-pass magnitude together with the full complex transfer, monotone phase winding, reciprocal-frequency phase symmetry, center group delay `4Q/omega0`, pole-zero mirroring, and replay. |
| `second_order_complex_zero_transfer_gate(summary)` | Gates a two-pole/two-complex-zero transfer using root-recovered frequencies and Q, DC/high-frequency gain ratio, finite real-axis dip, full complex response, and replay. |
| `balanced_three_phase_delta_load_gate(summary)` | Gates a balanced Y source feeding an equal delta-connected RL load using ABC sequence, sqrt(3) voltage/current identities, branch impedance, complex-power conservation, constant instantaneous three-phase power, and replay. |
| `ideal_transformer_identity_gate(summary)` | Gates a two-winding ideal transformer using signed voltage/current ratios, reflected load, source-resistance KVL, complex and instantaneous power conservation, and replay. |
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
   `mcp-server-radia-ltspice` configured (`mcp-ltspice` also works), ask the agent to modify the circuit.
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
   (see [BENCHMARKS.md](BENCHMARKS.md)).

### Same loop without MCP (just the CLI)

```bash
# 1-2. Extract netlist from LTspice schematic
radia-ltspice foo.asc -o foo.cir

# 3. Edit foo.cir by hand (or with any tool), then lint
radia-ltspice --check --strict foo.cir

# 4-5. Regenerate the schematic for LTspice
radia-ltspice foo.cir -o foo_v2.asc

# Or render to PDF/SVG via schemdraw — no LTspice install needed,
# useful for papers, slides, and web pages:
radia-ltspice foo.cir -o foo.py && python foo.py
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

round-trips through `radia-ltspice` cleanly:

```bash
radia-ltspice myckt.cir -o myckt.asc      # write LTspice schematic
radia-ltspice myckt.asc -o back.cir       # extract back
diff myckt.cir back.cir                     # node names may rename;
                                            # the .subckt block is byte-equal
```

The same holds for `radia-ltspice --check myckt.cir`: any drift
inside the subckt body would show up as `component count drift` or a
parser warning.

## Quality gates

Radia checks more than component counts. The public regression lane verifies
node-rename-invariant topology, grounded-pin placement, subcircuit retention,
controlled and behavioral sources, generated script execution, CLI error
behavior, and measure-log schemas. See [BENCHMARKS.md](BENCHMARKS.md) for the
reproducible public methodology.

### Third-party symbol libraries

For round-tripping schematics that use third-party LTspice libraries
(not bundled in `lib.zip`), point the
`LTSPICE_ASY_SEARCH_PATH` environment variable at the library's
`sym/` root directory. Multiple paths are separated by the OS path
separator (`;` on Windows, `:` on Linux/macOS):

```bash
# Linux / macOS
export LTSPICE_ASY_SEARCH_PATH="/path/to/library-a/sym:/path/to/library-b/sym"

# Windows (cmd.exe)
set LTSPICE_ASY_SEARCH_PATH=C:\path\to\library-a\sym;C:\path\to\library-b\sym
```

When the env var is set, both `asc → netlist` and `netlist → asc`
use the same `.asy` files, so node-to-pin topology survives the
round-trip for those vendor symbols.

If a public, redistributable circuit fails, attach the smallest reproducer to a
[Radia issue](https://github.com/ksugahar/Radia/issues).

## Test fixtures

Everything under `tests/ltspice/fixtures/` is author-authored: small RC / RLC /
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
| Test fixtures in `tests/ltspice/fixtures/` | © 2026 Mitsutoshi Sugahara | **MIT** |
| API reference in `docs/ltspice/pyltspice_api.md` | reformat of public PyLTSpice/spicelib API | fair-use technical citation |
| LTspice itself | © Analog Devices, Inc. | not redistributed (install separately) |
| LTspice bundled example circuits | © Analog Devices, Inc. | **excluded** |
| Textbook circuits / problem sets used during development | © respective publishers | **excluded** |
| Third-party GitHub circuits used during development | © respective authors | **excluded** |

## Author

[Mitsutoshi Sugahara (菅原光俊)](https://github.com/ksugahar) —
Department of Electric and Electronic Engineering, Kindai University

## License

MIT — see [component license](../../src/radia/ltspice/LICENSE).

