# Benchmarks

This document reports the converter's performance on a corpus of
real-world LTspice circuits.  All numbers come from LAB-local
benchmarks; the corpus itself is not redistributed (see [README §
Project history & scope](../README.md)).

## LTspice as ground-truth oracle — v0.4.0 (current)

When `LTspice.exe` is installed the converter now extracts via LTspice's
own `-netlist` by default. That also makes LTspice a **ground-truth
oracle** for auditing the pure-Python fallback: for each `.asc`, compare
`asc_to_netlist(use_ltspice=False)` against `asc_to_netlist(use_ltspice=True)`
under the node-rename-invariant topology signature.

On the bundled "Examples" corpus (100 files) the pure-Python extractor
matched LTspice's canonical topology only **55 %** at the start. Four
oracle-surfaced bugs were fixed:

| Fix | What LTspice does that pure-Python didn't |
|---|---|
| jumper net-tie | merges the two pins into one node, drops the symbol |
| inline subckt params | `X1 .. sub p=1 q=2` — `p=1`/`q=2` are params, not pins |
| special-function `A` devices | emits `A<n> ..` (S&H, varistor, PLL, Schmitt, ...) |
| `.asy` `SYMATTR Prefix` | classes a symbol by its own declared SPICE prefix (e.g. `xtal` → `C`) |

Result on Examples: strict topology match **55 % → 71 %**; counting the
op-amp `U`-vs-`X` class label (same wiring) as equivalent, **88 %** of
circuits have provably-correct wiring. The residual is the documented
multi-pin-vendor-symbol case (needs the symbol's `.asy` pin geometry) —
which the LTspice backend now handles for free. Every fix was gated on
zero self-round-trip regression across the full 4099-file Applications
corpus.

Reproduce on your own `.asc` corpus with the LAB harness
`bench/baseline.py` (gitignored, not shipped — point its `CORPUS_ROOTS`
at your files). Its `roundtrip_topology_match` column compares the
pure-Python self-round-trip; the oracle audit lives in the same harness.

## Topology (connectivity) preservation — v0.3.14

Component **count** and **GND-pin position** are necessary but not
sufficient round-trip metrics: a circuit can keep every component and
every ground connection yet have its internal wiring scrambled.
v0.3.14 adds a **node-rename-invariant topology signature**
(`ltspice_converter.topology`, Weisfeiler-Leman 1-WL on the
component↔node incidence graph) and measures it across the round-trip:

| Corpus | files | count match | **topology match** |
|---|---:|---:|---:|
| textbook              |  110 | 100 %  | **100.0 %** |
| LTspice Applications  | 4099 | 100 %  | **96.9 %** |
| LTspice Examples      |  100 | 99 %   | **92.0 %** |
| GitHub repos (unseen) |  720 | 100 %  | **79.3 %** |

The gap between the count column (saturated at ~100 %) and the topology
column is the *silent rewiring* the older metrics could not see. On the
most diverse, fully-unseen corpus (GitHub repos) roughly **1 in 5**
schematics has its connectivity altered by the round-trip while the
count check still reports clean.

**Every** topology failure across all four corpora involves an
`X`/`U`/`S`/`T` symbol (vendor IC, multi-pin opamp, voltage-controlled
switch, transmission line). **Zero** pure-standard-element circuits
(R/C/L/V/I/D/Q/M/J) drift. Two takeaways:

1. Standard schematics round-trip with **perfect** connectivity.
2. The drift is exactly the documented multi-pin-vendor-symbol weakness
   (a symbol whose `.asy` is not on the LTspice search path comes back
   with a different pin list). Set `LTSPICE_ASY_SEARCH_PATH` /
   `--asy-dir` to resolve it, or keep the original `.asc` as the source
   of truth.

`check_circuit` (MCP) and `ltspice-convert --check` now surface this as
a `topology drift` warning, so an AI agent's "ship" decision can
distinguish the safe regime from the unsafe one instead of trusting a
count check that is blind to rewiring. The check has **no false
alarms**: 1-WL never reports drift on a circuit that actually
round-tripped correctly (textbook corpus: 110/110 clean).

Reproduce against your own corpus with the LAB harness
`bench/baseline.py` (gitignored, not shipped) — its
`roundtrip_topology_match` column.

## Headline results — v0.3.13

Both round-trip arms now hit 100 % on all three real-world corpora
(component-count match), and the GND-pin topology proxy clears 98 %
on each.

| Metric | Corpus | files | v0.3.7 | **v0.3.13** |
|---|---|---:|---:|---:|
| `.asc → netlist → .asc` count | GitHub repos     |  720 |  96.7 % | **100.00 %** |
| `.asc → netlist → .asc` count | LTspice Examples |  100 |  99.0 % | **100.00 %** |
| `.asc → netlist → .asc` count | LTspice Applications | 4099 | 100.0 % | **100.00 %** |
| `netlist → schemdraw → netlist` count | GitHub repos | 720 |  ~85 % | **~100 %** (300-sample) |
| `netlist → schemdraw → netlist` count | LTspice Examples | 100 |  85.0 % | **100.00 %** |
| `netlist → schemdraw → netlist` count | LTspice Applications | 4099 | ~99 % | **100.00 %** (300-sample) |
| GND-pin position preservation | GitHub repos | 11075 GND-pins | 96.4 % | **98.93 %** |
| GND-pin position preservation | LTspice Examples |  1829 | 99.6 % | **99.73 %** |
| GND-pin position preservation | LTspice Applications | 50503 | 99.8 % | **99.78 %** |
| schemdraw script exec failure | GitHub repos     |  720 |  ~5 % |  **0.14 %** (1 file: matplotlib mathtext) |

### What landed in v0.3.8 - v0.3.13

| Version | Fix | Effect |
|---|---|---|
| **0.3.8** D2 | `NetlistParser` accepts modelless Q/J/M (LTspice default-NPN form); InstName prefix-fix triggers for any non-matching first letter | GitHub-corpus count 96.7 → 99.6 % |
| **0.3.9** E1 | MCP server gains `check_circuit` / `info_circuit` (round-trip lint + stats for AI agents); `asy_search_dirs` on the conversion tools | Agent self-validation loop |
| **0.3.10** D3 | `* @sym=` hint emitted **before** the component line (parser ties hints to NEXT component); 1-pin X-prefix vendor symbols emit `<name> <node> <subckt>` (was 2-token, dropped) | `.asc ↔ .cir` 100 % on all corpora |
| **0.3.11** F1 | Schemdraw extractor accepts K (mutual inductance) and A (digital) directives, not only `.`-prefixed; `.subckt` body re-emitted with real newlines | Examples schemdraw arm +13 pt |
| **0.3.12** G2 | Multi-pin SUBCIRCUIT FLAG fallback uses single-column monotonic-distance layout (was 4×N grid with tied distances); pin order survives round-trip without `.asy` | GitHub GND-pin 96.4 → 98.93 % (278 fewer drifts) |
| **0.3.13** H2 | `schemdraw_to_cir` recognises `elm.Coax()` as transmission line (was silently dropped) | Schemdraw arm 100 % on all corpora |

### Remaining residual (sub-percent)

- **GND-pin 0.2-1.1 %** drift on multi-pin vendor symbols whose `.asy`
  is not on the default LTspice search path.  Set
  `LTSPICE_ASY_SEARCH_PATH` (or pass `--asy-dir`) for the affected
  third-party libraries and these resolve.  Without the .asy, the
  v0.3.12 column fallback keeps **index order** but does not match
  the original wire-end coordinates.
- **0.14 %** schemdraw exec failure on a single GitHub-corpus file
  whose component label embeds `{Tperiod}`-style parameter braces
  that matplotlib's text renderer parses as LaTeX.  Not a converter
  bug per se; would require either escaping `{`/`}` in
  `_sanitize_label` (loses lint round-trip of `{PARAM}`
  references) or switching matplotlib's mathtext off globally.
  Deferred.



## Methodology

For each `.asc` in the sample:

1. **`asc_to_netlist`**: parse the .asc and emit a SPICE netlist
   (pure-Python path, no LTspice.exe).
2. **`netlist_to_asc`**: round-trip the netlist back through
   `NetlistToAsc().convert_string()`.
3. **Round-trip count match**: re-parse the round-tripped .asc and
   verify the **component count** is preserved.

We also measure two stricter checks:

- **Symbol-kind drift**: count how many original SYMBOL kinds (e.g.
  `ind2`, `schottky`, `pnp`, `polcap`) disappear from the regenerated
  .asc.  This captures whether the schematic *looks* the same when
  re-opened in LTspice.
- **GND-connectivity preservation**: rebuild the bipartite
  component-to-node graph and check whether each component's set of
  ground-touching pin positions matches between original and
  round-tripped netlists.  This is a node-rename-invariant proxy for
  topology preservation.

Per-source sample size is 100 (count) / 80 (kind + GND) circuits, drawn
deterministically (`random.seed(0)`).  `github_repos` is the most
diverse and the closest proxy for unseen real-world input;
`textbook` was iterated against during development and is a noisy
upper bound.

The script is `bench/baseline.py` (gitignored; depends on a
LAB-private corpus path).

## Textbook training-data corpus (LAB-private)

**Distribution policy**: This corpus is the converter's training
material.  **Only the converter itself is published**; the textbook
circuits stay on the author's LAB-private storage and are NOT
redistributed via this repo.  Migration of any `.cir` from the
corpus into `tests/fixtures/` would constitute redistribution and is
explicitly out of scope -- the public test suite uses only
author-authored minimal circuits the author wrote from scratch.

The numbers below describe how the public converter performs against
that private corpus.  You cannot reproduce them bit-for-bit on your
own machine, but the converter behaviour they measure is exactly
what `pip install git+https://github.com/ksugahar/ltspice-converter`
delivers.

A second benchmark consumes a small, curated set of circuits harvested
straight from electrical engineering textbooks. Each entry is a
self-contained `.cir` of 3–20 components that exercises a specific
SPICE feature (DC node analysis, RC/LCR/CR filter, diode companion
model, CMOS SUBCKT, ...).  Provenance for every circuit is recorded
in `circuit_db_<book>.json`:

```json
{
  "book": "MATLABで学ぶ 回路シミュレーションとモデリング",
  "author": "菅原光俊",
  "publisher": "鳥影社",
  "entries": [
    {
      "page": 158, "figure": "Fig 4.2.1",
      "chapter": 4, "section": "4.2 DC解析",
      "circuit_type": "linear DC node analysis (current source driven)",
      "cir_file": "circuits_matlab/textbook_p158_Fig_4-2-1.cir",
      "lint_status": "PASS (no warnings)",
      ...
    }, ...
  ]
}
```

The `.cir` files themselves stay LAB-private (publisher copyright) and
live under `<mcp-server>/circuit/textbook/circuits_<book>/`.

### Round-trip pass rate (28 textbook circuits)

| Book | Entries | PASS (clean) | PASS (with warn) | FAIL | Date |
|---|---:|---:|---:|---:|---|
| MATLAB で学ぶ                | 17 | 16 (94.1%) | 1 | 0 | 2026-05-19 v0.3.6 |
| Python で学ぶ (新版)         | 11 | 10 (90.9%) | 1 | 0 | 2026-05-19 v0.3.6 |
| **All textbook circuits**    | **28** | **26 (92.9%)** | **2** | **0** | 2026-05-19 v0.3.6 |
| MATLAB で学ぶ                | 17 | **17 (100%)** | 0 | 0 | 2026-05-19 **v0.3.7** |
| Python で学ぶ (新版)         | 11 | **11 (100%)** | 0 | 0 | 2026-05-19 **v0.3.7** |
| **All textbook circuits**    | **28** | **28 (100%)** | **0** | **0** | 2026-05-19 **v0.3.7** |

The two previously-warned cases (both `.SUBCKT INV` CMOS-inverter
copies) were cleared in v0.3.7 by fixing two C5 lint false-positives:

- **C5-fp-1**: `X<name> ... <subckt_name>` (subcircuit invocation)
  was being looked up against the `.model` table.  X-class refs now
  resolve against `.subckt` definitions instead.
- **C5-fp-2**: `.model` declarations inside a `.subckt ... .ends`
  body were flagged as orphans because the orphan check scanned the
  whole netlist while the reference check only saw top-level
  components.  Both passes are now top-level-only, so subckt-internal
  models stay in scope where they belong.

Three new pytest regression tests
(`test_check_x_subckt_invocation_not_treated_as_model_ref`,
`test_check_model_inside_subckt_not_flagged_as_orphan`,
`test_check_undefined_subckt_warning_when_truly_missing`) lock the
fix in place.

## Headline results — v0.3.1 (after C2: `.asy` lookup + isolation zone)

`C2` — when emitting a multi-pin SUBCIRCUIT, place the SYMBOL in an
isolated coordinate band (x ≥ 4096), then emit one FLAG per pin at
the **canonical offset** reported by `AsyParser.get_terminal_offsets`
(the same lookup the asc parser uses on re-extraction). On the way
back through the pipeline asc_parser finds the same offsets and the
node names match exactly. Without an available `.asy` file the
generator falls back to the v0.3.0 compact grid (count-preserving
but topology-lossy).

### GND-connectivity preservation (node-rename-invariant topology proxy)

| Source | v0.3.0 | v0.3.1 | Δ |
|---|---|---|---|
| LTspice Examples                | 43.8 % | **48.8 %** | +5 pt |
| LTspice Applications            | 13.8 % | **76.2 %** | **+62.4 pt** |
| GitHub repos                    | 30.0 % | **46.2 %** | +16.2 pt |
| **Mean (240 samples)**          | **29.2 %** | **57.1 %** | **+27.9 pt** |

Count preservation unchanged (still 99.25 % mean). The 62-point jump
on LTspice Applications is the .asy lookup hitting on common vendor
ICs (LTC*, ADP*, ADA*, AD*). The smaller bump on GitHub repos
reflects symbols whose `.asy` files live in third-party libraries
that are not on the default LTspice search path.

## Headline results — v0.3.0 (after B1 partial)

`B1` — extend `Component` with `extra_nodes: List[str]` and propagate all
pin connections through NetlistParser → AscGenerator. AscGenerator now
emits a compact 4×5 grid of `FLAG + WIRE` lines around each multi-pin
SUBCIRCUIT, preserving the full pin list in the regenerated `.asc`.

Result: count preservation unchanged (99.25%); GND-connectivity
**unchanged** at the corpus level. Diagnostic showed the geometric
placement chosen by the auto-layouter often coincides with another
component's WIRE endpoint, so `asc_parser._estimate_terminals`
selects those external endpoints over the dedicated multi-pin FLAGs.
A proper fix needs either (a) a dedicated isolation zone in the
layouter or (b) `.asy` file lookup (planned for Phase C1). For now
B1 preserves the data structurally (pin info is in the netlist) but
the geometric round-trip remains lossy.

## Headline results — v0.2.1 (after Phase A: A1-A4)

### schemdraw script execution (new check)

| Source | Failures (v0.2.0) | Failures (v0.2.1) |
|---|---|---|
| 150 samples (50 each of LTspice examples / applications / github_repos) | **7 / 150** | **0 / 150** |

Three classes of bug fixed in Phase A:
- `AttributeError: 'X.end' not defined` when an Opamp/multi-pin element
  is the last in a parallel group (now falls back to `d.move()`).
- `ValueError: Axis limits cannot be NaN or Inf` on empty drawings (now
  emits an invisible guard line).
- 4-pin BJT substrate variants (`npn3`, `pnp4`) were re-emitted as
  `SYMBOL <Q-model-name>` instead of `SYMBOL npn3` (model name was
  being confused with symbol kind). Fixed by extending `SYMBOL_TO_SPICE`
  to recognise substrate variants.

## Headline results — v0.2.0 (after P1 + P4)

### `.asc → netlist → .asc` round-trip (component count)

| Source | v0.1.0 | v0.2.0 | Δ |
|---|---|---|---|
| textbook (training-adjacent)    | 100%   | 100%   | 0   |
| LTspice Applications            |  98%   | **100%** | +2  |
| LTspice Examples                |  94%   |  **99%** | +5  |
| GitHub repos (unseen)           |  79%   |  **98%** | +19 |
| **Mean (400 samples)**          | **92.75%** | **99.25%** | **+6.5** |

### `netlist → schemdraw script → netlist` round-trip

| Source | v0.1.0 | v0.2.0 | Δ |
|---|---|---|---|
| textbook                        | 100%   | 100%   | 0   |
| LTspice Applications            |  96%   |  **98%** | +2  |
| LTspice Examples                |  74%   |  **80%** | +6  |
| GitHub repos                    |  60%   |  **82%** | +22 |
| **Mean (200 samples)**          | **82.5%** | **90.0%** | **+7.5** |

The schemdraw arm shares the `asc_to_netlist` stage with the
.asc-round-trip arm, so P1 vendor-symbol handling improved both.

### Symbol-kind preservation (240 samples, after P4)

Drift events (smaller is better):

| Symbol kind | v0.1.0 | v0.2.0 | Δ |
|---|---:|---:|---:|
| `ind2`                          | 91 |   0 | -91 |
| `schottky`                      | 56 |   0 | -56 |
| `pnp`                           | 47 |   4 | -43 |
| `cap` (variant losses)          | 23 |  34 |   ? |
| `polcap`                        | 11 |   0 | -11 |
| `npn3`                          | 13 |  13 |   0 |
| `LTspiceControlLibrary\*` (15 kinds) |  ~90 | ~85 |  -5 |
| **Files with any drift**        | 195 | 145 | -50 |

The `cap`/`res`/`voltage` rows are mostly count side-effects, not
actual variant losses — investigated separately.

### GND-connectivity preservation (new metric in v0.2.0)

This is the strictest check we run, and it shows the **biggest
remaining weakness** — topology drift in regenerated schematics:

| Source | GND-connectivity preserved |
|---|---|
| LTspice Examples                |  43.8% |
| LTspice Applications            |  13.8% |
| GitHub repos                    |  30.0% |

Even when the component count is exactly right, the regenerated .asc
often has different pin connections to ground than the original.  The
biggest contributor is **multi-pin vendor symbols** (LTC*/ADP*/ISO*
ICs with 5–20 pins): the converter currently preserves the
**count** of such symbols but reduces them to 2-pin generic blocks
because it does not have access to the symbol's `.asy` definition.

If your use case requires faithful topology preservation for vendor
ICs, the workaround today is:

- Keep the original `.asc` as the source of truth and use the converter
  for **netlist generation only** (asc → netlist, which is 100% across
  all sources).
- For the reverse direction (netlist → asc), use the converter only on
  circuits built from standard LTspice library symbols (R/C/L/V/I/D/Q/M/J).

## Known failure modes (after v0.2.0)

Even at 98–100% count-preservation pass rate, some classes of circuit
are weak:

- **Multi-pin vendor IC topology** — counts are preserved but pin
  layout collapses to 2 pins.  See [GND-connectivity preservation]
  (#gnd-connectivity-preservation-new-metric-in-v020) above.
- **`.subckt` body inclusion** — when a circuit embeds a `.subckt`
  block in the same .asc, only the invocation is preserved on
  round-trip; the body may be lost.
- **`LTspiceControlLibrary\\*`** math/logic blocks — currently treated
  as 5-pin opamps but the math operation identity (`and`, `mul`, etc.)
  is not fully recovered.
- **`npn3` (BJT with substrate symbol)** — recovered as `npn` (3-pin
  view), substrate connection lost.
- **schemdraw script execution errors on some inputs** — observed in
  ~6% of github_repos samples (XML parse errors, axis-limit NaN,
  unknown element types).  These produce no output but raise a
  Python exception, so they are at least visible to the caller.

These will be addressed in subsequent releases.  Run
`python bench/baseline.py` against your own corpus and file an
[issue](https://github.com/ksugahar/ltspice-converter/issues) with the
failing `.asc` if you hit one of these or a new class.

## Reproducibility

```bash
pip install ltspice-converter[test]
git clone https://github.com/ksugahar/ltspice-converter
cd ltspice-converter
# Edit bench/baseline.py:CORPUS_ROOTS to point at your own .asc corpus
python bench/baseline.py --per-source 100 --out bench/my_baseline.json
```
