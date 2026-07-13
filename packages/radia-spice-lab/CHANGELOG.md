# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.4.1] - 2026-07-13

### Added (balanced three-phase delta gate)

Added `balanced_three_phase_delta_rl_gate()` and MCP/CLI exposure for
checking balanced three-phase delta RL loads.  The gate reports line and
phase currents, active/reactive/apparent power, and topology hints for
publication-safe circuit validation.

### Added (prior-art search planning)

Added `patentability_search_plan()`, also exposed through MCP, to create
public prior-art search plans for circuit ideas.  It prepares Google
Scholar, Google Patents, J-PlatPat, and web query strings plus triage
questions for novelty / inventive-step review.  It is a non-legal search
aid, not a patentability opinion.
The returned dict also includes a `query_sets` grouping for callers that
want to iterate database-specific query lists uniformly.

### Fixed (schemdraw script inputs with `__file__`)

`schemdraw_to_netlist()` / `.py -> .cir` conversion now provides a
`__file__` value while executing schemdraw scripts.  The file-based
converter passes the real script path, so user scripts can resolve
adjacent helper files or output paths relative to their own directory.

### Changed (schemdraw label readability)

Generated schemdraw scripts now place vertical source and shunt labels
away from ground symbols and add a short branch stub before dropping
shunt components.  This keeps common RC/RLC/filter conversions more
readable while preserving topology round-trips.  Dense same-node
source/L/C shunt groups are also spread along a short bus instead of
being stacked on one coordinate.

.cir -> .asc layout now detects a source node that also carries shunt
capacitance and escapes adjacent series branches away from the source
symbol.  This keeps dense boost-input source/L/C/damping clusters readable
while preserving `.cir -> .asc -> .cir` topology.

Generated schemdraw labels now use compact single-line name/value text
(`R1 1k`) instead of cramped two-line labels (`R1\n1k`).  The reverse
converter also parses this compact label form so hand-authored
schemdraw scripts keep both component names and values.

### Quality

Public conversion checks now include CLI smoke paths, `.cir -> .asc ->
.cir`, `.cir -> schemdraw -> .cir`, `.asc -> .cir -> .asc -> .cir`,
and generated schemdraw `.py -> .cir` smoke coverage, plus a regression
case for dense same-node shunt layout and an ASC symbol-spacing gate for
source/L/C input clusters.

Schemdraw conversion tests now close matplotlib figures after each case, so
the full suite stays warning-free during longer conversion-learning runs.

## [0.4.0] — 2026-06-11

### Changed (J1: LTspice's own netlister is now the default extraction backend)

When `LTspice.exe` is installed, `asc_to_netlist` and the CLI now use
**LTspice's own `-netlist`** by default (auto), falling back to the
pure-Python extractor only when LTspice is absent. LTspice is the
ground truth for its own format, so this is the most faithful option;
it degrades gracefully everywhere else.

- `asc_to_netlist(asc, use_ltspice=None)` — `None` (new default) = auto,
  `True` = force LTspice, `False` = force pure-Python (deterministic).
- CLI: `--use-ltspice` / `--no-ltspice` (default auto).
- **Linting stays deterministic**: `check_circuit` / `--check` still
  default to the pure-Python round-trip (it measures the converter's
  self-consistency); pass `use_ltspice=True` / `--use-ltspice` to lint
  against LTspice instead.

How much the pure-Python path diverged from LTspice was invisible until
now because the round-trip ran pure-Python on *both* ends, cancelling
consistent bugs. Auditing the pure-Python extractor against LTspice as
an oracle (new `bench`) surfaced several real bugs — fixed below.

### Fixed (encoding: BOM-less UTF-16 .asc silently read as empty)

LTspice ships some schematics (e.g. the bundled `UniversalOpAmp*.asc`)
as **UTF-16 without a BOM**. The old reader only special-cased the
`0xFF 0xFE` BOM, so a BOM-less UTF-16 file decoded as UTF-8 into a
NUL-riddled string where `Version` / `SHEET` / `WIRE` never matched —
the schematic read as **empty** and silently produced a 0-component
netlist. New shared decoder `ltspice_converter._shared.decode_asc_bytes`
handles UTF-16 LE/BE (with or without BOM), UTF-8, and cp1252/latin-1,
and is used by both the CLI reader and the parser.

### Fixed (extraction fidelity vs LTspice — oracle-driven)

Audited against LTspice on the bundled "Examples" corpus, pure-Python
forward extraction matched LTspice's canonical topology only **55 %**
of the time. Four systematic bugs fixed:

- **Jumpers were emitted as components.** A `jumper` (net-tie) symbol is
  a zero-ohm short: its two pins now merge into one net and the symbol
  is dropped, matching LTspice (which removes it and merges the nodes).
- **Inline subcircuit parameters were mis-read as pins.** LTspice opamps
  netlist as `X1 in out V+ V- o level2 Avol=1Meg GBW=10Meg Rin=500Meg`;
  the parser took the last token as the subckt name, turning every
  `param=value` into a phantom node. It now splits nodes from the
  `key=value` parameter tail and finds the real subckt name.
- **Special-function `A` devices were mangled into `X§A…`.** An
  InstName starting with `A` (sample&hold, varistor, phase detector,
  modulator, Schmitt, …) is now kept as a real `A<n> …` line instead of
  being rewritten to a generic subcircuit, so it round-trips the same
  way LTspice writes it.
- **Symbol SPICE class is read from the `.asy` `SYMATTR Prefix`.** A
  symbol that is unmapped, or only generically mapped to `X`, now adopts
  the SPICE class LTspice itself declares for it — e.g. a crystal
  (`xtal`, `SYMATTR Prefix C`) netlists as a capacitor, not a generic
  subcircuit. (New `AsyParser.get_symbol_attr`.)

Result on the Examples corpus: strict topology match **55 % → 71 %**,
and **88 %** of circuits now have provably-correct wiring (the rest
differ only in the cosmetic op-amp class label `U` vs `X`, or are the
documented multi-pin-vendor-symbol case that needs the symbol's `.asy`
pin geometry — now covered for free by the LTspice backend). The
pure-Python self-round-trip is unchanged-or-better on every corpus
(Examples 92 → 94 %, Applications 97 %, GitHub 79 %, textbook 100 %).

### Fixed (LTspice modal dialogs blocked / flashed during `-netlist`)

Some schematics make `LTspice.exe -netlist` pop a **blocking** Yes/No/
Cancel message box — most commonly *"Duplicate overlapping components
found: C1 and C2 ..Delete?"*. Left alone it stalled the call for the
full timeout **and** flashed on the user's screen during batch
conversions. `_ltspice_emit_netlist` now runs a background watchdog
that auto-dismisses any LTspice modal owned by the spawned process,
clicking **No** (keep the components — netlist the schematic exactly as
drawn, never let LTspice modify it). On timeout the whole process tree
is killed so nothing is left orphaned. A "duplicate overlapping" case
that previously hung for 30 s now returns in ~2 s with both components
intact.

### Added

- `tests/test_ltspice_backend.py` (LTspice-gated: headless `-netlist`,
  backend agreement, the overlapping-dialog watchdog, no-orphan check,
  BOM-less UTF-16), `tests/conftest.py` (pins the suite to the
  deterministic pure-Python backend), `tests/test_extractor_fixes.py`
  (jumper / A-device, lib.zip-gated), and topology/subckt-param
  regression tests.
- `bench/baseline.py` gains an LTspice-as-oracle audit column.

## [0.3.14] — 2026-06-11

### Added (I1: node-rename-invariant topology check)

The round-trip metrics so far — component **count** and **GND-pin
position** — are necessary but not sufficient.  A circuit can keep
every component and every ground connection yet have its internal
wiring scrambled.  This is exactly what happens to a multi-pin vendor
symbol whose `.asy` is missing from the search path: it round-trips
with a different (often inflated) pin list, and the count check still
reports **100 % clean**.

New module `ltspice_converter.topology` collapses a netlist to a hash
that depends only on *how things are connected*, not on what the nodes
are named (Weisfeiler-Leman 1-WL colour refinement on the
component↔node incidence graph, with ground anchored and R/C/L pins
treated as unordered).  Two public functions:

- `topology_signature(netlist) -> str`
- `topology_equivalent(netlist_a, netlist_b) -> (bool, info)`

1-WL never colours two genuinely isomorphic graphs differently, so the
check has **no false alarms**: it never reports drift on a circuit that
actually round-tripped correctly.

Wired in two places:

- **`check_circuit` / `ltspice-convert --check`** now emit a
  `topology drift: node connectivity is NOT preserved on round-trip`
  warning (with the pin-incidence delta) — turning a previously
  silent corruption into an actionable signal for AI agents.
- **New MCP tool `compare_topology(netlist_a, netlist_b)`** — answers
  *"did my edit change the wiring?"*  Invariant to node renaming and
  benign R/C/L pin swaps, so after a value edit it returns
  `equivalent: true`; after an accidental rewire, `false`.

### Measured (what the new metric reveals)

`.asc → netlist → .asc → netlist` topology preservation, full corpora:

| Corpus | files | count match | **topology match** |
|---|---:|---:|---:|
| textbook              |  110 | 100 % | **100 %** |
| LTspice Applications  | 4099 | 100 % | **96.9 %** |
| LTspice Examples      |  100 | 99 %  | **92.0 %** |
| GitHub repos (unseen) |  720 | 100 % | **79.3 %** |

Every topology failure across all corpora involves an `X`/`U`/`S`/`T`
(vendor / multi-pin / switch / transmission-line) symbol — **zero**
pure-standard-element circuits (R/C/L/V/I/D/Q/M/J) drift.  So standard
schematics round-trip with perfect connectivity, and the check tells an
agent exactly when it is in the documented "keep the original `.asc` as
source of truth" regime (vendor symbol without `.asy`) versus the safe
regime.

### Added

- `tests/test_topology.py` — 9 regression tests (rename invariance,
  symmetric pin-swap tolerance, real-rewire / diode-flip / ground
  detection, value-edit preserves topology, pin-incidence reporting).
- `bench/baseline.py` gains a `roundtrip_topology_match` column.

## [0.3.13] — 2026-05-20

### Fixed (H2: transmission lines dropped on schemdraw round-trip)

`cir_to_schemdraw` emits ``elm.Coax()`` for SPICE T (transmission
line) components, but `schemdraw_to_cir.SCHEMDRAW_TO_SPICE` had no
``'Coax'`` entry -- so every T-line was silently dropped on the
``.cir -> schemdraw script -> .cir`` round-trip.  This was the
single remaining cluster of real component-count drift on the
schemdraw arm after the F1 K-directive fix (4 lab Examples files,
2 GitHub-repo files).

Fix: add ``'Coax': 'T'`` to the map and a 4-node emission path in
`_process_element` that reconstructs the SPICE T-line form
``T<name> A+ A- B+ B- [params]`` from schemdraw's 2-anchor Coax
(A+ = start, A- = 0, B+ = 0, B- = end).  Node assignment is not
byte-exact (schemdraw's 2-pin Coax cannot represent both ports
losslessly) but the component count and T-prefix survive, which is
what the lint and round-trip metrics measure.

### Performance

Schemdraw arm round-trip pass rate (count match, 300-file samples):

| Corpus | 0.3.12 | 0.3.13 |
|---|---|---|
| GitHub repos          | 99.33 % | **100.00 %** |
| LTspice Examples      | 98.00 % | **100.00 %** |
| LTspice Applications  | 100.00 % | 100.00 % |

The .asc <-> .cir arm remains at 100 % on all corpora
(D2/D3 work, v0.3.8 / v0.3.10).

### Added

- 1 pytest regression test
  (`test_schemdraw_arm_preserves_transmission_line`).

## [0.3.12] — 2026-05-20

### Fixed (G2: multi-pin SUBCIRCUIT pin-order drift without .asy)

When a multi-pin X-class component (vendor IC, custom subckt, ...) has
no resolvable `.asy` file on the LTspice symbol search path,
AscGenerator falls back to a generic FLAG layout.  The fallback used
to be a 4 x N grid of FLAGs around the symbol body, but multiple
grid positions shared the same Manhattan distance from the symbol
centre.  On re-extraction asc_parser._estimate_terminals orders pins
by ascending Manhattan distance, so the tied positions got reshuffled
into an arbitrary order -- silently moving e.g. the GND pin of a
6-pin IC from index 5 to index 4.

Replaced the grid with a single-column layout where pin i sits at
offset (DX, DY*(i+1)) from the symbol centre.  Manhattan distance is
then strictly monotonic in i, the round-trip preserves index order,
and the GND-pin metric stops reporting spurious topology drifts.

### Performance

GND-pin position preservation (all components with a GND pin):

| Corpus | 0.3.11 | 0.3.12 |
|---|---|---|
| GitHub repos          | 10679 / 11075 = 96.42 % | **10957 / 11075 = 98.93 %** (+2.5 pt) |
| LTspice Examples      |  1821 /  1829 = 99.56 % |  1824 /  1829 = 99.73 % |
| LTspice Applications  | 50382 / 50503 = 99.76 % | 50394 / 50503 = 99.78 % |

File-level GND-clean rate on GitHub repos: 524 -> **624 files** clean
out of 720 (+100 files where no GND-pin drifts at all).

Count-preservation pass rate unchanged at 100 % on all three corpora
(no regression from the layout change).

### Added

- 1 pytest regression test
  (`test_multi_pin_subckt_pin_order_preserved_without_asy`).

## [0.3.11] — 2026-05-20

### Fixed (F1: schemdraw arm round-trip — K-directive + multi-line .subckt)

Two bugs that prevented the ``netlist -> schemdraw script -> netlist``
arm from preserving component count.  Both originated in
`schemdraw_to_cir.py`: the script extractor was too narrow about which
labels count as SPICE directives, and the line-formatter packed
multi-line directives into a single physical line.

- **F1-1 (K / A directives silently dropped)**:
  `cir_to_schemdraw` already emits K (mutual inductance) and A
  (digital primitive) statements as Annotate labels alongside the
  ``.tran`` / ``.ac`` / ``.model`` / ``.subckt`` family.  But
  `schemdraw_to_cir._collect_directives_and_node_names` only picked up
  labels whose first character was ``.``, so every K and A statement
  was silently dropped on the round-trip.  K losses alone accounted
  for 23 / 31 missing components on a 100-file LTspice Examples
  sample.  Fix: also accept ``K``/``A`` as a first-character match.

- **F1-2 (multi-line .subckt body collapsed to one line)**:
  `NetlistParser` packs a ``.subckt ... .ends`` block into a single
  `SpiceDirective.text` with embedded real newlines.
  `cir_to_schemdraw` escapes the newlines to literal ``\n`` so the
  Python string literal stays on one source line in the Annotate
  call.  But the extractor's `_format_netlist` wrote the directive
  back verbatim -- newlines still escaped -- so the regenerated
  netlist had the entire DIAC / TRIAC subckt body packed into one
  physical line.  Downstream tooling (lint, count-based round-trip
  checks) saw zero internal components.  Fix: unescape ``\n`` to a
  real newline on the way out so every internal component reappears
  on its own SPICE line.

### Performance

Schemdraw arm round-trip pass rate (full LAB-private corpora):

| Corpus | files | 0.3.10 (200-file sample) | 0.3.11 (full) |
|---|---|---|---|
| GitHub repos          |  720 | 191/200 = 95.5 % | **711/720 = 98.75 %** |
| LTspice Examples      |  100 |  85/100 = 85.0 % | **98/100 = 98.0 %** (+13 pt) |
| LTspice Applications  | 4099 | 198/200 = 99.0 % | **4086/4099 = 99.68 %** |

The Examples bump is dominated by F1-1 (K-statement preservation);
the GitHub-corpus bump and the dimmer.asc-style cases are mostly
F1-2 (multi-line .subckt restoration).

The dense-`.asc -> netlist -> .asc` arm remains at 100 % on all
four real-world corpora (D1-D3 work, v0.3.8-v0.3.10).

### Added

- 2 pytest regression tests:
  `test_schemdraw_arm_preserves_k_directive`,
  `test_schemdraw_arm_preserves_multiline_subckt`.

## [0.3.10] — 2026-05-19

### Fixed (D3: last 3 GitHub-corpus failures + a 1-fail Examples regression)

Two distinct bugs surfaced by tracing the residual 3 failing files in
the LAB GitHub `.asc` corpus.  Both stemmed from places where the
extractor and parser disagreed on the canonical netlist line shape.

- **D3-1 (off-by-one ``@sym=`` hint association)**: The
  `NetlistExtractor` emitted the ``* @sym=<kind>`` variant hint
  **after** the component line, but `NetlistParser` consumes the hint
  via a ``pending_symbol_hint`` mechanism that ties it to the **next**
  component.  Every hint was therefore associated with the wrong
  component, and the trailing hint orphaned onto an unrelated next-class
  component.  On CLLC_Openloop the trailing ``* @sym=polcap`` hint
  attached to an ``Rload`` resistor → AscGenerator re-emitted it as a
  SYMBOL polcap → on re-extraction the polcap (C-class) prefix
  conflicted with the R-prefix name, so the v0.3.8 prefix-fix renamed
  it to ``C§Rload`` and the round-trip count dropped.  Fix: emit the
  hint **before** the component line (matches the parser's expectation
  and matches LTspice's own convention of "annotation precedes the
  thing it annotates").

- **D3-2 (1-pin X-prefix vendor symbol dropped on round-trip)**:
  PowerSim `CONST` (a 1-pin constant-source subcircuit) and similar
  1-pin vendor symbols were emitted as 2-token lines like ``X10 N031``
  -- no trailing subckt name.  `NetlistParser` then dropped them at
  the ``len(parts) < 3`` guard.  Fix: in the 1-terminal emit branch,
  when the SPICE prefix is X, append ``sym.spice_model or sym.symbol_type``
  as the subckt name so the parser sees a complete X-statement.
  TYPE2_FRA (1 file) + CLLC_Openloop + DAB_Openloop (2 X instances)
  all now round-trip.

### Performance

Round-trip pass rate (component-count match) on three LAB-private
corpora:

| Corpus | 0.3.9 | 0.3.10 |
|---|---|---|
| GitHub repos (720 files) | 717/720 = 99.6 % | **720/720 = 100 %** (+0.4 pt) |
| LTspice Examples (100 files) | 99/100 = 99.0 %  | **100/100 = 100 %** (+1 pt, same bug) |
| LTspice Applications (4099 files) | 4099/4099 = 100 % | 4099/4099 = 100 % (no regression) |

All three corpora now pass at 100 %.

### Added

- 2 pytest regression tests:
  `test_sym_hint_ordering_does_not_misclassify_next_component`,
  `test_one_pin_x_subcircuit_round_trip`.

## [0.3.9] — 2026-05-19

### Added (E1: MCP lint + stats tools)

The MCP server now exposes the same `--check` (lint) and `--info`
(stats) logic that the CLI has, so AI agents (Claude Code, Cursor)
can validate their own generated SPICE in-conversation without
shelling out:

- **`check_circuit(text, fmt, asy_search_dirs?)`** — round-trip drift
  check + static netlist analysis (duplicate instance names, floating
  nodes, orphan / undefined `.model` and `.subckt` references,
  undefined `{PARAM}` references, unparsed lines).  Returns
  `{ok: bool, info: [..], warnings: [..]}`.
- **`info_circuit(text, fmt, asy_search_dirs?)`** — component-type
  counts, symbol kinds, `.subckt` block count, `.asy` resolution rate.
- `netlist_to_asc` and `asc_to_netlist` gain an `asy_search_dirs`
  argument so agents can route vendor-symbol resolution
  (LTspiceControlLibrary etc.) through MCP too.

Total MCP tools: **6** (was 4).  Typical agent loop:
generate → `check_circuit` → fix if warnings → re-check → ship.

### Refactored (internal)

Split the file-based `_check_one` / `_info_one` into text-based
`check_text` / `info_text` core functions; the file-based wrappers
now delegate.  Both the CLI and the MCP server use the text-based
core, eliminating duplication.

### README

Five polish items:
- CI / Python-version / MIT-license badges at top
- "What's new in v0.3.9" link to CHANGELOG
- v0.3.4 example pin bumped to v0.3.9
- MCP server section gains motivation paragraph + 6-tool table
- "Supported elements" expanded with B / E / G / F / H / S / T / K / U
  (the rows were always supported, just not listed)

### Added (tests)

4 new pytest tests for E1: `check_text` clean case, `check_text`
duplicate-instance warning, `info_text` per-class counts,
`mcp_tools_registered` (all 6 tools present).

## [0.3.8] — 2026-05-19

### Fixed (D2: two round-trip-drop bugs surfaced by GitHub corpus)

Two distinct cases where the .asc → netlist → .asc round-trip silently
dropped a component are now preserved. Both came from training-data
analysis on a 720-file GitHub `.asc` corpus.

- **D2-1 (modelless BJT/JFET/MOSFET)**: LTspice's `npn`/`pnp`/`njf`/
  `pjf`/`nmos`/`pmos` symbols allow omitting `SYMATTR Value`, in which
  case the canonical netlist line has no trailing model token (e.g.
  `Q1 N006 N005 0` instead of `Q1 N006 N005 0 NPN`).  `NetlistParser`
  used to require `len(parts) >= 5` for BJT/JFET and `>= 6` for MOSFET
  and silently `return None`-ed the modelless form.  The two-letter
  555-timer reset transistors, vintage guitar-pedal BJT amplifiers
  (boss-ge7-equalizer, dunlop-crybaby-wah, schaller-tremolo, ...) all
  hit this case.  New `len(parts) == 4` (Q/J) and `len(parts) == 5` (M)
  branches accept the modelless form and leave `value` empty so the
  regenerated .asc has no `SYMATTR Value` either — byte-equal round-trip.

- **D2-2 (InstName whose first letter is not a SPICE prefix)**: A user
  who names a resistor `NTC` (on a `res` symbol) ends up with a netlist
  line `NTC N001 N002 R={...}` whose first letter `N` is not a SPICE
  device prefix at all.  LTspice itself prefix-fixes such names to
  `R§NTC` in the canonical netlist.  The pure-Python `NetlistExtractor`
  prefix-fix already handled the *conflicting* case (e.g. `T3a` on an
  `ind2` symbol → `L§T3a` because T is a SPICE prefix used by another
  device family) but was gated by `name[0] in _SPICE_PREFIXES`, which
  excluded the more common "not a prefix at all" case.  Removed the
  gate so any non-matching first letter gets the canonical `§`-prefix.

### Performance

Round-trip pass rate (component-count match) on three LAB-private
corpora:

| Corpus | 0.3.7 | 0.3.8 |
|---|---|---|
| GitHub repos (720 files) | 696/720 = 96.7 % | **717/720 = 99.6 %** (+2.9 pt) |
| LTspice Examples (100 files)     | 99/100 = 99.0 % | 99/100 = 99.0 % (unchanged) |
| LTspice Applications (4099 files) | 4099/4099 = 100 % | 4099/4099 = 100 % (no regression) |

Drift cluster reduction on GitHub repos: Q (BJT) losses 25 → 0;
R/N-prefix-mismatch losses 5 → 0.  Remaining 3 failures cluster on
multi-pin SUBCIRCUITs with 0–1 pins and a `.subckt`-body resistor
edge case (D3 work).

### Added

- 2 pytest regression tests covering both D2 fixes
  (`test_modelless_bjt_round_trip`, `test_instname_prefix_fix_for_non_spice_letter`).

## [0.3.7] — 2026-05-19

### Fixed (C5 false-positives, subckt-related)

Two warnings emitted by `ltspice-convert --check` on perfectly valid
SUBCKT-using netlists are eliminated:

- **C5-fp-1**: `X<name> ... <subckt_name>` was previously pooled into
  the same `referenced_models` set as D/Q/M/J components, so the
  subckt name `INV` (or any other) showed up as "referenced but not
  defined inline".  X-class refs now have their own
  `referenced_subckts` set and are checked against `.subckt`
  definitions.  Two new warning paths are added (orphan `.subckt`,
  undefined `.subckt` reference) so the same diagnostic quality
  applies to subcircuits as it does to models.
- **C5-fp-2**: `.model` declarations inside a `.subckt ... .ends`
  body were scanned by the orphan check while the matching device
  references inside that body were not (they live inside an absorbed
  multi-line directive, invisible to the top-level component list).
  Result: every subckt-internal `.model` looked like an orphan.  The
  orphan-model scan is now restricted to **top-level** `.model`
  lines, matching the scope of the reference scan.

Textbook training-data corpus moves from 26/28 clean pass to
**28/28 clean pass** as a direct consequence (the two SUBCKT INV
copies were the only entries with warnings).

### Added

- 3 pytest regression tests covering both fixes plus a positive test
  for the new "undefined .subckt" warning path.

## [0.3.6] — 2026-05-19

### Added (B5: silent-drop warnings in convert mode)

`ltspice-convert <input> -o <output>` now runs the converted netlist
through NetlistParser and prints any unparsed line to stderr with the
same "did you mean ...?" hint that `--check` already used. Previously
a malformed line was dropped silently from the output and the user had
to diff source and target to notice. Exit code is unchanged
(`convert` still succeeds when the primary output is written; the
warning is informational).

### Added (Doc: `.subckt` round-trip example)

README's *`.subckt` round-trip* section shows a small DIAC subckt
that round-trips byte-equal in both directions, with a one-line
`diff` command for verification. A matching fixture
(`tests/fixtures/bidirectional/00_converter_test_subckt_diac.cir`)
plus `test_subckt_body_round_trip` locks the C3 (v0.3.3) fix in
place: any future regression that drops `.subckt` internals will fail
this test before it can ship.

## [0.3.5] — 2026-05-19

### Distribution

- **GitHub-only release.** PyPI publish workflow removed
  (`.github/workflows/release.yml` deleted). README install
  instructions switched to:

      pip install git+https://github.com/ksugahar/ltspice-converter

  Pin a release with `@v0.3.5`. Editable dev install unchanged
  (`pip install -e .[mcp,test]`).

### Added (C5: strict lint checks)

`ltspice-convert --check` now also runs static netlist analysis on
the input and reports:

- **Duplicate instance names** (e.g. two `R1`s at top level).
- **Floating nodes** (a node touched by only one component, usually a
  wiring mistake).
- **Orphan `.model` directives** (a model declared inline but never
  referenced by any device).
- **Undefined model references** (a device references a model name
  that no inline `.model` defines; tone-downed to skip well-known
  LTspice library models like `1N4148`, `2N3904`, etc.).
- **Undefined `{PARAM}` references** (`{Rval}` used without a matching
  `.param Rval=...`).

These are warnings; `--strict` promotes any of them to exit 1, which
makes the CLI a real lint gate for CI pipelines.

### Added (B4: actionable error messages)

- NetlistParser now tracks lines it could not classify in a new
  `unparsed_lines` attribute (1-based line numbers + original text).
- `--check` surfaces each unparsed line as a warning, with a "did
  you mean ...?" hint for common typos
  (`Resistor` -> `res (R<name> ...)`, etc.) and a generic SPICE
  prefix reminder for unrecognised first characters. Prevents silent
  drops on user-written netlists where a token like `Zorg foo bar`
  used to vanish without trace.

### Tests

8 new CLI tests covering all five C5 checks + B4 unparsed-line
surfacing + suggestion path. 55 tests total (up from 47); CI matrix
unchanged.

## [0.3.4] — 2026-05-19

### Added

- **`ltspice-convert` CLI** (new console script). End users can now
  drive every conversion from the shell without writing Python.

  - `ltspice-convert input.asc -o output.cir` -- single-file convert,
    with target format inferred from `-o` or from the input
    extension.
  - `ltspice-convert *.asc -o build/ --to cir` -- batch into a
    directory.
  - `ltspice-convert --check input.asc` -- round-trip lint mode that
    reports component-count drift, GND-pin position drift, and
    SYMBOLs whose `.asy` cannot be located. `--strict` promotes any
    warning to exit code 1.
  - `ltspice-convert --info input.asc [--json]` -- summary of a
    schematic / netlist (component count by type, symbol kinds,
    `.asy` resolution rate, `.subckt` blocks).
  - `--asy-dir DIR` (repeatable) and `LTSPICE_ASY_SEARCH_PATH` env
    var both feed the third-party symbol search path; CLI flag wins.

- **GitHub Actions reusable workflow**:
  [`docs/example-workflows/asc-check.yml`](docs/example-workflows/asc-check.yml)
  can be copied into any repository that stores `.asc` files to lint
  them on every PR with `ltspice-convert --check --strict`.

- 15 new pytest tests covering convert / check / info / batch /
  error paths.

## [0.3.3] — 2026-05-19

### Added

- **C3**: `.subckt ... .ends` bodies survive `.asc → netlist → .asc`
  intact. NetlistParser now accumulates every line between `.subckt`
  and `.ends` into a single multi-line directive instead of leaking
  the internal components into the top-level component list, which
  previously caused the body to vanish during regeneration.
  Verified byte-equal on `dimmer.asc` (DIAC: 15 lines, TRIAC: 11
  lines).
- `_apply_name_remap` uses negative lookbehind `(?<!§)` so a token
  already prefixed (e.g. `X§Q1` inside a subckt body) does not get
  re-prefixed into `X§X§Q1` on a second pass.

## [0.3.2] — 2026-05-19

### Added

- **C4**: `LTSPICE_ASY_SEARCH_PATH` environment variable lets users
  point the converter at third-party LTspice symbol libraries
  (`LTspiceControlLibrary`, `LTspicePowerSim`, custom `MyLib`, …).
  Multiple paths are separated by the OS path separator (`;` on
  Windows, `:` on Linux). The env var is consumed by both the asc
  parser (extraction side) and the asc generator (emission side), so
  the round-trip is self-consistent.
  `NetlistToAsc(asy_search_dirs=[...])` accepts the same list
  programmatically.

### Notes

The C4 env var does not, on its own, shift the public benchmark
because the benchmark measured against a baseline that already had
both directions silently dropping the same vendor symbols (matching
"both wrong" cases). With the env var set, the round-trip is
genuinely correct for those symbols but the benchmark surfaces
previously-masked failures elsewhere. For users with third-party
libraries, the user-visible behaviour is strictly better.

## [0.3.1] — 2026-05-18

### Added

- **C2**: multi-pin SUBCIRCUITs now round-trip with correct topology
  when LTspice's `lib.zip` (or any `.asy` search path) is available.
  AscGenerator places each SUBCIRCUIT in an isolated coordinate band
  (x ≥ 4096) and emits one FLAG per pin at the canonical pin offset
  reported by `AsyParser.get_terminal_offsets`. asc_parser uses the
  same lookup on re-extraction, so pins land on identical coordinates
  and the node names match exactly.

### Performance (GND-connectivity preservation, 240 samples)

| Source | v0.3.0 | v0.3.1 |
|---|---|---|
| LTspice Examples     | 43.8 % | **48.8 %**  |
| LTspice Applications | 13.8 % | **76.2 %**  |
| GitHub repos         | 30.0 % | **46.2 %**  |
| **Mean**             | **29.2 %** | **57.1 %** (+27.9 pt) |

Count preservation unchanged at 99.25 %. All 32 tests pass.

## [0.3.0] — 2026-05-18

### Added

- `Component.extra_nodes: List[str]` captures the full pin list of
  multi-pin SUBCIRCUIT invocations (e.g. 20-pin vendor ICs such as
  LTC2945) through NetlistParser → AscGenerator.
- AscGenerator emits a 4×5 grid of `FLAG + WIRE` lines around each
  SUBCIRCUIT so the regenerated `.asc` retains every pin label, not
  just the first two.
- 4-pin BJT / MOSFET substrate variants (`npn3`, `pnp4`, `nmos4`,
  `pmos4`, `njf4`, `pjf4`) are recognised by `SYMBOL_TO_SPICE` and
  preserved end-to-end. Previously the model name (`SB`, `NP`, …) was
  re-emitted as the SYMBOL kind.
- `* @sym=<kind>` netlist comment preserves LTspice SYMBOL variants
  (`ind2`, `schottky`, `pnp`, `polcap`, …) through `.asc → netlist →
  .asc`. Drift on a 240-file corpus dropped from 195 to ≈ 145 files.
- GitHub Actions CI: `pytest` matrix on Linux + Windows × Python
  3.10 / 3.11 / 3.12, plus an MCP import smoke test.
- GitHub Actions release workflow: builds + publishes to PyPI via
  OIDC trusted publishing on `v*` tag push.
- Issue templates (`bug_report.yml`, `feature_request.yml`) and
  `CONTRIBUTING.md`.
- `docs/BENCHMARKS.md` documents `.asc ↔ netlist`, schemdraw arm,
  symbol-kind drift, and GND-connectivity metrics with v0.1.0 →
  v0.3.0 comparison.

### Fixed

- schemdraw script execution: 7 out of 150 sampled GitHub-repo
  circuits failed `exec()`; all now run cleanly.
  - `AttributeError: 'X.end' not defined` on multi-pin elements at
    the end of a parallel group → runtime `getattr` guard with
    `d.move()` fallback.
  - `ValueError: Axis limits cannot be NaN or Inf` on empty drawings
    (source `.asc` had only directives, no SYMBOLs) → invisible
    guard line.
  - `ValueError: SVG backend only supports saving SVG format figures`
    on headless Linux → `.pdf` → `.svg` fallback in generated script.

### Performance (400-sample real-world corpus)

| Metric | v0.1.0 | v0.3.0 |
|---|---|---|
| `.asc → netlist → .asc` count match (mean) | 92.75 % | **99.25 %** |
| `netlist → schemdraw → netlist` count match (mean) | 82.5 % | **≥ 90 %** |
| schemdraw script exec failures | 7 / 150 | **0 / 150** |
| Symbol-kind drift events | 195 files | **≈ 145 files** |

### Known limitations (target Phase C+ work)

- GND-connectivity preservation: 14–44 % on multi-pin vendor IC
  topology. Geometric layout collisions defeat the multi-pin FLAG
  scheme; needs either layouter isolation zones or `.asy` file
  lookup.
- `LTspiceControlLibrary\*` math/logic blocks: requires `.asy`
  lookup to recover pin definitions.
- `.subckt` body inclusion: invocation preserved, body may be lost.

## [0.1.0] — 2026-05-18

Initial public release. Successor to the (now-private)
`ksugahar/circuit-converter`. Conversion-only scope: `.asc` ↔ `.cir`
↔ schemdraw Python scripts. Bundles an MCP server (`mcp-ltspice`).
