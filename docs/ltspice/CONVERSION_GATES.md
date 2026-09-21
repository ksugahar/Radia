# LTspice conversion quality gates

Radia treats circuit conversion as an electrical-contract problem, not a text
formatting problem. A clean round-trip must preserve connectivity and device
meaning even when anonymous node names, drawing coordinates, or statement
ordering change.

## Public regression lane

The redistributable tests and fixtures live in `tests/ltspice`. Run them from
the Radia repository root:

```bash
python -m pip install -e ".[ltspice]"
python -m pytest tests/ltspice -q
```

The lane checks:

- `.cir -> .asc -> .cir` and `.asc -> .cir -> .asc` conversions.
- `.cir -> schemdraw Python -> .cir` script execution and recovery.
- Node-rename-invariant topology signatures.
- Ground and multi-pin connectivity, not only component counts.
- Subcircuits, controlled sources, behavioral expressions, switches, coupled
  inductors, transmission lines, and representative power circuits.
- BOM-less UTF-16 and other LTspice text-decoding cases.
- CLI exit codes, diagnostics, JSON summaries, and MCP wrappers.
- LTspice `.measure` and stepped-measure parsing schemas.

All committed circuit fixtures are small, author-created regression inputs.
Third-party circuit collections and vendor example libraries are not part of
the public test corpus.

## Conversion checks

Use the canonical CLI to gate a circuit before committing it:

```bash
radia-ltspice --check --strict path/to/circuit.asc
radia-ltspice --info --json path/to/circuit.asc
```

`--check` reports independent signals:

1. Parse and conversion success.
2. Component and directive retention.
3. Ground-pin consistency.
4. Node-rename-invariant topology equivalence.
5. Symbol-resolution warnings.

A component-count match does not override a topology mismatch.

## Backend contract

For `.asc -> netlist`, the default `auto` mode uses LTspice's own netlist
export when LTspice is installed. The pure-Python extractor remains available
for deterministic, cross-platform operation:

```bash
radia-ltspice input.asc -o output.cir --use-ltspice
radia-ltspice input.asc -o output.cir --no-ltspice
```

The two routes share the same output contract and topology checks. Tests that
must run without LTspice set `LTSPICE_NETLIST_PREFER=0` explicitly.

## Third-party symbols

Set `LTSPICE_ASY_SEARCH_PATH` or pass one or more `--asy-dir` arguments when a
schematic depends on `.asy` files outside the standard LTspice library:

```bash
radia-ltspice --asy-dir path/to/vendor/sym --check --strict circuit.asc
```

Missing symbol metadata is a warning or strict-mode failure because pin order
cannot be inferred safely from appearance alone.

## Adding a regression

When a conversion bug is fixed:

1. Reduce the circuit to the smallest author-created reproducer.
2. Add it under `tests/ltspice/fixtures`.
3. Add a focused assertion for the lost device or connection.
4. Run the full `tests/ltspice` lane.
5. Keep any nonredistributable source circuit outside the repository.

The implementation is canonical under `src/radia/ltspice`; compatibility
imports must never carry a second copy of converter logic.
