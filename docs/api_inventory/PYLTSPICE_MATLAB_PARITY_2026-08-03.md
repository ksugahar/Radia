# PyLTSpice 6.0.1 to MATLAB compatibility inventory

Reference package: PyLTSpice 6.0.1 with spicelib 1.6.3. The reference was
installed in an isolated temporary environment for public-signature inspection;
its source is not copied into Radia.

## Definition of complete

The port is complete only when every public PyLTSpice top-level class and method
is classified in `matlab/pyltspice_api_compatibility.json`, every `native-matlab`
or `native-mex` entry has a MATLAB regression test, and the checked inventory
contains no `missing` entry.  snake_case compatibility aliases and the idiomatic
MATLAB camelCase names are both part of the contract.

## Current implementation

| Family | MATLAB surface | Status |
|---|---|---|
| RAW read | `radia.ltspice.RawRead`, `Trace` | Native MATLAB; real/complex, ASCII/binary, step slicing, export/table/CSV/Excel |
| RAW write | `radia.ltspice.RawWrite` | Native MATLAB ASCII writer; real/complex round-trip |
| Netlist editing | `radia.ltspice.SpiceEditor`, `SpiceCircuit` | Native MATLAB core editing and query operations; hierarchy and component-attribute parity remain incomplete |
| ASC editing | `radia.ltspice.AscEditor`, `SchematicEditor` | Native MATLAB core value/directive/component operations; graphical wire/move/scale parity remains incomplete |
| Simulation | `radia.ltspice.SimRunner`, `RunTask`, `SimCommander`, `LTspice` | Native MATLAB synchronous execution and batch/parallel Radia execution; asynchronous callback/task filtering remains incomplete |
| Log reading | `radia.ltspice.LTSpiceLogReader`, `LogReader` | Native MATLAB scalar measure API; stepped condition filtering and complex dataset transforms remain incomplete |
| Simulink | `radia.simulink.ltspiceSFunction` | Readable Level-2 MATLAB S-Function; arbitrary fixed-at-compile vector input width through `InputNames` |

## MEX promotion boundary

Do not make the complete Level-2 wrapper a MEX S-Function.  Profile warmed runs
first.  Binary RAW decoding/encoding, very large trace slicing, or persistent
native circuit state may be promoted to independently callable MEX functions
when measured transfer and parsing cost is material.  Such state uses checked
opaque `uint64` handles and remains owned by the Level-2 MATLAB wrapper.

## Remaining gates

1. Generate the exhaustive method manifest from the pinned Python packages and
   fail CI when the upstream signature changes.
2. Complete hierarchical netlist and graphical ASC editing.
3. Complete asynchronous `SimRunner`/`RunTask` behavior and callback semantics.
4. Complete multi-plot RAW, stepped-log query, operating-point, alias,
   back-annotation, and property readers.
5. Run MATLAB and Simulink E2E tests on LAB and 100号機 with the current ADI
   LTspice release.

Until all five gates pass, documentation and release notes must say
"PyLTSpice-compatible core" rather than "complete PyLTSpice port".
