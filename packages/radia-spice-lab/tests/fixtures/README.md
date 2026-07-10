# Test Fixtures

All circuits in this directory are author-authored minimal examples
designed to exercise the converter round-trip.  Nothing here is
derived from LTspice's bundled `Educational/` or `Applications/`
examples, from textbooks, or from third-party GitHub repositories.

## Layout

- **`bidirectional/`** — minimal `.asc` + `.cir` pairs.  Used by the
  pytest round-trip suite to verify that converting from one form to
  the other and back is information-preserving.
- **`converter/01_rc_lowpass/` ... `07_horizontal_cap/`** — seven
  per-circuit demo directories.  Each contains:
  - `test_<name>.asc` — original LTspice schematic
  - `test_<name>.cir` — equivalent SPICE netlist (when present)
  - `test_<name>.py` — generator/demo script (historical reference;
    the bare imports inside reflect the converter's pre-package layout
    and are kept for transparency, not as runnable examples)
  - `test_<name>.pdf` — schemdraw rendering for visual reference
- **`converter/netlist_to_asc/`** — additional `.cir` → `.asc` round-trip
  fixtures with their generator scripts.

Pytest is configured (`norecursedirs` in `pyproject.toml`) to skip this
directory during test collection so the historical generator scripts
do not pollute the test suite.

To exercise the converter on these fixtures from a fresh installation:

```python
from pathlib import Path
import ltspice_converter as lc

base = Path("tests/fixtures/bidirectional")
for asc_path in base.glob("*.asc"):
    asc_text = asc_path.read_text(encoding="utf-8")
    netlist = lc.asc_to_netlist(asc_text)
    print(asc_path.name, "->", len(netlist), "byte netlist")
```
