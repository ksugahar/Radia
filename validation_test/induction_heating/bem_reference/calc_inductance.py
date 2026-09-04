"""Compatibility launcher for the canonical IH inductance application.

The implementation lives in :mod:`radia.panels.calc_inductance`. Validation
must execute the shipped application backend instead of maintaining a second
copy of several thousand lines of BEM/PEEC coupling code.
"""

from radia.panels.calc_inductance import *  # noqa: F401,F403
from radia.panels.calc_inductance import main as _main


if __name__ == "__main__":
    raise SystemExit(_main())
