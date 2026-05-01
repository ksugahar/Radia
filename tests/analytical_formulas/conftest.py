"""Make the analytical_formulas test directory standalone-runnable.

Adds ``src/`` to ``sys.path`` so the tests can be invoked from a fresh
checkout without a global ``pip install -e`` step.
"""

import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
