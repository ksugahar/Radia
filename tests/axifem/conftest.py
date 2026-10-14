"""Make the pure-Python axifem reference implementation importable
from `tests/axifem/test_*.py` modules.

The reference lives under `_reference_python/` (single source of truth for
the pure-NumPy/SciPy Henrotte/Meeker prototype that the C++
`axifem` module was ported from).  Adding the directory to
sys.path here lets every test module write a plain
`from axifem_core import ...` instead of fragile relative-path
sys.path tricks.

The underscore prefix on `_reference_python/` keeps pytest's
default collector from descending into it (the modules themselves
are not pytest modules).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REF_DIR = os.path.join(HERE, "_reference_python")
if REF_DIR not in sys.path:
    sys.path.insert(0, REF_DIR)
