"""Compatibility imports for the historical eddy-current validation scripts.

The production implementation lives in :mod:`radia.ngsbem_eddy`. Validation
must exercise that implementation rather than maintain a drifting copy of the
surface-current, FEM-BEM, and SIBC solvers.
"""

from radia.ngsbem_eddy import *  # noqa: F401,F403
