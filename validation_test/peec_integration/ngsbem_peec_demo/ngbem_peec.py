"""Compatibility imports for the historical PEEC validation scripts.

The production implementation lives in :mod:`radia.ngsbem_peec`. Keeping a
second solver implementation under ``validation_test`` caused the validation
lane to exercise stale NGSolve BEM APIs instead of the code shipped to users.
"""

from radia.ngsbem_peec import *  # noqa: F401,F403
