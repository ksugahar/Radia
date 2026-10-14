"""Compatibility imports for the historical PEEC bridge scripts.

The production Loop-Star bridge lives in :mod:`radia.ngsbem_interface`.
Validation uses that implementation so topology and solver fixes cannot drift.
"""

from radia.ngsbem_interface import *  # noqa: F401,F403
