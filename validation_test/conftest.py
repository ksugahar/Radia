"""Pytest configuration for the heavy Radia validation suite.

`tests/` is the lightweight debug/CI gate.  `validation_test/` is the
operator-triggered suite for solver-heavy, GUI, Cubit, golden, benchmark, and
cross-validation checks.  Reuse the root test configuration so imports,
markers, and optional-dependency handling stay identical.
"""

from tests.conftest import *  # noqa: F401,F403


def pytest_collection_modifyitems(config, items):
    """Mark every test in this tree as validation, then reuse root markers."""
    import pytest
    from tests.conftest import pytest_collection_modifyitems as _root_modify

    _root_modify(config, items)
    for item in items:
        item.add_marker(pytest.mark.validation)
