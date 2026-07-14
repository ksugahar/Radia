"""Pytest configuration for the heavy Radia validation suite.

`tests/` is the lightweight debug/CI gate.  `validation_test/` is the
operator-triggered suite for solver-heavy, GUI, Cubit, golden, benchmark, and
cross-validation checks.  Reuse the root test configuration so imports,
markers, and optional-dependency handling stay identical.
"""

from tests.conftest import *  # noqa: F401,F403


def _is_compute_host():
    import platform

    node = platform.node().lower().split(".", 1)[0]
    return node == "mdx" or node.startswith("mdx-") or node == "hibino" or node.startswith("hibino-")


def pytest_collection_modifyitems(config, items):
    """Mark validation tests and keep compute-only jobs on mdx/hibino."""
    import pytest
    from tests.conftest import pytest_collection_modifyitems as _root_modify

    _root_modify(config, items)
    compute_skip = pytest.mark.skip(reason="compute_host validation runs only on mdx or hibino")
    for item in items:
        item.add_marker(pytest.mark.validation)
        if item.get_closest_marker("compute_host") and not _is_compute_host():
            item.add_marker(compute_skip)
