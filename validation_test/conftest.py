"""Pytest configuration for the heavy Radia validation suite.

`tests/` is the lightweight debug/CI gate.  `validation_test/` is the
operator-triggered suite for solver-heavy, GUI, Cubit, golden, benchmark, and
cross-validation checks.  Reuse the root test configuration so imports,
markers, and optional-dependency handling stay identical.
"""

from tests.conftest import *  # noqa: F401,F403


def _load_slow_nodeids():
    """Load measured slow node IDs owned by the validation suite."""
    from pathlib import Path

    path = Path(__file__).with_name("slow_nodeids.txt")
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return set()
    return {
        line.split("#", 1)[0].strip().replace("\\", "/")
        for line in lines
        if line.split("#", 1)[0].strip()
    }


SLOW_NODEIDS = _load_slow_nodeids()


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
        nodeid = item.nodeid.replace("\\", "/")
        item.add_marker(pytest.mark.validation)
        if nodeid in SLOW_NODEIDS:
            item.add_marker(pytest.mark.slow)
        if item.get_closest_marker("compute_host") and not _is_compute_host():
            item.add_marker(compute_skip)
