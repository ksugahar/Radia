"""Shared execution-boundary fixtures for radia-mcp validation tests."""

from __future__ import annotations

import pytest


@pytest.fixture(scope="module")
def ngsolve_taskmanager():
    """Own NGSolve parallel state for modules that explicitly opt in."""
    try:
        import ngsolve as ng
    except ImportError:
        yield
        return

    with ng.TaskManager():
        yield
