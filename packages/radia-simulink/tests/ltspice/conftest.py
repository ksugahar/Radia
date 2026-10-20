"""Shared pytest configuration.

The converter now defaults to the *auto* extraction backend (it uses
LTspice's own ``-netlist`` when LTspice.exe is installed). That is the
right production default, but it would make the test suite slow,
machine-dependent, and non-deterministic. This autouse fixture pins the
default back to the pure-Python extractor for every test by setting
``LTSPICE_NETLIST_PREFER=0``; tests that specifically exercise the
LTspice backend opt in explicitly and are skipped when LTspice is
absent (see ``test_ltspice_backend.py``).
"""
from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def _deterministic_pure_python_backend(monkeypatch):
    monkeypatch.setenv("LTSPICE_NETLIST_PREFER", "0")
