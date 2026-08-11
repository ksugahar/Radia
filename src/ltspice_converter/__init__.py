"""Compatibility import for the former standalone LTspice package.

Use :mod:`radia.ltspice` for new code. This module contains no circuit
implementation; it forwards the old namespace to Radia's canonical package.
"""
from __future__ import annotations

from radia import ltspice as _impl
from radia.ltspice import *  # noqa: F401,F403

__all__ = _impl.__all__
__version__ = _impl.__version__
__path__ = _impl.__path__


def __getattr__(name: str):
    return getattr(_impl, name)
