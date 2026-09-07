"""Compatibility import; the solver-neutral implementation owns this API."""
from ..optimization.line_search import armijo_backtracking

__all__ = ["armijo_backtracking"]
