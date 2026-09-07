"""Compatibility imports; linear CG belongs to linear algebra, not nonlinear CG."""
from ..matrix_solvers.krylov import _as_matvec, linear_conjugate_gradient

__all__ = ["linear_conjugate_gradient"]
