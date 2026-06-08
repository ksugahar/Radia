"""Compatibility shim (productionization M1, 2026-06-08).

The HDiv-type VIM nonlinear solver was PROMOTED to the production package ``radia.hdiv_vim`` (module
``radia.hdiv_vim._nonlinear``).  This thin shim re-exports it so existing examples/tests that still do
``import hdiv_demag_tet_nonlinear as nl`` keep working unchanged.

NEW CODE should use the production API:  ``from radia.hdiv_vim import solve_nonlinear_newton, ...``
(see docs/hdiv_vim/README.md).  This shim is transitional; the examples will migrate to the API in
M1-continued.
"""
from radia.hdiv_vim import _nonlinear as _m

# re-export every public + single-underscore name (tests use nl._scalar_fixed_point, _bh_table_funcs, ...)
globals().update({k: getattr(_m, k) for k in dir(_m) if not k.startswith("__")})
