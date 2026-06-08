"""Compatibility shim (productionization M1, 2026-06-08).

The HDiv-type VIM CORE was PROMOTED to the production package ``radia.hdiv_vim`` (module
``radia.hdiv_vim._core``).  This thin shim re-exports it so existing examples/tests that still do
``import hdiv_demag_tet as tet`` keep working unchanged.

NEW CODE should use the production API:  ``from radia.hdiv_vim import build_demag, demag_factor, ...``
(see docs/hdiv_vim/README.md).  This shim is transitional; the examples will migrate to the API in
M1-continued.
"""
from radia.hdiv_vim import _core as _m

# re-export every public + single-underscore name (tests use tet._bary_tri, tet.C_TRI, tri_potential, ...)
globals().update({k: getattr(_m, k) for k in dir(_m) if not k.startswith("__")})
