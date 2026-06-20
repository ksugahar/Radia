"""Test isolation for the HDiv/yano demag-backend GLOBAL.

``rad.set_demag_backend(...)`` mutates a module-global (``radia._demag_backend``).
A test that forces a backend must restore it, or it LEAKS into later tests -- which
is exactly how ``tests/feec/`` went chronically CI-red:

  ``test_hdiv_vim_hex_wedge`` forced ``"hdiv"`` and "restored" it with
  ``prev = rad.set_demag_backend("hdiv"); finally: rad.set_demag_backend(prev)`` --
  but ``set_demag_backend`` returns the NEW value (not the previous), so ``prev``
  was ``"hdiv"`` and the restore RE-APPLIED ``"hdiv"``.  pytest runs ``hex_wedge``
  before ``newton_vs_radia`` + ``pm_mixing`` (alphabetical), which then ran
  ``rad.Solve`` (default -> leaked ``"hdiv"``) on a raw-tet iron and hit
  ``demag_backend='hdiv' needs a mesh-backed soft iron``.

This autouse fixture pins the global to the ``"auto"`` default around EVERY feec
test, so a forced backend can never leak across tests (defence in depth -- the
hex_wedge save/restore is also corrected to capture the previous value directly).
"""
import pytest


@pytest.fixture(autouse=True)
def _isolate_demag_backend():
    try:
        import radia as rad
    except Exception:                      # radia unavailable -> feec tests skip anyway
        yield
        return
    rad.set_demag_backend("auto")          # clean start, regardless of any prior leak
    try:
        yield
    finally:
        rad.set_demag_backend("auto")      # clean end, regardless of test outcome
