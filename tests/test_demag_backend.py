"""The runtime demag-backend selector + yano-type deprecation (radia.set_demag_backend).

The yano-type MSC demag backend (rad.Solve on hex/wedge soft iron) is DEPRECATED, superseded by the FEEC
HDiv-VIM (radia.vim).  It is exposed as an OPTIONAL, selectable runtime backend so it can be removed
cleanly: "yano" (legacy default, one-time DeprecationWarning) vs "hdiv" (forbids the yano path until the
HDiv-VIM is wired into rad.Solve).  These tests lock the API gate WITHOUT a real solve."""
import pytest
import radia as rad


@pytest.fixture(autouse=True)
def _restore_backend():
    """Always leave the global backend at the default so other tests are unaffected."""
    prev = rad.get_demag_backend()
    emitted = rad._yano_deprecation_emitted
    yield
    rad.set_demag_backend(prev)
    rad._yano_deprecation_emitted = emitted


def test_default_backend_is_yano():
    assert rad.get_demag_backend() == "yano"


def test_set_get_round_trip_and_return_previous():
    prev = rad.set_demag_backend("hdiv")
    assert prev == "yano"
    assert rad.get_demag_backend() == "hdiv"


def test_invalid_backend_raises():
    with pytest.raises(ValueError):
        rad.set_demag_backend("bogus")


def test_solverconfig_sets_backend():
    rad.SolverConfig(demag_backend="hdiv")
    assert rad.get_demag_backend() == "hdiv"


def test_hdiv_backend_requires_registered_mesh():
    """With 'hdiv' selected, rad.Solve on a container that has NO HDiv-registered soft-iron body
    (here a bogus handle) raises NotImplementedError redirecting to radia.vim.soft_iron_from_mesh --
    the FEEC HDiv-VIM needs the NGSolve mesh association (a registered iron dispatches the real
    HDiv-VIM solve; see tests/feec/test_hdiv_radsolve_dispatch.py)."""
    rad.set_demag_backend("hdiv")
    with pytest.raises(NotImplementedError):
        rad.Solve(0, 1e-4, 10, 0)


def test_yano_backend_emits_deprecation_once(monkeypatch):
    """The yano path emits a one-time DeprecationWarning.  Stub the C++ Solve so the gate/warning is
    exercised without a real solve."""
    calls = []
    monkeypatch.setattr(rad, "_cpp_Solve", lambda *a, **k: calls.append(a) or (None, 0))
    rad.set_demag_backend("yano")
    rad._yano_deprecation_emitted = False
    with pytest.warns(DeprecationWarning, match="yano-type MSC demag backend"):
        rad.Solve(0, 1e-4, 10, 0)
    # second call: no second warning (once per session), still dispatches
    import warnings
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        rad.Solve(0, 1e-4, 10, 0)
    assert not any(issubclass(w.category, DeprecationWarning) for w in rec)
    assert len(calls) == 2
