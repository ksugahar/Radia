"""Guards: affine sampling owns coordinates before releasing mapped points.

Assumes: NGSolve owns geometry; these are lifetime/stability contracts, not
an independent validation of the HDiv formulation or a new accuracy gate.
"""
from types import SimpleNamespace
import weakref

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from radia.vim import _vim  # noqa: E402


class _MappedPoint:
    def __init__(self, coordinates, copies):
        self.coordinates = coordinates
        self.copies = copies

    @property
    def point(self):
        return _BorrowedPoint(self)


class _BorrowedPoint:
    def __init__(self, mapped):
        self.owner = weakref.ref(mapped)

    def _mapped(self):
        mapped = self.owner()
        assert mapped is not None, "mapped owner died before coordinate copy"
        return mapped

    def __iter__(self):
        return iter(self._mapped().coordinates)

    def __array__(self, dtype=None, copy=None):
        mapped = self._mapped()
        assert copy is True, "borrowed coordinates require an explicit owning copy"
        mapped.copies.append(copy)
        return np.array(mapped.coordinates, dtype=dtype, copy=True)


def _rule(monkeypatch):
    rule = [SimpleNamespace(point=p) for p in ((0., 0., 0.), (1., 0., 0.), (0., 1., 0.))]
    monkeypatch.setattr(_vim.ng, "IntegrationRule", lambda *args: rule)
    return rule


def test_production_affine_copies_while_mapped_owner_is_alive(monkeypatch):
    rule = _rule(monkeypatch)
    copies = []

    def trafo(ip):
        u, v, _ = ip.point
        return _MappedPoint(np.array([u, v, u + v]), copies)

    origin, jacobian = _vim._ngsolve_affine(trafo, ng.TRIG, 2)
    assert origin.shape == (3,) and jacobian.shape == (3, 2)
    assert np.isfinite(origin).all() and np.isfinite(jacobian).all()
    assert copies == [True] * (2 * len(rule)), f"copy events: {copies!r}"


def test_production_affine_keeps_the_finite_consecutive_fit_guard(monkeypatch):
    rule = _rule(monkeypatch)
    copies = []
    calls = 0

    def trafo(ip):
        nonlocal calls
        shift = (calls // len(rule)) % 2
        calls += 1
        u, v, _ = ip.point
        return _MappedPoint(np.array([u + shift, v, u + v]), copies)

    with pytest.raises(RuntimeError, match="unstable after 16"):
        _vim._ngsolve_affine(trafo, ng.TRIG, 2)
    assert calls == 16 * len(rule), f"mapped-point calls: {calls}"
