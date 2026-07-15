"""Fast API-contract tests for HDiv-VIM ChargeGram dispatch.

The pure 2D / hex / wedge shortcuts must preserve the public ``ChargeGram`` build
parameters.  They are the production paths for motor cross-sections and
structured cube benchmarks, so silently falling back to helper defaults can turn
a requested fast HACApK build into a much tighter one.
"""

import pytest

pytest.importorskip("ngsolve")
import radia.vim._vim as V  # noqa: E402


class _Element:
    def __init__(self, n_vertices):
        self.vertices = [object()] * n_vertices


class _Mesh:
    def __init__(self, dim, n_vertices=None):
        self.dim = dim
        self._elements = [] if n_vertices is None else [_Element(n_vertices)]

    def Elements(self, _kind):
        return self._elements

    def GetCurveOrder(self):
        return 1


class _FES:
    globalorder = 1

    def __init__(self, mesh):
        self.mesh = mesh


@pytest.fixture(autouse=True)
def _mock_native_topology_query(monkeypatch):
    monkeypatch.setattr(
        V, "_volume_vertex_counts",
        lambda mesh: frozenset(len(element.vertices) for element in mesh.Elements(None)))


@pytest.mark.parametrize(
    ("mesh", "helper_name"),
    [
        (_Mesh(2), "_build_charge_gram_2d"),
        (_Mesh(3, 8), "_build_charge_gram_hex"),
        (_Mesh(3, 6), "_build_charge_gram_wedge"),
    ],
)
def test_chargegram_shortcuts_forward_hacapk_params(monkeypatch, mesh, helper_name):
    calls = []

    def fake_helper(_fes, **kwargs):
        calls.append(kwargs)
        return "B", "G", "M"

    monkeypatch.setattr(V, helper_name, fake_helper)
    monkeypatch.setattr(V, "_configure_cpp_operator", lambda *args: args)

    result = V.build_charge_gram(_FES(mesh), eps=1.0e-4, leafsize=24, eta=1.75)

    assert result == ("B", "G", "M")
    assert len(calls) == 1
    assert calls[0]["eps"] == pytest.approx(1.0e-4)
    assert calls[0]["leafsize"] == 24
    assert calls[0]["eta"] == pytest.approx(1.75)
    assert calls[0]["materialize_mass"] is True


def test_chargegram_native_solve_skips_python_mass_materialization(monkeypatch):
    calls = []

    def fake_hex(_fes, **kwargs):
        calls.append(kwargs)
        return "B", "G", None

    monkeypatch.setattr(V, "_build_charge_gram_hex", fake_hex)
    monkeypatch.setattr(V, "_configure_cpp_operator", lambda *args: args)

    result = V.build_charge_gram(
        _FES(_Mesh(3, 8)), _materialize_mass=False)

    assert result == ("B", "G", None)
    assert calls[0]["materialize_mass"] is False


def test_chargegram_hex_shortcut_forwards_image_params(monkeypatch):
    calls = []

    def fake_hex(_fes, **kwargs):
        calls.append(kwargs)
        return "B", "G", "M"

    monkeypatch.setattr(V, "_build_charge_gram_hex", fake_hex)
    monkeypatch.setattr(V, "_configure_cpp_operator", lambda *args: args)

    V.build_charge_gram(
        _FES(_Mesh(3, 8)),
        eps=1.0e-4,
        leafsize=24,
        eta=1.75,
        image_masks=[1],
        image_signs=[-1],
    )

    assert calls == [{
        "eps": 1.0e-4,
        "leafsize": 24,
        "eta": 1.75,
        "image_masks": [1],
        "image_signs": [-1],
        "materialize_mass": True,
        "build_hmatrix": True,
    }]
