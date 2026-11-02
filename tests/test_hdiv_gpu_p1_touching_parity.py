"""Device flat-Gram P1 self/touching rule against the CPU ChargeGram entry.

Guards: the device flatEntry takes the same rule as the CPU
HighOrderTetEntryStrategy on identical staged inputs: the endpoint-graded P1
outer rule for self and vertex-touching pairs of affine P1 hosts (BDM1), the
established rule for every pair of a BDM2 host, and the far rule first.
Dropping the device touching branch reproduces the fixed-rule face self entry
(about 0.5 % off); widening it to BDM2 hosts moves those entries.
Assumes: the CPU entry is the reference; both sides share the pairEntry
arithmetic and production potentials, so this is parity, not independent
accuracy. No bitwise claim. PARITY_REL is a PROPOSED detector, not an adopted
criterion. Optional GPU lane: skips visibly without the native binding.
"""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
cuda = pytest.importorskip("radia._hdiv_cuda", reason="Optional isolated native GPU binding")
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402
from radia import hdiv_gpu  # noqa: E402
from radia.vim._vim import _prepare_gpu_charge_inputs, _rp  # noqa: E402

PARITY_REL = 1e-10
_GEOMETRY = ("cell_vertices", "face_vertices", "charge_host", "charge_kind", "charge_exponents",
             "tet_rule", "tri_rule", "tet_low_rule", "tri_low_rule", "far_factor")


def _staged(tmp_path, order):
    mesh = MakeStructured3DMesh(hexes=False, nx=1, ny=1, nz=1)
    path = tmp_path / f"cube-{order}.vol"
    mesh.ngmesh.Save(str(path))
    with ng.TaskManager():
        mesh = ng.Mesh(str(path))
        fes = ng.HDiv(mesh, order=order)
        _, inputs = _prepare_gpu_charge_inputs(fes, backend="gpu")
        tables, geometry = hdiv_gpu._device_charge_tables(path, fes, cuda)
    arguments = hdiv_gpu._device_gram_arguments(tables, inputs, geometry, None)
    return inputs, {name: arguments[name] for name in _GEOMETRY}


def _touching(device, i, j):
    """Direct-pair host vertices with bitwise-equal coordinates (CPU P1TouchHosts)."""
    def vertices(c):
        kind, host = device["charge_kind"][c], device["charge_host"][c]
        table = device["face_vertices"].reshape(-1, 3, 3) if kind else device["cell_vertices"].reshape(-1, 4, 3)
        return {tuple(v) for v in table[host]}
    return bool(vertices(i) & vertices(j))


@pytest.mark.parametrize("order", [1, 2])
def test_device_flat_entries_match_cpu_charge_gram(tmp_path, record_property, order):
    inputs, device = _staged(tmp_path, order)
    gram = _rp._ChargeGramHMatrix(**inputs)
    n = len(device["charge_host"])
    rows, cols = (a.ravel().astype(np.uint64) for a in np.meshgrid(np.arange(n), np.arange(n), indexing="ij"))
    actual = cuda._diagnostic_flat_charge_entries(**device, rows=rows, cols=cols, backend="gpu").reshape(n, n)
    expected = np.array([[gram.entry(i, j) for j in range(n)] for i in range(n)])
    scale = np.sqrt(np.outer(np.abs(np.diag(expected)), np.abs(np.diag(expected))))
    relative = np.abs(actual - expected)/scale
    touching = np.array([[_touching(device, i, j) for j in range(n)] for i in range(n)])
    worst = np.unravel_index(np.argmax(relative), relative.shape)
    record_property("charges", n)
    record_property("max_relative", float(relative.max()))
    # Both classes must be present, or the counterexample is unobservable.
    assert touching.any() and (~touching).any(), touching
    assert relative.max() <= PARITY_REL, (
        f"order={order} pair={worst} touching={touching[worst]} actual={actual[worst]!r} "
        f"expected={expected[worst]!r} relative={relative[worst]:.3e} limit={PARITY_REL:.0e}")


def test_diagnostic_entries_fail_loudly(tmp_path):
    _, device = _staged(tmp_path, 1)
    n = len(device["charge_host"])
    pair = np.zeros(1, dtype=np.uint64)
    with pytest.raises(ValueError, match="backend"):
        cuda._diagnostic_flat_charge_entries(**device, rows=pair, cols=pair, backend="cpu")
    with pytest.raises((ValueError, RuntimeError), match="out of range"):
        cuda._diagnostic_flat_charge_entries(**device, rows=np.array([n], dtype=np.uint64), cols=pair, backend="gpu")
    with pytest.raises(ValueError, match="shape"):
        cuda._diagnostic_flat_charge_entries(**device, rows=pair, cols=np.zeros(2, dtype=np.uint64), backend="gpu")
    with pytest.raises(TypeError):
        cuda._diagnostic_flat_charge_entries(**device, rows=pair.astype(np.int64), cols=pair, backend="gpu")
