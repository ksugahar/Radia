"""Fast contract for the reusable C++ ChargeGram build used by Optuna tuning."""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve.meshes import MakeStructured3DMesh

from radia.vim import ChargeGram


def test_charge_gram_rebuild_preserves_configured_demag_action():
    mesh = MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=2)
        _, gram, _ = ChargeGram(
            fes, eps=1.0e-12, leafsize=16, eta=2.0,
            _build_hmatrix=False)

    probe = np.random.default_rng(20260826).standard_normal(fes.ndof)
    gram.build_hmatrix(eps=1.0e-12, leaf=16, eta=2.0)
    first = np.asarray(gram.apply_configured_demag(probe, True), dtype=float)
    first_stats = dict(gram.stats())

    gram.build_hmatrix(eps=1.0e-12, leaf=16, eta=2.0)
    second = np.asarray(gram.apply_configured_demag(probe, True), dtype=float)
    second_stats = dict(gram.stats())

    assert first_stats["n_dof"] == second_stats["n_dof"] == gram.ndof()
    np.testing.assert_array_equal(second, first)
