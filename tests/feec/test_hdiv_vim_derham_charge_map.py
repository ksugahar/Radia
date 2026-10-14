"""Discrete de Rham gate for the production HDiv-VIM charge map.

HCurl is not a state unknown in the material solve, but its curl image must be
charge-free.  This test exercises NGSolve's mapped HCurl/HDiv spaces and the
same public ChargeGram construction used by production; it does not recreate
element orientations or Piola transforms in Python.
"""

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry, Pnt  # noqa: E402

from radia.vim import ChargeGram  # noqa: E402


@pytest.mark.parametrize("order", [1, 2])
def test_boundary_compatible_hcurl_curl_is_in_charge_map_kernel(order):
    """The production B map annihilates the discrete curl image for BDM1/2."""
    mesh = ng.Mesh(
        OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=1.5)
    )

    with ng.TaskManager():
        hcurl = ng.HCurl(mesh, order=order, dirichlet=".*")
        hdiv = ng.HDiv(mesh, order=order)

        vector_potential = ng.GridFunction(hcurl)
        free = np.asarray(list(hcurl.FreeDofs()), dtype=bool)
        coefficients = vector_potential.vec.FV().NumPy()
        coefficients[:] = 0.0
        coefficients[free] = np.random.default_rng(20260825 + order).standard_normal(
            int(free.sum())
        )

        magnetization = ng.GridFunction(hdiv)
        magnetization.Set(ng.curl(vector_potential))
        charge_map, gram, _ = ChargeGram(
            hdiv,
            eps=1e-10,
            leafsize=256,
            ho_far_factor=float("inf"),
        )

    m = np.asarray(magnetization.vec.FV().NumPy(), dtype=float)
    charge = np.asarray(charge_map @ m, dtype=float)
    demag = np.asarray(
        gram.apply_configured_demag(np.ascontiguousarray(m)), dtype=float
    )
    assert np.linalg.norm(m) > 1e-8
    assert np.linalg.norm(charge) / np.linalg.norm(m) < 5e-11
    assert np.linalg.norm(demag) / np.linalg.norm(m) < 5e-11
