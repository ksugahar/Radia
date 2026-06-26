"""Golden test: SYMMETRY MODELS (1/2, 1/4, 1/8) for the HDiv-VIM demag via the IMAGE method.

The loop handling on cut meshes is already automatic (test_hdiv_vim_symmetry_loops.py).  This locks the
SYMMETRY-MODEL DEMAG VALUE: reflecting the reduced model's spherical-cap surface charge over the
reduction planes (sign = (-1)^#z-reflections, since sigma = n_z) reconstructs the full sphere, so the
1/2, 1/4, 1/8 models reproduce the FULL-sphere demag from ~1/2, 1/4, 1/8 the surface DOF.

See examples/vim/hdiv_demag_symmetry_image.py.
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "vim"))
import hdiv_demag_symmetry_image as sym  # noqa: E402


@pytest.fixture(scope="module")
def res():
    return sym.run(h=0.4)


def test_symmetry_models_reproduce_full_demag(res):
    """1/2, 1/4, 1/8 each reproduce the full-sphere demag (within 1%) via the image method."""
    for r in res["models"]:
        assert abs(r["match_vs_full"]) < 1e-2, \
            f"{r['model']} demag {r['demag']:.5f} off full {res['full_demag']:.5f} by {100*r['match_vs_full']:.2f}%"


def test_symmetry_models_use_fewer_dof(res):
    """The point of a symmetry model: fewer independent surface DOF (cap tris) the more it is reduced."""
    by = {r["model"]: r["cap_el"] for r in res["models"]}
    assert by["eighth"] < by["quarter"] < by["half"], \
        f"cap-tri count should shrink half->quarter->eighth: {by}"
    assert by["eighth"] < 0.4 * by["half"], f"1/8 should use far fewer DOF than 1/2: {by}"
