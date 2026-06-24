"""Golden test: the HDiv-VIM loop machinery is AUTOMATIC on symmetry-reduced (cut) meshes -- 1/2, 1/4,
1/8 models (verify-first, 2026-06-08).

The painful part of symmetry models in MSC/six-face surface-charge is the LOOP handling: the loop-star basis must be
constructed with a cohomology-aware `installCycle` on the cut domain (the symmetry planes introduce new
cycles).  In the HDiv-VIM the loops are simply ker(B) (B = charge map), so on ANY cut/reduced mesh they
are field-null BY CONSTRUCTION (N = B^T G B => N.loop = 0 for loop in ker B) with NO hand-crafted basis
and NO cohomology bookkeeping -- the loop count just adapts to the reduced topology.

This locks that on 1/2, 1/4, 1/8 reduced meshes the loops exist and are field-null to machine precision,
automatically.  (Scope: this is the LOOP machinery on cut meshes -- the demag VALUE of a symmetry MODEL
additionally needs the image-augmented Gram (the image method, like Radia IMA), which is separate and not
covered here.)
"""
import os
import sys

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "vim"))
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, OrthoBrick, Pnt  # noqa: E402

from conftest import hdiv_vim_dense_N_and_loops  # noqa: E402


def _reduced_sphere(box_lo):
    """Sphere intersected with a box -> full / 1/2 / 1/4 / 1/8 octant by the box lower corner."""
    g = CSGeometry()
    g.Add(Sphere(Pnt(0, 0, 0), 1.0) * OrthoBrick(Pnt(*box_lo), Pnt(2, 2, 2)))
    with ng.TaskManager():
        return ng.Mesh(g.GenerateMesh(maxh=0.5))


@pytest.mark.parametrize("tag,box_lo", [
    ("full", (-2, -2, -2)),
    ("half", (0, -2, -2)),
    ("quarter", (0, 0, -2)),
    ("eighth", (0, 0, 0)),
])
def test_loops_field_null_on_symmetry_reduced_mesh(tag, box_lo):
    """On full/1/2/1/4/1/8 meshes the ker(B) loops are field-null (N = B^T G B with the C++ Gram), so
    N @ loop = 0 BY CONSTRUCTION -- no hand-crafted loop-star, no cohomology bookkeeping."""
    with ng.TaskManager():
        N, loops, n_loop = hdiv_vim_dense_N_and_loops(_reduced_sphere(box_lo))
    assert n_loop > 0, f"{tag}: expected loops (ker B) on the reduced mesh"
    Nn = np.linalg.norm(N, 2)
    loop_res = max(np.linalg.norm(N @ loops[k]) for k in range(n_loop)) / Nn
    assert loop_res < 1e-8, f"{tag}: loops not field-null (residual {loop_res:.2e}) -- ker(B) should be exact"
