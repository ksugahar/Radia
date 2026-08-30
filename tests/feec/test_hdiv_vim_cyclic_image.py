"""Cyclic (N-fold rotational) images in the HDiv-VIM charge Gram.

The image method historically folded MIRROR images (a 1/2, 1/4, 1/8 reduced model reproduces the full
model).  A cyclic image instead maps the source body by a ROTATION about +z, so one sector of a
rotationally symmetric machine reproduces the whole ring.  Two things make it work:

* charges are SCALARS under a rotation that carries the magnetization along with the geometry
  (sigma(R x) = R M . R n = M . n), so a cyclic image needs no factor of its own and an alternating
  N/S pole pattern is just the sign (-1)^k on the existing image_signs;
* the eval point takes the INVERSE transform (|x - T(y)| = |T^-1(x) - y|).  A mirror is an involution
  so the historical code could not distinguish the two; a rotation is not.

ANCHOR (test_pi_rotation_matches_mirror_pair): a rotation by pi about +z is exactly the composition of
the x- and y-mirrors, i.e. the existing IMA mask 3.  The new rotation path must therefore reproduce the
already-golden mirror path to round-off on the same mesh -- a direct kernel-vs-kernel check that needs no
new reference geometry.

Why the physical ring lives in validation_test/: the reduced-vs-full RING comparison needs two meshes and
a real solve; see validation_test/feec/test_hdiv_vim_cyclic_ring.py.
"""
import math

import numpy as np
import pytest

pytest.importorskip("ngsolve")
import radia as rad  # noqa: E402
import radia.vim as vim  # noqa: E402
import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

A = 0.01
MU_R = 50.0
H0 = 1.0e4
HCF = ng.CoefficientFunction((0.0, 0.0, H0))


def _offset_box(n=2):
    """A small hex block OFF the z axis, so a rotation about z genuinely moves it."""
    return MakeStructured3DMesh(
        hexes=True, nx=n, ny=n, nz=n,
        mapping=lambda X, Y, Z: (1.2*A + A*X, -0.5*A + A*Y, -0.5*A + A*Z))


def _mean_magnetization(image_masks, image_signs, image_rot_angle, *, order=1):
    rad.UtiDelAll()
    with ng.TaskManager():
        mesh = _offset_box()
        fes = ng.HDiv(mesh, order=int(order))
        B, G, M_mass = vim.ChargeGram(
            fes, image_masks=image_masks, image_signs=image_signs,
            image_rot_angle=image_rot_angle)
        # The Gram entries are what the image path changes; probe them through the demag operator
        # applied to a fixed unit charge vector so the comparison needs no solve.
        rng = np.random.default_rng(20260811)
        x = rng.standard_normal(B.shape[1])
        return np.asarray(G.apply_configured_demag(x))


@pytest.mark.parametrize("order", (1, 2))
def test_pi_rotation_matches_mirror_pair(order):
    """R_pi about +z == the x-mirror composed with the y-mirror (IMA mask 3), to round-off.

    This pins the cyclic rotation kernel to the already-golden mirror kernel: same mesh, same sign,
    two different code paths through ImageEvalPoint.
    """
    mirror = _mean_magnetization([3], [1.0], [], order=order)
    rotate = _mean_magnetization([0], [1.0], [math.pi], order=order)
    scale = float(np.linalg.norm(mirror))
    assert scale > 0.0
    rel = float(np.linalg.norm(rotate - mirror) / scale)
    assert rel < 1e-10, f"pi-rotation image != mask-3 mirror image (rel {rel:.3e})"


def test_pi_rotation_matches_mirror_pair_negative_sign():
    """Same anchor with the antisymmetric sign -- the sign rides the image unchanged."""
    mirror = _mean_magnetization([3], [-1.0], [])
    rotate = _mean_magnetization([0], [-1.0], [math.pi])
    rel = float(np.linalg.norm(rotate - mirror) / float(np.linalg.norm(mirror)))
    assert rel < 1e-10, f"antisymmetric pi-rotation != mask-3 mirror (rel {rel:.3e})"


def test_rotation_changes_the_gram():
    """Guard against a no-op: a non-involutive rotation must NOT equal the direct Gram."""
    direct = _mean_magnetization([], [], [])
    quarter = _mean_magnetization(
        [0, 0, 0], [1.0, 1.0, 1.0],
        [0.5*math.pi, math.pi, 1.5*math.pi])
    rel = float(np.linalg.norm(quarter - direct) / float(np.linalg.norm(direct)))
    assert rel > 1e-3, "4-fold cyclic images left the Gram unchanged -- the fold is a no-op"


def test_unclosed_rotation_list_is_rejected():
    """The Gram symmetrizes each image with its transpose, so a list that is not closed under
    inversion would be silently symmetrized into different physics.  Fail loud instead."""
    with pytest.raises(Exception) as excinfo:
        _mean_magnetization([0], [1.0], [0.5*math.pi])       # +90 deg without its -90 deg partner
    assert "inver" in str(excinfo.value).lower()


def test_identity_image_is_rejected():
    """mask 0 with angle 0 is the identity -- it would silently double the direct term."""
    with pytest.raises(Exception) as excinfo:
        _mean_magnetization([0], [1.0], [0.0])
    assert "identity" in str(excinfo.value).lower()


def test_length_mismatch_is_rejected():
    with pytest.raises(Exception):
        _mean_magnetization([0, 0], [1.0, 1.0], [math.pi])


def test_solve_cyclic_alternating_requires_even_pole_count():
    rad.UtiDelAll()
    with ng.TaskManager():
        mesh = _offset_box()
        with pytest.raises(ValueError) as excinfo:
            vim.Solve(mesh, mu_r=MU_R, H_ext=HCF, order=1,
                      image_cyclic=5, image_cyclic_alternating=True)
    assert "EVEN" in str(excinfo.value)
