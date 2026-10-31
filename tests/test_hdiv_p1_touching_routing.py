"""Production-entry routing of affine P1 self/touching ChargeGram pairs.

Guards: the direct value, the image (IMA) value and the TET shape derivative
all take the endpoint-graded P1 outer rule for self/touching host pairs.  A
bypass of any one path brings back the fixed-rule truncation (about 0.5 % on a
face self entry, about 1.3 % on its derivative) or breaks the exact on-plane
image identity.
Assumes: the closed-form triangle potential below is independent of the
production kernels; the 1e-5 relative detector is a PROPOSED surrogate
threshold, not an adopted acceptance criterion.
"""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402
from radia.vim._vim import _prepare_gpu_charge_inputs, _rp  # noqa: E402

DETECTOR_REL = 1e-5


def _triangle_potential(V, x):
    """INT_T 1/|x-y| dS for x in the triangle plane (closed form)."""
    n = np.cross(V[1] - V[0], V[2] - V[0])
    n /= np.linalg.norm(n)
    total = 0.0
    for i in range(3):
        a, b, c = V[i], V[(i + 1) % 3], V[(i + 2) % 3]
        t = (b - a) / np.linalg.norm(b - a)
        m = np.cross(t, n)
        if np.dot(m, c - a) > 0:
            m = -m
        h = np.dot(a - x, m)
        sa, sb = np.dot(a - x, t), np.dot(b - x, t)
        total += h * np.log((sb + np.linalg.norm(b - x)) / (sa + np.linalg.norm(a - x)))
    return total


def _triangle_self_energy(V, order=48):
    q, w = np.polynomial.legendre.leggauss(order)
    q, w = (q + 1) / 2, w / 2
    e = V[1:] - V[0]
    area = np.linalg.norm(np.cross(e[0], e[1]))
    total = sum(wi * wj * (1 - qi) * _triangle_potential(V, V[0] + qi * e[0] + (1 - qi) * qj * e[1])
                for qi, wi in zip(q, w) for qj, wj in zip(q, w))
    return total * area / (4 * np.pi)


@pytest.fixture(scope="module")
def charge_inputs():
    mesh = MakeStructured3DMesh(hexes=False, nx=1, ny=1, nz=1)
    with ng.TaskManager():
        _, kw = _prepare_gpu_charge_inputs(ng.HDiv(mesh, order=1), backend="gpu")
    return kw


# Charges of the first cell and the faces sharing a vertex with it: self,
# touching and far pairs, small enough for a surrogate run.
def _subset(kw):
    cells = np.asarray(kw["cell_verts"]).reshape(-1, 4, 3)
    faces = np.asarray(kw["face_verts"]).reshape(-1, 3, 3)
    kind, host = np.asarray(kw["charge_kind"]), np.asarray(kw["charge_host"])
    corner = {tuple(v) for v in cells[0]}
    keep = [i for i in range(len(kind)) if (kind[i] == 0 and host[i] == 0) or
            (kind[i] == 1 and sum(tuple(v) in corner for v in faces[host[i]]) >= 2)]
    return np.array(keep)


def _gram(kw, **extra):
    gram = _rp._ChargeGramHMatrix(**{**kw, **extra})
    idx = _subset(kw)
    return gram, np.array([[gram.entry(i, j) for j in idx] for i in idx])


def _constant_faces(kw):
    kind = np.asarray(kw["charge_kind"])
    expo = np.asarray(kw["charge_expo"]).reshape(-1, 3)
    return [i for i in range(len(kind)) if kind[i] == 1 and not expo[i].any()]


def test_direct_face_self_entry_takes_graded_rule(charge_inputs):
    kw = charge_inputs
    gram = _rp._ChargeGramHMatrix(**kw)
    faces = np.asarray(kw["face_verts"]).reshape(-1, 3, 3)
    host = np.asarray(kw["charge_host"])
    i = _constant_faces(kw)[0]
    actual = gram.entry(i, i)
    expected = _triangle_self_energy(faces[host[i]])
    rel = abs(actual - expected) / expected
    assert rel < DETECTOR_REL, f"face self entry {actual!r} vs closed form {expected!r} (rel {rel:.3e})"


def test_on_plane_image_doubles_direct_entries(charge_inputs):
    kw = charge_inputs
    direct = _rp._ChargeGramHMatrix(**kw)
    image = _rp._ChargeGramHMatrix(**kw, image_masks=np.array([1], dtype=np.int32),
                                   image_signs=np.array([1.0]))
    faces = np.asarray(kw["face_verts"]).reshape(-1, 3, 3)
    kind, host = np.asarray(kw["charge_kind"]), np.asarray(kw["charge_host"])
    on_plane = [i for i in range(len(kind)) if kind[i] == 1 and np.all(faces[host[i], :, 0] == 0.0)]
    assert on_plane, "fixture lost its faces on the x=0 mirror plane"
    for i in on_plane:
        assert image.entry(i, i) == pytest.approx(2 * direct.entry(i, i), rel=1e-12, abs=0.0), (
            f"charge {i}: image {image.entry(i, i)!r} vs 2*direct {2 * direct.entry(i, i)!r}")


@pytest.mark.parametrize("images", [None, ([1], [1.0])], ids=["direct", "image"])
def test_shape_derivative_matches_scaling_and_central_difference(charge_inputs, images):
    kw = dict(charge_inputs)
    if images:
        kw.update(image_masks=np.array(images[0], dtype=np.int32), image_signs=np.array(images[1]))
    gram, G = _gram(kw)
    cells = np.asarray(kw["cell_verts"]).reshape(-1, 4, 3)
    faces = np.asarray(kw["face_verts"]).reshape(-1, 3, 3)
    idx = _subset(kw)
    kind = np.asarray(kw["charge_kind"])[idx]
    # Uniform scaling: a surface pair is homogeneous of degree 3, a volume pair 5.
    D = np.asarray(gram.tet_charge_gram_directional_derivative(
        np.ascontiguousarray(cells), np.ascontiguousarray(faces)))[np.ix_(idx, idx)]
    expected = np.where(kind == 1, 3.0, 5.0) * np.diag(G)
    np.testing.assert_allclose(np.diag(D), expected, rtol=1e-10, atol=0.0)

    # Nonlinear field tangential on the x=0 mirror plane (keeps the symmetry).
    def field(x):
        return np.stack([x[..., 0] * np.sin(x[..., 0] + 2 * x[..., 1]),
                         np.cos(x[..., 1]) - x[..., 2] ** 2, x[..., 0] * x[..., 2] + 0.3], -1)

    vc, vf = field(cells), field(faces)
    D = np.asarray(gram.tet_charge_gram_directional_derivative(
        np.ascontiguousarray(vc), np.ascontiguousarray(vf)))[np.ix_(idx, idx)]
    h = 1e-5
    plus = _gram(kw, cell_verts=np.ascontiguousarray((cells + h * vc).ravel()),
                 face_verts=np.ascontiguousarray((faces + h * vf).ravel()))[1]
    minus = _gram(kw, cell_verts=np.ascontiguousarray((cells - h * vc).ravel()),
                  face_verts=np.ascontiguousarray((faces - h * vf).ravel()))[1]
    fd = (plus - minus) / (2 * h)
    err = np.max(np.abs(D - fd)) / np.max(np.abs(fd))
    assert err < DETECTOR_REL, f"derivative vs central difference: max rel {err:.3e}"


def test_image_value_and_derivative_are_continuous_off_the_plane(charge_inputs):
    # Translate the mesh off the mirror plane by d: near images then stay on one
    # rule, so the translation derivative matches a central difference taken
    # anywhere in d (a geometric rule switch shows up here as a jump).
    images = dict(image_masks=np.array([1], dtype=np.int32), image_signs=np.array([1.0]))
    cells = np.asarray(charge_inputs["cell_verts"]).reshape(-1, 4, 3)
    faces = np.asarray(charge_inputs["face_verts"]).reshape(-1, 3, 3)
    shift = np.array([1.0, 0.0, 0.0])

    def moved(d):
        return dict(charge_inputs, **images, cell_verts=np.ascontiguousarray((cells + d * shift).ravel()),
                    face_verts=np.ascontiguousarray((faces + d * shift).ravel()))

    d, h = 0.37267799625, 1e-5
    idx = _subset(charge_inputs)
    gram, _ = _gram(moved(d))
    D = np.asarray(gram.tet_charge_gram_directional_derivative(
        np.ascontiguousarray(np.broadcast_to(shift, cells.shape)),
        np.ascontiguousarray(np.broadcast_to(shift, faces.shape))))[np.ix_(idx, idx)]
    fd = (_gram(moved(d + h))[1] - _gram(moved(d - h))[1]) / (2 * h)
    err = np.max(np.abs(D - fd)) / np.max(np.abs(fd))
    assert err < DETECTOR_REL, f"translation derivative vs central difference: max rel {err:.3e}"
