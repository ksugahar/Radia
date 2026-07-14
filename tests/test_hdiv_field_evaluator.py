"""Fast contracts for the persistent C++ HDiv exterior-field evaluator."""

import numpy as np
import pytest

rp = pytest.importorskip("radia._radia_pybind")


def _empty_images():
    return np.empty(0, dtype=np.int32), np.empty(0, dtype=float)


def _cloud(xyz, strength, *, masks=None, signs=None, theta=0.05,
           tree_rtol=1.0e-5):
    if masks is None:
        masks, signs = _empty_images()
    return rp._HDivFieldEvaluator.from_cloud(
        np.ascontiguousarray(xyz, dtype=float),
        np.ascontiguousarray(strength, dtype=float),
        np.ascontiguousarray(masks, dtype=np.int32),
        np.ascontiguousarray(signs, dtype=float),
        16, theta, 2, 1, tree_rtol, 8,
    )


def test_numpy_direct_buffer_and_ima_one_call_match_reference():
    rng = np.random.default_rng(12)
    xyz = rng.uniform(-0.5, 0.5, (80, 3))
    strength = rng.normal(size=80)
    obs = rng.uniform(1.0, 2.0, (12, 3))
    plain = _cloud(xyz, strength)
    direct = np.asarray(plain.field(obs, "direct"))

    delta = obs[:, None, :] - xyz[None, :, :]
    r2 = np.einsum("nqi,nqi->nq", delta, delta)
    reference = np.einsum("q,nqi,nq->ni", strength, delta, r2**-1.5)
    assert np.allclose(direct, reference, rtol=2e-14, atol=2e-14)

    masks = np.array([1, 4, 5], dtype=np.int32)
    signs = np.array([1.0, -1.0, -1.0])
    images = _cloud(xyz, strength, masks=masks, signs=signs)
    got = np.asarray(images.field(obs, "direct"))
    expected = direct.copy()
    for mask, sign in zip(masks, signs):
        axes = [axis for axis in range(3) if mask & (1 << axis)]
        reflected = obs.copy()
        reflected[:, axes] *= -1.0
        contribution = np.asarray(plain.field(reflected, "direct"))
        contribution[:, axes] *= -1.0
        expected += sign*contribution
    assert np.allclose(got, expected, rtol=2e-14, atol=2e-14)
    assert images.stats()["image_count"] == 3
    assert images.candidate_algorithm_for(10_000_000) == "direct"


def test_quadrupole_tree_matches_direct_for_far_cloud():
    rng = np.random.default_rng(8)
    xyz = rng.uniform(-0.5, 0.5, (2000, 3))
    strength = rng.uniform(0.1, 1.0, 2000)
    obs = rng.uniform(2.0, 3.0, (200, 3))*rng.choice([-1.0, 1.0], (200, 3))
    evaluator = _cloud(xyz, strength)
    direct = np.asarray(evaluator.field(obs, "direct"))
    tree = np.asarray(evaluator.field(obs, "tree"))
    relative = np.linalg.norm(tree-direct, axis=1)/np.linalg.norm(direct, axis=1)
    assert relative.max() < 5.0e-6
    assert evaluator.stats()["tree_nodes"] > 1


def test_auto_rejects_inaccurate_signed_tree():
    rng = np.random.default_rng(2)
    xyz = rng.uniform(-0.5, 0.5, (1200, 3))
    strength = rng.normal(size=1200)
    strength -= strength.mean()
    obs = rng.uniform(1.5, 2.5, (64, 3))*rng.choice([-1.0, 1.0], (64, 3))
    # Deliberately loose geometry opening plus strict probe tolerance: forced
    # tree is inaccurate for this cancellation-heavy cloud, so auto must use
    # the exact direct path.
    evaluator = _cloud(xyz, strength, theta=0.2, tree_rtol=1.0e-10)
    direct = np.asarray(evaluator.field(obs, "direct"))
    tree = np.asarray(evaluator.field(obs, "tree"))
    automatic = np.asarray(evaluator.field(obs, "auto"))
    assert np.linalg.norm(tree-direct) > 1.0e-8*np.linalg.norm(direct)
    assert np.array_equal(automatic, direct)


def test_field_input_shape_is_checked():
    evaluator = _cloud(np.array([[0.0, 0.0, 0.0]]), np.array([1.0]))
    with pytest.raises(RuntimeError, match="shape"):
        evaluator.field(np.zeros(3), "direct")
