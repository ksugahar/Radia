"""Fast regression for the flat BDM2 analytic host-block ChargeGram path."""

import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


rp = pytest.importorskip("radia._radia_pybind")
ROOT = Path(__file__).resolve().parents[1]


def _g01(order):
    x, w = np.polynomial.legendre.leggauss(order)
    return 0.5 * (x + 1.0), 0.5 * w


def _tet_ref(order):
    x, w = _g01(order)
    points, weights = [], []
    for a, wa in zip(x, w):
        for b, wb in zip(x, w):
            for c, wc in zip(x, w):
                points.append((a, b * (1 - a), c * (1 - a) * (1 - b)))
                weights.append(wa * wb * wc * (1 - a) ** 2 * (1 - b))
    return np.asarray(points), np.asarray(weights)


def _tri_ref(order):
    x, w = _g01(order)
    points, weights = [], []
    for u, wu in zip(x, w):
        for v, wv in zip(x, w):
            points.append((u, v * (1 - u)))
            weights.append(wu * wv * (1 - u))
    return np.asarray(points), np.asarray(weights)


def _oracle(
    extra_unsupported_charge=False, *, image_masks=(), image_signs=(),
    image_rot_angle=(),
):
    vertices = np.asarray(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.2, 1.0, 0.0], [0.3, 0.4, 0.9]]
    )
    faces = np.asarray(
        [vertices[[0, 2, 1]], vertices[[0, 1, 3]], vertices[[1, 2, 3]], vertices[[2, 0, 3]]]
    )
    cell_exponents = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    face_exponents = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (2, 0, 0), (1, 1, 0), (0, 2, 0)]
    host = [0] * len(cell_exponents)
    kind = [0] * len(cell_exponents)
    exponents = list(cell_exponents)
    for face in range(4):
        host.extend([face] * len(face_exponents))
        kind.extend([1] * len(face_exponents))
        exponents.extend(face_exponents)
    if extra_unsupported_charge:
        # A degree-three face charge makes the constructor retain the scalar fallback for the whole oracle.
        # It is not sampled below; its only job is to provide an independent old-path reference.
        host.append(0)
        kind.append(1)
        exponents.append((3, 0, 0))
    tet_points, tet_weights = _tet_ref(3)
    tri_points, tri_weights = _tri_ref(3)
    gram = rp._ChargeGramHMatrix(
        cell_verts=vertices.ravel().tolist(),
        face_verts=faces.ravel().tolist(),
        n_el=1,
        charge_host=host,
        charge_kind=kind,
        charge_expo=np.asarray(exponents, dtype=np.int32).ravel().tolist(),
        ref_tet_pts=tet_points.ravel().tolist(),
        ref_tet_w=tet_weights.tolist(),
        ref_tri_pts=tri_points.ravel().tolist(),
        ref_tri_w=tri_weights.tolist(),
        image_masks=np.asarray(image_masks, dtype=np.int32),
        image_signs=np.asarray(image_signs, dtype=float),
        build=False,
    )
    if image_rot_angle:
        gram.set_image_rotations(np.asarray(image_rot_angle, dtype=float))
    return gram


def test_rt2_host_block_matches_scalar_analytic_entries():
    """Cell/face block contraction must preserve every old scalar entry to roundoff."""
    block = _oracle()
    scalar = _oracle(extra_unsupported_charge=True)
    n = 28
    block_matrix = np.asarray([[block.entry(i, j) for j in range(n)] for i in range(n)])
    scalar_matrix = np.asarray([[scalar.entry(i, j) for j in range(n)] for i in range(n)])
    scale = np.max(np.abs(scalar_matrix))
    assert np.max(np.abs(block_matrix - scalar_matrix)) <= 2e-15 * scale
    assert np.max(np.abs(block_matrix - block_matrix.T)) <= 2e-15 * scale


_CYCLIC_IMAGE_CHILD = r"""
import json
import math
import runpy
import sys

import numpy as np

scope = runpy.run_path(sys.argv[1])
gram = scope["_oracle"](
    image_masks=(0, 0, 0),
    image_signs=(-1.0, 1.0, -1.0),
    image_rot_angle=(0.5 * math.pi, math.pi, 1.5 * math.pi),
)
n = 28
matrix = np.asarray([[gram.entry(i, j) for j in range(n)] for i in range(n)])
print("RESULT " + json.dumps({"matrix": matrix.tolist(), "stats": dict(gram.stats())}))
"""


def _cyclic_image_matrix(*, disable_image_block):
    env = dict(os.environ)
    for name in tuple(env):
        if name.startswith("RADIA_HDIV_") or name == "RADIA_HACAPK_MATVEC_MKL_THREADS":
            env.pop(name)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(ROOT / "src"), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    if disable_image_block:
        env["RADIA_HDIV_DISABLE_HO_IMAGE_BLOCK"] = "1"
    else:
        env.pop("RADIA_HDIV_DISABLE_HO_IMAGE_BLOCK", None)
    proc = subprocess.run(
        [sys.executable, "-c", _CYCLIC_IMAGE_CHILD, str(Path(__file__).resolve())],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=180,
    )
    for line in proc.stdout.splitlines():
        if line.startswith("RESULT "):
            return json.loads(line[7:])
    raise AssertionError(
        "cyclic-image child produced no result:\n%s\n%s"
        % (proc.stdout[-2000:], proc.stderr[-2000:])
    )


def test_rt2_cyclic_image_host_block_matches_scalar_entries():
    """Vectorized cyclic-image blocks must preserve the scalar image oracle."""
    blocked_result = _cyclic_image_matrix(disable_image_block=False)
    scalar_result = _cyclic_image_matrix(disable_image_block=True)
    blocked = np.asarray(blocked_result["matrix"])
    scalar = np.asarray(scalar_result["matrix"])

    scale = np.max(np.abs(scalar))
    assert scale > 0.0
    assert np.max(np.abs(blocked - scalar)) <= 4e-15 * scale
    assert np.max(np.abs(blocked - blocked.T)) <= 2e-15 * scale
    assert blocked_result["stats"]["ho_image_block_enabled"] == 1.0
    assert blocked_result["stats"]["release_claim_eligible"] == 1.0
    assert scalar_result["stats"]["ho_image_block_enabled"] == 0.0
    assert scalar_result["stats"]["nonproduction_numerical_override_active"] == 1.0
    assert scalar_result["stats"]["release_claim_eligible"] == 0.0
