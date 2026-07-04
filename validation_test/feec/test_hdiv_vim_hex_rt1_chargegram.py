"""Smoke validation for the HEX RT1 charge-Gram H-matrix overload.

The C++ overload is intentionally low-level: Python passes the Q2 hex lattice,
Q1 charge metadata, and quadrature tables directly to ``_ChargeGramHMatrix``.
This test keeps the smallest possible case live so the pybind signature,
entry oracle, and H-matrix matvec cannot drift independently.
"""

from __future__ import annotations

import math

import pytest

import radia._radia_pybind as _rp


def _unit_hex_q2_nodes() -> list[float]:
    nodes: list[float] = []
    for iz in range(3):
        for iy in range(3):
            for ix in range(3):
                nodes.extend([ix / 2.0, iy / 2.0, iz / 2.0])
    return nodes


def _make_hex_rt1_gram(*, build: bool, image_masks=None, image_signs=None):
    return _rp._ChargeGramHMatrix(
        hex_cell_nodes=_unit_hex_q2_nodes(),
        quad_face_nodes=[],
        n_el=1,
        n_bf=0,
        charge_host=[0],
        charge_kind=[0],
        charge_expo=[0, 0, 0],
        sym_tet_pts=[0.25, 0.25, 0.25],
        sym_tet_w=[1.0 / 6.0],
        sym_tri_pts=[1.0 / 3.0, 1.0 / 3.0],
        sym_tri_w=[0.5],
        gl_out=[0.5],
        gw_out=[1.0],
        gl_in=[0.5],
        gw_in=[1.0],
        far_tet_pts=[0.25, 0.25, 0.25],
        far_tet_w=[1.0 / 6.0],
        far_tri_pts=[1.0 / 3.0, 1.0 / 3.0],
        far_tri_w=[0.5],
        near_grade=1.5,
        far_inner_factor=4.0,
        image_masks=list(image_masks or []),
        image_signs=list(image_signs or []),
        eps=1e-6,
        leaf=4,
        eta=2.0,
        build=build,
    )


def test_hex_rt1_entry_oracle_is_positive():
    gram = _make_hex_rt1_gram(build=False)

    entry = gram.entry(0, 0)

    assert math.isfinite(entry)
    assert entry > 0.0


def test_hex_rt1_hmatrix_matvec_matches_entry_oracle():
    oracle = _make_hex_rt1_gram(build=False)
    gram = _make_hex_rt1_gram(build=True)

    entry = oracle.entry(0, 0)
    y = gram.matvec_sym([1.0])

    assert gram.ndof() == 1
    assert len(y) == 1
    assert y[0] == entry


def test_hex_rt1_image_fold_increases_positive_entry_and_matvec():
    """IMA image folding is part of the HEX RT1 Gram constructor, not a Python-side postprocess.

    On a unit hex with one positive volume-charge basis, a +x mirror adds a finite positive image
    interaction on top of the direct self entry.  The built H-matrix matvec must see exactly the same folded
    entry as the entry oracle.
    """
    direct = _make_hex_rt1_gram(build=False)
    image_oracle = _make_hex_rt1_gram(build=False, image_masks=[1], image_signs=[1.0])
    image = _make_hex_rt1_gram(build=True, image_masks=[1], image_signs=[1.0])

    direct_entry = direct.entry(0, 0)
    image_entry = image_oracle.entry(0, 0)
    y = image.matvec_sym([1.0])

    assert image_entry > direct_entry
    assert len(y) == 1
    assert y[0] == image_entry


def test_hex_rt1_image_masks_and_signs_must_match():
    with pytest.raises(ValueError, match="image_masks.*image_signs|image_signs.*image_masks"):
        _make_hex_rt1_gram(build=False, image_masks=[1], image_signs=[])
