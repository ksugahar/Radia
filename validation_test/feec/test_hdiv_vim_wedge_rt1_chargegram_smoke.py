"""Smoke validation for the WEDGE RT1 charge-Gram H-matrix overload.

This mirrors the HEX low-level gate: Python passes the prism Q2 lattice,
monomial charge metadata, quadrature tables, and optional IMA image masks
directly to ``_ChargeGramHMatrix``.  The smallest one-prism case keeps the
pybind signature, entry oracle, H-matrix matvec, and image-fold path from
drifting independently.
"""

from __future__ import annotations

import math

import pytest

import radia._radia_pybind as _rp


def _unit_wedge_q2_nodes() -> list[float]:
    tri6 = [(1.0, 0.0), (0.0, 1.0), (0.0, 0.0),
            (0.5, 0.5), (0.0, 0.5), (0.5, 0.0)]
    nodes: list[float] = []
    for iz in range(3):
        z = iz / 2.0
        for x, y in tri6:
            nodes.extend([x, y, z])
    return nodes


def _make_wedge_rt1_gram(*, build: bool, image_masks=None, image_signs=None):
    return _rp._ChargeGramHMatrix(
        wedge_cell_nodes=_unit_wedge_q2_nodes(),
        face_nodes=[],
        face_type=[],
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
        near_grade=0.6,
        far_inner_factor=1.5,
        image_masks=list(image_masks or []),
        image_signs=list(image_signs or []),
        eps=1e-6,
        leaf=4,
        eta=2.0,
        build=build,
    )


def test_wedge_rt1_entry_oracle_is_positive():
    gram = _make_wedge_rt1_gram(build=False)

    entry = gram.entry(0, 0)

    assert math.isfinite(entry)
    assert entry > 0.0


def test_wedge_rt1_hmatrix_matvec_matches_entry_oracle():
    oracle = _make_wedge_rt1_gram(build=False)
    gram = _make_wedge_rt1_gram(build=True)

    entry = oracle.entry(0, 0)
    y = gram.matvec_sym([1.0])

    assert gram.ndof() == 1
    assert len(y) == 1
    assert y[0] == entry


def test_wedge_rt1_image_fold_increases_positive_entry_and_matvec():
    """A +x IMA mirror adds a finite positive image interaction to the direct one-prism entry."""
    direct = _make_wedge_rt1_gram(build=False)
    image_oracle = _make_wedge_rt1_gram(build=False, image_masks=[1], image_signs=[1.0])
    image = _make_wedge_rt1_gram(build=True, image_masks=[1], image_signs=[1.0])

    direct_entry = direct.entry(0, 0)
    image_entry = image_oracle.entry(0, 0)
    y = image.matvec_sym([1.0])

    assert image_entry > direct_entry
    assert len(y) == 1
    assert y[0] == image_entry


def test_wedge_rt1_image_masks_and_signs_must_match():
    with pytest.raises(ValueError, match="image_masks.*image_signs|image_signs.*image_masks"):
        _make_wedge_rt1_gram(build=False, image_masks=[1], image_signs=[])
