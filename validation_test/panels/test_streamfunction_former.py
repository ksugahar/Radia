"""Small contracts for the stream-function printable-former STEP output."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest


REPO = Path(__file__).resolve().parents[2]
PANELS = REPO / "src" / "radia" / "panels"
sys.path.insert(0, str(PANELS))

import calc_streamfunction as csf  # noqa: E402


def test_decimate_polyline_keeps_endpoints_and_reduces_dense_spine():
    poly = np.array([[0.0, 0.0, 0.0], [0.001, 0.0, 0.0], [0.003, 0.0, 0.0], [0.010, 0.0, 0.0]])

    decimated = csf._decimate_polyline(poly, 0.004)

    assert np.allclose(decimated[0], poly[0])
    assert np.allclose(decimated[-1], poly[-1])
    assert len(decimated) < len(poly)


def test_write_former_step_smoke(tmp_path):
    pytest.importorskip("netgen.occ")
    out = tmp_path / "former.step"
    chain = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.02, 0.0, 0.001],
            [0.03, 0.01, 0.0],
        ],
        float,
    )

    info = csf._write_former_step(
        chain,
        wire_diam=1.0e-3,
        clearance=2.0e-4,
        margin=2.0e-3,
        wall=1.0e-3,
        decimate=0.0,
        filename=str(out),
    )

    assert out.exists()
    assert out.stat().st_size > 0
    assert info["channel_segments"] == 2
    assert info["channel_volume_m3"] > 0.0
    assert info["former_volume_m3"] < info["plate_volume_m3"]
