"""The accepted beak-fin delivery gate must rest on a tracked fixture.

The gate closed on a 48 mm geometry, but that STEP lived only in a scratch
directory, so the acceptance could not be reproduced from the repository
alone.  These checks pin the tracked fixture to the geometry the generator
produces and to the basis the metrics file declares, so the two cannot drift
apart silently.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "coil_from_cad" / "fixtures" / "beak_fin_48mm.step"
METRICS = (ROOT / "validation_test" / "induction_heating"
           / "beak_fin_delivery_metrics_150kHz.json")
ACCEPTED_LENGTH_MM = 48.0


def _make_beak_fin(length_mm):
    sys.path.insert(0, str(ROOT / "validation_test" / "induction_heating"))
    from make_beak_fin_step import make_beak_fin

    return make_beak_fin(length_mm)


def test_accepted_fixture_is_tracked():
    assert FIXTURE.is_file(), (
        "the accepted 48 mm delivery-gate basis is not in the repository; "
        "the acceptance cannot be reproduced without it")


def test_metrics_declare_the_fixture_length():
    doc = json.loads(METRICS.read_text(encoding="utf-8"))
    assert doc["acceptance"]["accepted"] is True
    assert doc["length_resolved_2026_09_21"]["geometry_length_mm"] == (
        ACCEPTED_LENGTH_MM)
    assert "48 mm" in doc["acceptance"]["basis"]


def test_fixture_matches_the_generator():
    """Byte equality is not expected -- a STEP writer stamps its own header --
    so compare the solid the two produce."""
    build123d = pytest.importorskip("build123d")

    tracked = build123d.import_step(str(FIXTURE)).solids()
    assert len(tracked) == 1, "the beak fin is one solid by construction"
    tracked = tracked[0]
    generated = _make_beak_fin(ACCEPTED_LENGTH_MM)

    assert tracked.volume == pytest.approx(generated.volume, rel=1e-9)
    assert tracked.area == pytest.approx(generated.area, rel=1e-9)
    got, want = tracked.bounding_box(), generated.bounding_box()
    for axis in ("X", "Y", "Z"):
        assert getattr(got.min, axis) == pytest.approx(
            getattr(want.min, axis), abs=1e-9)
        assert getattr(got.max, axis) == pytest.approx(
            getattr(want.max, axis), abs=1e-9)


def test_fixture_is_the_declared_length():
    build123d = pytest.importorskip("build123d")

    box = build123d.import_step(str(FIXTURE)).solids()[0].bounding_box()
    assert box.max.Z - box.min.Z == pytest.approx(ACCEPTED_LENGTH_MM, abs=1e-9)
