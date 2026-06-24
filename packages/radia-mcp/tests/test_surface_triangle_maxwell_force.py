"""Readable surface-triangle Maxwell traction gates."""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    MU0,
    air_gap_maxwell_pressure,
    maxwell_traction_summary,
    surface_triangle_maxwell_traction_summary,
)


def _sum_vectors(rows):
    n = len(rows[0])
    return [sum(row[i] for row in rows) for i in range(n)]


def test_surface_triangle_normal_field_distributes_p1_force_load():
    tri = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    row = surface_triangle_maxwell_traction_summary(tri, (0.0, 0.0, 1.0))
    pressure = air_gap_maxwell_pressure(1.0)

    assert row["area"] == pytest.approx(0.5)
    assert row["unit_normal"] == pytest.approx([0.0, 0.0, 1.0])
    assert row["traction_Pa"] == pytest.approx([0.0, 0.0, pressure])
    assert row["integrated_force_N"] == pytest.approx([0.0, 0.0, 0.5 * pressure])
    assert _sum_vectors(row["nodal_force_loads_N"]) == pytest.approx(row["integrated_force_N"])
    assert row["nodal_force_loads_N"][0] == pytest.approx([0.0, 0.0, pressure / 6.0])


def test_reversing_triangle_orientation_reverses_integrated_force():
    tri = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    rev = list(reversed(tri))
    row = surface_triangle_maxwell_traction_summary(tri, (0.0, 0.0, 1.0))
    rev_row = surface_triangle_maxwell_traction_summary(rev, (0.0, 0.0, 1.0))

    assert rev_row["area"] == pytest.approx(row["area"])
    assert rev_row["unit_normal"] == pytest.approx([0.0, 0.0, -1.0])
    assert rev_row["integrated_force_N"] == pytest.approx(
        [-value for value in row["integrated_force_N"]]
    )


def test_surface_triangle_oblique_field_matches_patch_traction_summary():
    tri = [(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    B = (3.0, 4.0, 0.0)
    row = surface_triangle_maxwell_traction_summary(tri, B)
    patch = maxwell_traction_summary(B, (1.0, 0.0, 0.0), area_m2=row["area"])

    assert row["unit_normal"] == pytest.approx([1.0, 0.0, 0.0])
    assert row["traction_Pa"] == pytest.approx(patch["traction_Pa"])
    assert row["integrated_force_N"] == pytest.approx(patch["force_N"])
    assert row["traction_Pa"] == pytest.approx([-3.5 / MU0, 12.0 / MU0, 0.0])


def test_surface_triangle_maxwell_traction_rejects_bad_inputs():
    with pytest.raises(ValueError):
        surface_triangle_maxwell_traction_summary(
            [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (2.0, 2.0, 2.0)],
            (0.0, 0.0, 1.0),
        )
    with pytest.raises(ValueError):
        surface_triangle_maxwell_traction_summary(
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
            (0.0, 1.0),
        )


if __name__ == "__main__":
    test_surface_triangle_normal_field_distributes_p1_force_load()
    test_reversing_triangle_orientation_reverses_integrated_force()
    test_surface_triangle_oblique_field_matches_patch_traction_summary()
    test_surface_triangle_maxwell_traction_rejects_bad_inputs()
    print("[OK] surface-triangle Maxwell traction helpers validated.")
