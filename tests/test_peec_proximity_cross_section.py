"""Cross-section contracts for the perimeter PEEC proximity iteration."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from radia.peec_proximity import solve_proximity_iterative  # noqa: E402


def test_explicit_rectangular_model_uses_the_supplied_perimeter():
    paths = [
        [((0.0, -1.0e-3, 0.0), (0.1, -1.0e-3, 0.0))],
        [((0.0, 1.0e-3, 0.0), (0.1, 1.0e-3, 0.0))],
    ]
    result = solve_proximity_iterative(
        R_f=np.diag([1.0e-3, 1.0e-3]),
        L_f=np.diag([1.0e-7, 1.0e-7]),
        filament_paths=paths,
        frequency=100.0e3,
        sigma=5.8e7,
        wire_radius_m=0.0,
        n_peri=2,
        max_iter=1,
        self_impedance_per_m=2.0e-3 + 1.0e-3j,
        dc_resistance_per_m=1.0e-3,
        perimeter_m=12.0e-3,
        internal_impedance_model="rectangular-dowell",
    )
    assert result["perimeter_m"] == pytest.approx(12.0e-3)
    assert result["internal_impedance_model"] == "rectangular-dowell"
    assert np.all(np.isfinite(result["Zs_fil"]))


def test_explicit_model_requires_a_physical_perimeter():
    with pytest.raises(ValueError, match="perimeter_m"):
        solve_proximity_iterative(
            R_f=np.diag([1.0e-3]),
            L_f=np.diag([1.0e-7]),
            filament_paths=[[
                ((0.0, 0.0, 0.0), (0.1, 0.0, 0.0)),
            ]],
            frequency=100.0e3,
            sigma=5.8e7,
            wire_radius_m=0.0,
            n_peri=1,
            self_impedance_per_m=2.0e-3 + 1.0e-3j,
            dc_resistance_per_m=1.0e-3,
            perimeter_m=0.0,
        )


@pytest.mark.parametrize(
    ("paths", "n_peri", "message"),
    [
        ([], 1, "one path per PEEC filament"),
        ([[((0.0, 0.0, 0.0), (0.1, 0.0, 0.0))]], 2,
         "must be an integer multiple of n_peri"),
    ],
)
def test_proximity_rejects_inconsistent_filament_bundle(paths, n_peri, message):
    with pytest.raises(ValueError, match=message):
        solve_proximity_iterative(
            R_f=np.diag([1.0e-3]),
            L_f=np.diag([1.0e-7]),
            filament_paths=paths,
            frequency=100.0e3,
            sigma=5.8e7,
            wire_radius_m=1.0e-3,
            n_peri=n_peri,
        )
