import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
EXPECTED = {
    "bending_endpack_saturation_opt.json",
    "bidirectional_coordinate_transform_2d.json",
    "clebsch_dipole_saturation_3d_throat.json",
    "clebsch_pole_shape_optimization_2d.json",
    "combined_function_frenet_sweep.json",
    "endpack_cobake.json",
    "endpack_cobake_loft.json",
    "endpack_spectrometer_saturation.json",
    "endpack_two_plane.json",
    "ffag_sector_two_plane.json",
    "leaf_coupling_perturbation_3d.json",
    "leaf_coupling_perturbation_3d_sweep.json",
    "scaling_ffag_sector_saturation.json",
    "twist_rate_leaf_coupling.json",
    "twisting_quadrupole_pole.json",
    "weakform_pullback_kata.json",
}


def _reject_non_finite(token):
    raise ValueError(f"non-finite JSON token: {token}")


def test_demo_record_inventory_is_complete_and_strict_json():
    actual = {path.name for path in HERE.glob("*.json")}
    assert actual == EXPECTED

    for path in sorted(HERE.glob("*.json")):
        data = json.loads(
            path.read_text(encoding="utf-8"),
            parse_constant=_reject_non_finite,
        )
        assert isinstance(data, dict), path
        assert data, path
