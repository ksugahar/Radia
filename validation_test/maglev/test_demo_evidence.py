import importlib.util
import json
import math
from pathlib import Path


EVIDENCE = Path(__file__).resolve().parent / "demos"
DOCS = Path(__file__).resolve().parents[2] / "docs" / "maglev" / "demos"
EXPECTED = {
    "cube_alpha_sweep_results.json",
    "cube_alpha_tensor_results.json",
    "cuboid_vector_bulk_results.json",
    "physical_tensor_rom_fem.json",
    "physical_tensor_rom_sphere.json",
    "rotating_magnet_eddy_results.json",
    "ellipsoid/ellipsoid_alpha_omega_axisym_results.json",
    "ellipsoid/ellipsoid_alpha_tensor_3d_results.json",
    "ellipsoid/ellipsoid_alpha_tensor_results.json",
    "sphere/coil_maglev_equilibrium_results.json",
    "sphere/coil_sphere_eddy_force_results.json",
    "sphere/maglev_sphere_force_results.json",
    "team28/team28_cln_sweep_results.json",
}


def _load(name):
    return json.loads((EVIDENCE / name).read_text(encoding="utf-8"))


def _numbers(value):
    if isinstance(value, dict):
        for child in value.values():
            yield from _numbers(child)
    elif isinstance(value, list):
        for child in value:
            yield from _numbers(child)
    elif value is None:
        return
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        yield float(value)


def test_maglev_demo_evidence_inventory_is_complete_and_finite():
    actual = {path.relative_to(EVIDENCE).as_posix() for path in EVIDENCE.rglob("*.json")}
    assert actual == EXPECTED
    for name in EXPECTED:
        assert all(math.isfinite(value) for value in _numbers(_load(name))), name


def test_docs_json_routes_to_validation_and_copied_runs_stay_local(tmp_path):
    spec = importlib.util.spec_from_file_location(
        "maglev_validation_output", DOCS / "_validation_output.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    assert module.validation_output("root.json", DOCS) == EVIDENCE / "root.json"
    assert module.validation_output("sphere.json", DOCS / "sphere") == (
        EVIDENCE / "sphere" / "sphere.json"
    )
    assert module.validation_output("copied.json", tmp_path) == tmp_path / "copied.json"


def test_mixed_galerkin_and_rom_claims_are_preserved():
    sweep = _load("cube_alpha_sweep_results.json")
    assert len(sweep["f_hz"]) == len(sweep["re_alpha_over_V"]) == 73
    assert sweep["re_alpha_over_V"][-1] > 0.99

    tensor = _load("cube_alpha_tensor_results.json")
    assert tensor["n_foster"] > 0
    assert tensor["mimo_transfer_worst_relerr"] < 0.01

    bulk = _load("cuboid_vector_bulk_results.json")["runs"]["h_0.18mm"]["leading_tau"]
    assert bulk["z"]["tau_us"] > bulk["y"]["tau_us"] > bulk["x"]["tau_us"]
    assert max(abs(bulk[axis]["err_pct"]) for axis in "xyz") < 3.0

    sphere = _load("physical_tensor_rom_sphere.json")
    assert sphere["rom"]["band_fit_relerr"] < 0.001
    fem = _load("physical_tensor_rom_fem.json")["rom"]
    assert all(fem[axis]["band_fit_relerr"] < 0.05 for axis in "xyz")
    assert fem["z"]["dominant_tau_us"][0] > fem["x"]["dominant_tau_us"][0]


def test_shape_anisotropy_and_force_claims_are_preserved():
    axisym = _load("ellipsoid/ellipsoid_alpha_omega_axisym_results.json")
    assert axisym["sphere_validation_worst_relerr"] < 0.02

    analytic = _load("ellipsoid/ellipsoid_alpha_tensor_results.json")
    assert abs(sum(analytic["demag_N"]) - 1.0) < 1.0e-10
    assert analytic["anisotropy_kappaz_over_kappax"] > 2.0

    full3d = _load("ellipsoid/ellipsoid_alpha_tensor_3d_results.json")
    assert full3d["sphere_worst_relerr"] < 0.04
    assert full3d["sphere_isotropy_worst"] < 0.001
    assert full3d["reproduces_analytic_ordering"] is True

    equilibrium = _load("sphere/coil_maglev_equilibrium_results.json")
    assert equilibrium["stability"]["stable"] is True
    assert abs(equilibrium["equilibrium"]["residual_percent"]) < 0.01

    force = _load("sphere/coil_sphere_eddy_force_results.json")
    assert force["key_check"]["ratio_moves_toward_one_as_aL_shrinks"] is True
    assert force["key_check"]["abs_dev_1mm"] < force["key_check"]["abs_dev_5mm"]


def test_dynamic_and_team28_claims_are_preserved():
    moving = _load("rotating_magnet_eddy_results.json")["rows"]
    assert len(moving) >= 4
    assert moving[0]["source_only_F_err"] < 0.001
    assert moving[-1]["source_only_F_err"] > 0.1
    assert max(row["cln_F_err"] for row in moving) < 0.003

    team28 = _load("team28/team28_cln_sweep_results.json")
    assert len(team28["dZ_mm"]) == len(team28["fz_cln_N"]) == 25
    assert team28["max_abs_cln_minus_lab_N"] < 0.001
    assert abs(
        team28["equilibrium_abs_height_mm"]["cln"]
        - team28["published_ref"]["levitation_height_mm"]
    ) < 0.6
