"""The gap-family report: synthetic three-engine results, no solver."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

path = Path(__file__).resolve().parents[1] / "validation_test/c_type_three_engine/analyze_gap_family.py"
spec = importlib.util.spec_from_file_location("c_type_gap_family_report", path)
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)

X = np.linspace(-0.02, 0.02, 17)
POINTS = [[float(x), 0.0, 0.0] for x in X]
CORE = np.abs(X) <= 0.010 + 1e-14
HDIV = np.stack([np.zeros_like(X), np.zeros_like(X), 1.0 + 0.5 * X ** 2], axis=1)


def _payload(name, layers, *, order=2.0, hdiv_shift=0.0, mesh_sha="m"):
    # each FEM route differs from HDiv by an amount that decays as N^-order
    fem = {}
    for engine, amplitude in (("reduced_a", 3.0e-3), ("mixed_total_reduced_omega", 5.0e-3)):
        offset = amplitude * (6.0 / layers) ** order
        fem[engine] = (HDIV + np.array([0.0, 0.0, offset])).tolist()
    fields = {"hdiv_mmm": (HDIV + hdiv_shift).tolist(), **fem}
    core = {engine: np.asarray(field)[CORE] for engine, field in fields.items()}
    pairs = {}
    for left in report.ENGINES:
        for right in report.ENGINES:
            if left < right:
                pairs[f"{left}__vs__{right}"] = {
                    "relative_rms": report.convergence._relative_rms(core[left], core[right])}
    return {
        "schema": report.LEVEL_SCHEMA, "mode": "linear", "machine": "fake",
        "radia_version": "0", "comparison_contract": {"c": 1},
        # the real contracts carry the level's mesh hashes; they must not
        # make the levels look like different experiments
        "engine_checkpoint_contracts": {
            engine: {"mode": "linear", "fem_order": 1, "linear_solver": "direct",
                     "kelvin_domain_vol_sha256": "k" + name,
                     "iron_vol_sha256": "i", "observation_points": POINTS,
                     "implementation_sha256": {"x": name}}
            for engine in report.ENGINES},
        "mesh_result_sha256": mesh_sha + name,
        "observation_points_m": POINTS, "gap_core_half_length_m": 0.010,
        "median_plane_projected_fields_T": fields,
        "pairwise_median_projected_gap_core": pairs,
        "maximum_gap_core_pairwise_relative_rms": max(p["relative_rms"] for p in pairs.values()),
        "engines": {engine: {"ndof": 10, "mesh_elements": 5, "runtime_s": 1.0}
                    for engine in report.ENGINES},
    }


def _manifest(names=("n06", "n12", "n24"), layers=(6, 12, 24)):
    return {"schema": report.FAMILY_SCHEMA, "passed": True, "gap_refinement_ratio": 2.0,
            "levels": [{"name": n, "gap_layers": l, "gap_elements": 100 * l,
                        "mesh_result_sha256": "m" + n} for n, l in zip(names, layers)]}


def _family(**kwargs):
    manifest = _manifest()
    results = {row["name"]: _payload(row["name"], row["gap_layers"], **kwargs)
               for row in manifest["levels"]}
    paths = {name: path for name in results}
    return manifest, results, paths


def test_second_order_gap_convergence_is_recovered():
    manifest, results, paths = _family(order=2.0)
    out = report.analyze(manifest, results, paths, hdiv_identity_tolerance=1e-9)
    assert out["passed"]
    assert out["hdiv_identity"]["identical"]
    for engine in report.FEM_ENGINES:
        row = out["fem_convergence_in_gap_layers"][engine]
        assert row["contracting"]
        assert row["observed_order"] == pytest.approx(2.0, abs=1e-6)
        assert out["distance_to_hdiv_shrinks_every_level"][engine]


def test_hdiv_that_moves_with_the_gap_fails_the_gate():
    manifest, results, paths = _family()
    results["n12"] = _payload("n12", 12, hdiv_shift=1e-6)
    out = report.analyze(manifest, results, paths, hdiv_identity_tolerance=1e-9)
    assert not out["passed"]
    assert not out["checks"]["hdiv_mmm_unchanged_by_gap_refinement"]


def test_result_on_another_mesh_contract_is_rejected():
    manifest, results, paths = _family()
    results["n24"]["mesh_result_sha256"] = "other"
    with pytest.raises(RuntimeError, match="mesh contract"):
        report.analyze(manifest, results, paths, hdiv_identity_tolerance=1e-9)


def test_levels_under_different_contracts_are_rejected():
    manifest, results, paths = _family()
    results["n24"]["radia_version"] = "1"
    with pytest.raises(RuntimeError, match="contract or Radia version"):
        report.analyze(manifest, results, paths, hdiv_identity_tolerance=1e-9)


def test_a_level_solved_with_another_linear_solver_is_rejected():
    manifest, results, paths = _family()
    results["n24"]["engine_checkpoint_contracts"]["reduced_a"]["linear_solver"] = "ams"
    with pytest.raises(RuntimeError, match="contract or Radia version"):
        report.analyze(manifest, results, paths, hdiv_identity_tolerance=1e-9)


def test_scale_family_manifest_is_refused():
    manifest, results, paths = _family()
    manifest["schema"] = "radia.validation.c-type-cubit-mesh-family.v1"
    with pytest.raises(RuntimeError, match="gap family"):
        report.analyze(manifest, results, paths, hdiv_identity_tolerance=1e-9)
