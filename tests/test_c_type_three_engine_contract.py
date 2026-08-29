import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SUITE = ROOT / "validation_test" / "c_type_three_engine"


def _load_convergence_module():
    path = SUITE / "run_mesh_convergence.py"
    spec = importlib.util.spec_from_file_location("c_type_mesh_convergence", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_c_type_cad_is_cubit_authority():
    journal = (SUITE / "cad" / "c_type_iron.jou").read_text(encoding="utf-8")
    assert "create volume loft" in journal
    assert "volume all copy reflect z" in journal
    assert "netgen.occ" not in journal


def test_three_engine_runner_has_shared_physics_contract():
    runner = (SUITE / "run_three_engine.py").read_text(encoding="utf-8")
    assert 'rad.RadiaField(coil, "h")' in runner
    assert 'rad.RadiaField(coil, "b")' in runner
    assert "HDiv-MMM" in runner
    assert "HCurl reduced-A" in runner
    assert "H1 Omega-reduced-Omega" in runner
    assert "fixed_mesh_equality_claimed" in runner
    assert "pairwise_raw_full_tube" in runner
    assert "pairwise_median_projected_gap_core" in runner
    assert "reflection_diagnostics" in runner
    assert "median_plane_one_sided_tangential_relative_rms" in runner
    assert "not a reflection-symmetry error" in runner
    assert 'kelvin_vol = mesh_dir / "kelvin_domain.vol"' in runner
    assert "has_kelvin_identification" in runner
    assert '"--hdiv-gram-eps"' in runner
    assert 'default=1.0e-14' in runner
    assert '"--primary-only"' in runner
    assert '"--reduced-a-relax"' in runner
    assert '"relax": float(options.reduced_a_relax)' in runner
    assert 'primary_pair_name = "hdiv_mmm__vs__omega_reduced_omega"' in runner
    assert '"primary_gap_core_relative_rms"' in runner
    assert '"bh_interpolation_shared"' in runner
    assert '"hdiv_gram_eps": float(options.hdiv_gram_eps)' in runner
    assert '"nonlinear_converged": nonlinear_converged' in runner
    assert 'nonlinear_stats.get("nonlinear_converged_final_stage", False)' in runner
    assert '".reduced-a-checkpoint.json"' in runner
    assert '".omega-checkpoint.json"' in runner
    assert '"implementation_sha256"' in runner
    assert '"radia_pybind"' in runner
    assert '"coil_builder"' in runner
    assert '"primary_accuracy_passed"' in runner
    assert '"all_pairwise_within_tolerance"' in runner
    assert 'dirichlet="GND"' in runner
    assert "finite_outer_air_box_forbidden" in runner
    assert "outer_boundary" not in runner


def test_omega_and_hdiv_share_the_nonlinear_bh_interpolation_contract():
    source = (
        ROOT / "src" / "radia" / "scalar_potential_solver.py"
    ).read_text(encoding="utf-8")

    assert "def _build_bh_interpolator(" in source
    assert "PchipInterpolator" in source
    assert "B_max + MU_0 * (value - H_max)" in source
    nonlinear = source.split("def solve_nonlinear(", 1)[1].split(
        "def solve_nonlinear_newton(", 1
    )[0]
    assert "B_of_H = _build_bh_interpolator(bh)" in nonlinear
    assert "interp1d" not in nonlinear

    from radia.scalar_potential_solver import MU_0, _build_bh_interpolator
    from radia.vim._nonlinear import _bh_table_funcs

    table = np.loadtxt(SUITE.parents[1] / "src" / "radia" / "panels" /
                       "samples" / "em_sample_bh.txt")
    omega_b_of_h = _build_bh_interpolator(table)
    _, hdiv_b_of_h, _, h_max, _ = _bh_table_funcs(table[:, 0], table[:, 1])
    h_values = np.geomspace(max(table[1, 0] * 1e-3, 1e-8), 4.0 * h_max, 2000)
    hdiv_values = np.where(
        h_values <= h_max,
        hdiv_b_of_h(h_values),
        table[-1, 1] + MU_0 * (h_values - h_max),
    )
    omega_values = np.asarray([omega_b_of_h(value) for value in h_values])
    np.testing.assert_array_equal(omega_values, hdiv_values)


def test_reduced_a_nonlinear_direct_path_is_symmetric_and_residual_checked():
    source = (
        ROOT / "src" / "radia" / "vector_potential_solver.py"
    ).read_text(encoding="utf-8")
    nonlinear = source.split("def solve_nonlinear(", 1)[1].split(
        "def solve_hysteresis(", 1
    )[0]
    assert "BilinearForm(fes, symmetric=True)" in nonlinear
    assert "eps = 1e-6" in nonlinear
    assert "inverse='pardisospd'" in nonlinear
    assert "maximum_linear_relative_residual" in nonlinear
    assert "reduced-A linear solve failed its relative-residual" in nonlinear
    assert "reduced-A linear solve produced a non-finite residual" in nonlinear
    assert "_build_nu_of_b_interpolator(bh_data)" in nonlinear
    assert "interp1d" not in nonlinear


def test_reduced_a_inverts_the_shared_pchip_b_of_h_law():
    from radia.scalar_potential_solver import MU_0, _build_bh_interpolator
    from radia.vector_potential_solver import _build_nu_of_b_interpolator

    table = np.loadtxt(
        SUITE.parents[1]
        / "src"
        / "radia"
        / "panels"
        / "samples"
        / "em_sample_bh.txt"
    )[:, :2]
    b_of_h = _build_bh_interpolator(table)
    nu_of_b, b_max = _build_nu_of_b_interpolator(table)
    h_samples = np.geomspace(max(table[1, 0] * 1e-3, 1e-8), 4 * table[-1, 0], 200)
    for h_value in h_samples:
        b_value = b_of_h(h_value)
        recovered_h = nu_of_b(b_value) * b_value
        np.testing.assert_allclose(recovered_h, h_value, rtol=2e-12, atol=1e-8)
    assert b_max == table[-1, 1]
    b_above = 1.5 * b_max
    expected_h = table[-1, 0] + (b_above - b_max) / MU_0
    np.testing.assert_allclose(nu_of_b(b_above) * b_above, expected_h)


def test_mesh_builder_keeps_hdiv_air_mesh_free():
    builder = (SUITE / "build_cubit_meshes.py").read_text(encoding="utf-8")
    helper = (SUITE / "cubit_reflection_mesh.py").read_text(encoding="utf-8")
    assert '"iron.vol"' in builder
    assert '"kelvin_domain.vol"' in builder
    assert "HDiv-MMM must not acquire an artificial air domain" in builder
    assert "create sphere radius" in helper
    assert "add_kelvin_cubit" in builder
    assert "outer_boundary" not in builder
    assert 'export netgen "{output.as_posix()}" order 1 overwrite' in builder
    assert "_reflection_inventory" in builder
    assert "missing_reflected_volume_elements" in builder
    assert "_kelvin_identification_inventory" in builder
    assert "maximum_pair_translation_error_m" in builder
    assert "_kelvin_fes_inventory" in builder
    assert '"slaved_free_dofs"' in builder
    assert '"trace_norm_ratio"' in builder


def test_reflection_builder_copies_meshed_volumes():
    helper = (SUITE / "cubit_reflection_mesh.py").read_text(encoding="utf-8")
    assert "copy reflect" in helper
    mesh_command = 'cubit.cmd(f"mesh volume {iron_up} {air_up}")'
    reflect_call = "iron_down = _reflect_meshed_volume(cubit, iron_up)"
    assert mesh_command in helper
    assert reflect_call in helper
    assert helper.rindex(mesh_command) < helper.rindex(reflect_call)
    assert "meshing two pre-reflected geometric halves independently is forbidden" in helper.lower()


def test_nonlinear_omega_uses_the_periodic_kelvin_h1_factory():
    source = (
        ROOT / "src" / "radia" / "scalar_potential_solver.py"
    ).read_text(encoding="utf-8")
    picard = source.split("def _picard_loop(", 1)[1].split(
        "def _eval_bh_curve(", 1
    )[0]

    assert "fes = self._make_h1_space(dirichlet)" in picard
    assert "fes = H1(" not in picard


def test_tracked_kelvin_results_preserve_pass_and_failed_gate_evidence():
    results = SUITE / "results"
    mesh = json.loads(
        (results / "lab_20260829_mesh.json").read_text(encoding="utf-8")
    )
    linear = json.loads(
        (results / "lab_20260829_linear_order2.json").read_text(encoding="utf-8")
    )
    nonlinear = json.loads(
        (results / "lab_20260829_nonlinear_order1.json").read_text(
            encoding="utf-8"
        )
    )
    nonlinear_order2 = json.loads(
        (results / "lab_20260829_nonlinear_order2_primary.json").read_text(
            encoding="utf-8"
        )
    )

    assert mesh["passed"] is True
    assert mesh["kelvin_identification"]["pair_count"] == 462
    assert mesh["kelvin_identification"]["maximum_pair_translation_error_m"] < 1e-12
    assert mesh["kelvin_fes"]["slaved_free_dofs"] == 462
    assert abs(mesh["kelvin_fes"]["trace_norm_ratio"] - 1.0) < 1e-12
    assert linear["passed"] is True
    assert linear["maximum_gap_core_pairwise_relative_rms"] < 0.03
    assert nonlinear["nonlinear_converged"] is True
    assert nonlinear["passed"] is False
    assert nonlinear["maximum_gap_core_pairwise_relative_rms"] > 0.03
    assert nonlinear_order2["passed"] is True
    assert nonlinear_order2["primary_accuracy_passed"] is True
    assert nonlinear_order2["all_pairwise_within_tolerance"] is True
    assert nonlinear_order2["nonlinear_converged"] is True
    assert nonlinear_order2["comparison_contract"]["primary_pair"] == [
        "hdiv_mmm",
        "omega_reduced_omega",
    ]
    assert nonlinear_order2["primary_gap_core_relative_rms"] < 0.003
    assert nonlinear_order2["comparison_contract"]["bh_interpolation_shared"]


def test_tracked_four_level_accuracy_certificate_is_complete():
    results = SUITE / "results"
    manifest = json.loads(
        (results / "cubit_20260830_mesh_family.json").read_text(encoding="utf-8")
    )
    certificate = json.loads(
        (
            results
            / "mdx_hibino_20260830_nonlinear_order2_convergence.json"
        ).read_text(encoding="utf-8")
    )

    assert manifest["passed"] is True
    assert [row["name"] for row in manifest["levels"]] == [
        "coarse",
        "medium",
        "fine",
        "finer",
    ]
    assert certificate["schema"] == (
        "radia.validation.c-type-absolute-accuracy-certificate.v1"
    )
    assert certificate["passed"] is True
    assert certificate["claim"]["analytic_absolute_truth_claimed"] is False
    assert certificate["convergence_levels"] == ["medium", "fine", "finer"]
    assert all(certificate["checks"].values())
    assert certificate["fine_mesh_maximum_pairwise_relative_rms"] < 0.005
    assert certificate["combined_relative_numerical_uncertainty"] < 0.01
    assert certificate["reproducibility"]["reference_machine"].casefold() == "mdx"
    assert certificate["reproducibility"]["replicate_machine"].casefold() == (
        "hibino"
    )
    assert certificate["reproducibility"]["maximum_relative_rms"] < 1e-9
    assert all(
        row["convergence_levels"] == ["medium", "fine", "finer"]
        for row in certificate["engine_convergence"].values()
    )


def test_mesh_family_builder_uses_one_geometric_cubit_sequence():
    source = (SUITE / "build_mesh_family.py").read_text(encoding="utf-8")
    assert 'DEFAULT_LEVELS = (' in source
    assert '("coarse", 1.25)' in source
    assert '("medium", 1.00)' in source
    assert '("fine", 0.80)' in source
    assert '("finer", 0.64)' in source
    assert 'from build_cubit_meshes import CUBIT_DEFAULT, build' in source
    assert '"geometric_size_sequence": True' in source
    assert '"element_counts_strictly_increase"' in source
    assert "netgen.occ" not in source


def test_accuracy_certificate_requires_all_three_converged_routes():
    source = (SUITE / "run_mesh_convergence.py").read_text(encoding="utf-8")
    assert 'ENGINES = ("hdiv_mmm", "reduced_a", "omega_reduced_omega")' in source
    assert 'set(payload.get("engines", {})) != set(ENGINES)' in source
    assert 'if not payload.get("nonlinear_converged", False)' in source
    assert '"mesh_convergence_passed"' in source
    assert '"fine_mesh_cross_formulation_passed"' in source
    assert '"combined_uncertainty_passed"' in source
    assert '"independent_host_reproducibility_passed"' in source
    assert '"independent_host_replicate_required": True' in source
    assert 'reproducibility_passed = False' in source
    assert 'accuracy replication must use an independent host' in source
    assert 'replicate result uses a different Radia version' in source
    assert 'replicate result uses different implementation inputs' in source
    assert '"analytic_absolute_truth_claimed": False' in source
    assert '"convergence_levels": level_names[-3:]' in source
    assert 'subprocess.Popen(' in source
    assert '"--primary-only"' not in source


def test_accuracy_certificate_uses_the_last_three_of_four_mesh_levels():
    module = _load_convergence_module()
    level_names = ["coarse", "medium", "fine", "finer"]

    def payload(value):
        return {
            "observation_points_m": [[0.0, 0.0, 0.0]],
            "gap_core_half_length_m": 0.01,
            "median_plane_projected_fields_T": {
                "hdiv_mmm": [[0.0, 0.0, value]],
            },
        }

    payloads = [payload(value) for value in (-20.0, 1.0, 1.05, 1.075)]
    result = module._engine_convergence(
        payloads, level_names, "hdiv_mmm", refinement_ratio=1.25
    )

    assert result["convergence_levels"] == ["medium", "fine", "finer"]
    assert result["penultimate_increment_absolute_rms_T"] == pytest.approx(0.05)
    assert result["last_increment_absolute_rms_T"] == pytest.approx(0.025)
    assert result["contracting"] is True
    assert result["observed_order"] > 0.5


def test_accuracy_replication_rejects_same_host_and_software_drift():
    module = _load_convergence_module()
    base = {
        "machine": "mdx",
        "mesh_result_sha256": "mesh",
        "comparison_contract": {"order": 2},
        "radia_version": "4.95.71",
        "engine_checkpoint_contracts": {"hdiv_mmm": {"implementation": "a"}},
        "observation_points_m": [[0.0, 0.0, 0.0]],
        "gap_core_half_length_m": 0.01,
        "median_plane_projected_fields_T": {
            engine: [[0.0, 0.0, 1.0]] for engine in module.ENGINES
        },
    }
    with pytest.raises(RuntimeError, match="independent host"):
        module._replicate_metrics(base, dict(base))

    different_host = {**base, "machine": "hibino"}
    drifted = {**different_host, "radia_version": "4.95.70"}
    with pytest.raises(RuntimeError, match="different Radia version"):
        module._replicate_metrics(base, drifted)

    metrics = module._replicate_metrics(base, different_host)
    assert metrics["maximum_relative_rms"] == 0.0
