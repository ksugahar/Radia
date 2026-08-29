import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUITE = ROOT / "validation_test" / "c_type_three_engine"


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
    assert '"hdiv_gram_eps": float(options.hdiv_gram_eps)' in runner
    assert '"nonlinear_converged": nonlinear_converged' in runner
    assert 'nonlinear_stats.get("nonlinear_converged_final_stage", False)' in runner
    assert '".reduced-a-checkpoint.json"' in runner
    assert '".omega-checkpoint.json"' in runner
    assert '"implementation_sha256"' in runner
    assert '"radia_pybind"' in runner
    assert 'dirichlet="GND"' in runner
    assert "finite_outer_air_box_forbidden" in runner
    assert "outer_boundary" not in runner


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
