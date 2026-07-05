from radia_mcp.radia_ngsolve.knowledge.loop_learning import (
    TOPICS,
    get_loop_learning_documentation,
)
from radia_mcp.radia_ngsolve.server import loop_learning


def test_loop_learning_topics_cover_current_loop_lessons():
    assert {
        "overview",
        "dual_lane",
        "mesh_geometry_vol",
        "force_moment",
        "motor_airgap_torque",
        "fem_bem_trace_orientation",
        "fem_bem_solver_report",
        "bem_demag_source_mesh",
        "acoustic_impedance_power",
        "rf_acoustic_passivity",
        "source_native_seed_queue",
        "autonomous_basic_learning",
        "em_force_target",
        "artifact_feedback",
        "mcp_closure",
    }.issubset(TOPICS)

    dual = get_loop_learning_documentation("dual_lane")
    assert "one artifact teaches twice" in dual
    assert "Public/open lane" in dual
    assert "Source-tool lane" in dual
    assert "private MCP or converter" in dual

    mesh = get_loop_learning_documentation("mesh_geometry_vol")
    assert "volumeelements > 0" in mesh
    assert "triangle surface elements" in mesh
    assert "tetrahedron volume elements" in mesh
    assert "register material volume blocks" in mesh
    assert "Netgen/OCC is enough for tet-only meshes" in mesh
    assert "hex+pyramid+tet" in mesh
    assert "semantic inventory gate" in mesh

    force = get_loop_learning_documentation("force_moment")
    assert "Lorentz force" in force
    assert "coenergy" in force
    assert "absolute tolerance near zero crossings" in force
    assert "model-input artifact id" in force
    assert "loaded-solution id" in force

    motor = get_loop_learning_documentation("motor_airgap_torque")
    assert "tau(theta) = Br(theta)*Bt(theta)/mu0" in motor
    assert "T = r^2*L*integral tau(theta) dtheta" in motor
    assert "phi = pi/2" in motor
    assert "air_gap_shear_torque_from_angle_samples" in motor
    assert "project/model" in motor
    assert "model input package" in motor

    orientation = get_loop_learning_documentation("fem_bem_trace_orientation")
    assert "normal_flux_artifact_id" in orientation
    assert "normal_flux_digest" in orientation
    assert "normal_flux_convention" in orientation
    assert "netgen_vol_first_order_fem_bem_trace_package_handoff" in orientation

    solver_report = get_loop_learning_documentation("fem_bem_solver_report")
    assert "linear_solver_report_artifact_id" in solver_report
    assert "linear_solver_report_digest" in solver_report
    assert "linear_solver_residual_norm" in solver_report
    assert "result_artifact_id" in solver_report
    assert "run_started_at" in solver_report
    assert "tool_version" in solver_report
    assert "notebook_source_artifact_id" in solver_report
    assert "notebook_source_digest" in solver_report
    assert "notebook_source_path" in solver_report
    assert "parameter_set_artifact_id" in solver_report
    assert "objective_observable_id" in solver_report
    assert "require_parameter_set_artifact=True" in solver_report
    assert "timing_breakdown_s" in solver_report
    assert "require_linear_solver_report=True" in solver_report

    bem = get_loop_learning_documentation("bem_demag_source_mesh")
    assert "surface_mesh_digest" in bem
    assert "surface_row_count" in bem
    assert "source_balance_digest" in bem
    assert "pm_demag_margin_screening_package_gate" in bem

    acoustic = get_loop_learning_documentation("acoustic_impedance_power")
    assert "R = (Zs - Z0)/(Zs + Z0)" in acoustic
    assert "absorption = 1 - |R|^2" in acoustic
    assert "P_boundary" in acoustic
    assert "acoustic_impedance_reflection_summary" in acoustic

    rf = get_loop_learning_documentation("rf_acoustic_passivity")
    assert "S^H S" in rf
    assert "frequency_grid_digest" in rf
    assert "model_input_artifact_id" in rf
    assert "model_input_digest" in rf
    assert "model_input_path" in rf
    assert "solver_result_artifact_provenance_timing_gate" in rf
    assert "timing_breakdown_s" in rf
    assert "run_date_utc" in rf
    assert "sweep_axis_digest" in rf
    assert "solver_configuration_digest" in rf
    assert "relative_tolerance" in rf
    assert "Purely reactive impedance" in rf

    source_native = get_loop_learning_documentation("source_native_seed_queue")
    assert "source-native example" in source_native
    assert "Generated scripts" in source_native
    assert "replay harnesses" in source_native
    assert "source_native_example" in source_native
    assert "learning_lanes" in source_native
    assert "candidate" in source_native
    assert "encoded and verified MCP changes" in source_native

    autonomous = get_loop_learning_documentation("autonomous_basic_learning")
    assert "process every queued slot" in autonomous
    assert "computed/reference/tolerance/pass" in autonomous
    assert "source-tool" in autonomous
    assert "solver-ready queue" in autonomous
    assert "build_autonomous_basic_learning_artifact" in autonomous

    em_force = get_loop_learning_documentation("em_force_target")
    assert "force_torque_motor" in em_force
    assert "parallel-wire Lorentz force" in em_force
    assert "magnetic air-gap pressure" in em_force
    assert "IPM dq torque" in em_force
    assert "build_em_force_target_artifact" in em_force
    assert "not as a claim" in em_force

    artifact_feedback = get_loop_learning_documentation("artifact_feedback")
    assert "cross_validation_artifact_to_mcp_feedback_gate" in artifact_feedback
    assert "solver_result_artifact_provenance_timing_gate" in artifact_feedback
    assert "learning_lanes.public" in artifact_feedback
    assert "notebook_source_artifact_id" in artifact_feedback
    assert "verification.public" in artifact_feedback
    assert "public-safe lesson" in artifact_feedback
    assert "learned" in artifact_feedback


def test_loop_learning_closure_prevents_overclaiming():
    doc = get_loop_learning_documentation("mcp_closure")

    assert "collected" in doc
    assert "encoded" in doc
    assert "verified" in doc
    assert "learned" in doc
    assert "If only cross-validation files were written" in doc
    assert "Apply the labels per lane" in doc
    assert "Apply the labels per slot" in doc
    assert "shared_solver_session_health_gate" in doc
    assert "shared-engine eval status" in doc
    assert "matlab.engine.find_matlab()" in doc
    assert "find_matlab()" in doc
    assert "needs_attention" in doc
    assert "started_new_process=false" in doc
    assert "killed_process=false" in doc
    assert "visible shared engine name is not enough" in doc
    assert "successful solver-session attach" in doc
    assert "solver-native preflight verdict" in doc
    assert "session-health" in doc
    assert "evidence, not as physics validation" in doc

    overview = get_loop_learning_documentation("overview")
    assert "every slot boundary" in overview
    assert "Do not wait until a full loop is over" in overview


def test_loop_learning_mcp_tool_dispatches_without_private_provenance():
    doc = loop_learning("all")

    assert "public-safe curated corpus" not in doc
    assert "public-safe curated corpus" not in doc
    assert ("_cross" + "val") not in doc
    assert "learned" in doc
