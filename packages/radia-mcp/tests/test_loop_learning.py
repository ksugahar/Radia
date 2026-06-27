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
        "acoustic_impedance_power",
        "rf_acoustic_passivity",
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

    force = get_loop_learning_documentation("force_moment")
    assert "Lorentz force" in force
    assert "coenergy" in force
    assert "absolute tolerance near zero crossings" in force

    motor = get_loop_learning_documentation("motor_airgap_torque")
    assert "tau(theta) = Br(theta)*Bt(theta)/mu0" in motor
    assert "T = r^2*L*integral tau(theta) dtheta" in motor
    assert "phi = pi/2" in motor
    assert "air_gap_shear_torque_from_angle_samples" in motor

    acoustic = get_loop_learning_documentation("acoustic_impedance_power")
    assert "R = (Zs - Z0)/(Zs + Z0)" in acoustic
    assert "absorption = 1 - |R|^2" in acoustic
    assert "P_boundary" in acoustic
    assert "acoustic_impedance_reflection_summary" in acoustic

    rf = get_loop_learning_documentation("rf_acoustic_passivity")
    assert "S^H S" in rf
    assert "Purely reactive impedance" in rf


def test_loop_learning_closure_prevents_overclaiming():
    doc = get_loop_learning_documentation("mcp_closure")

    assert "collected" in doc
    assert "encoded" in doc
    assert "verified" in doc
    assert "learned" in doc
    assert "If only cross-validation files were written" in doc
    assert "Apply the labels per lane" in doc
    assert "Apply the labels per slot" in doc

    overview = get_loop_learning_documentation("overview")
    assert "every slot boundary" in overview
    assert "Do not wait until a full loop is over" in overview


def test_loop_learning_mcp_tool_dispatches_without_private_provenance():
    doc = loop_learning("all")

    assert "W:\\" not in doc
    assert "S:\\" not in doc
    assert "_crossval" not in doc
    assert "learned" in doc
