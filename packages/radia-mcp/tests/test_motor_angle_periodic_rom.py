from radia_mcp.motor.angle_periodic_rom_knowledge import (
    SECTIONS,
    get_angle_periodic_rom_knowledge,
)


def test_angle_periodic_motor_rom_knowledge_is_complete_and_fail_loud():
    assert set(SECTIONS) == {
        "architecture",
        "face_policy",
        "angle_rom",
        "time_domain",
        "ports",
        "mesh_gate",
        "validation",
        "limits",
    }
    assert "cycle basis" in get_angle_periodic_rom_knowledge("face_policy")
    assert "positive-real CLN" in get_angle_periodic_rom_knowledge("time_domain")
    assert "C ABI version 1" in get_angle_periodic_rom_knowledge("ports")
    assert "Unknown topic" in get_angle_periodic_rom_knowledge("missing")
    assert len(get_angle_periodic_rom_knowledge("all")) > 2000
