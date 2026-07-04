from radia_mcp.gmsh.gmsh_knowledge import get_gmsh_documentation


def test_gmsh_knowledge_records_zup_geo_opt_and_animation_export():
    cli = get_gmsh_documentation("cli")
    opt_file = get_gmsh_documentation("opt_file")
    animation = get_gmsh_documentation("animation")

    assert "PowerShell `-string` quoting" in cli
    assert "Unknown variable 'W'" in cli
    assert "Print 'W:/path/frame.png'; Exit;" in cli

    assert "Z-up x-z plane post view" in opt_file
    assert "General.RotationX = -68" in opt_file
    assert "General.RotationZ = 0" in opt_file
    assert "case.geo.opt" in opt_file

    assert "Programmatic PNG/GIF Export" in animation
    assert "gmsh.fltk.initialize()" in animation
    assert 'Path("C:/temp")' in animation
    assert "PostProcessing.AnimationCycle = 0" in animation
