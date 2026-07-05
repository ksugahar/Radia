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
    assert "gmsh_post_display_contract" in opt_file
    assert "write_gmsh_post_launch_artifact" in opt_file
    assert "writeGmshPostLaunchArtifact" in opt_file
    assert "General.Clip0A/B/C/D" in opt_file

    assert "Programmatic PNG/GIF Export" in animation
    assert "gmsh.fltk.initialize()" in animation
    assert 'Path("C:/temp")' in animation
    assert "PostProcessing.AnimationCycle = 0" in animation
    assert "gmsh_animation_export.ipynb" in animation
    assert "writes GIF/MP4 movies" in animation


def test_gmsh_knowledge_records_cubit_mesh_export_geo_companion_contract():
    policy = get_gmsh_documentation("policy")

    assert 'cubit-mesh-export` `export gmsh "case.msh"`' in policy
    assert "case.geo.opt" in policy
    assert "case.msh.opt" in policy
    assert "A plain `case.opt` is not auto-loaded" in policy
    assert "gmsh_post_display_contract" in policy
    assert "cut-plane metadata" in policy
