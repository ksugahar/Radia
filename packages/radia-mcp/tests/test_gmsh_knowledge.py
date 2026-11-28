from radia_mcp.gmsh.gmsh_knowledge import get_gmsh_documentation


def test_magnetic_movie_recipe_is_discoverable_without_optional_runtime():
    recipe = get_gmsh_documentation(" magnetic_movies ")
    assert recipe in get_gmsh_documentation("all")
    assert "magnetic_movies" in get_gmsh_documentation("unknown")
    for contract in (
        "signed instantaneous vector", "not particle trajectories",
        "Equal crossing flux does not give uniform screen density",
        "full azimuth", "reference mesh lookup", "0.2 ms",
        "decoded movie frame", "not installed-wheel tools",
        "gmsh_streamlines", "verify_movie.py",
    ):
        assert contract in recipe


def test_gmsh_knowledge_records_zup_geo_opt_and_animation_export():
    cli = get_gmsh_documentation("cli")
    opt_file = get_gmsh_documentation("opt_file")
    animation = get_gmsh_documentation("animation")

    assert "PowerShell `-string` quoting" in cli
    assert "Unknown variable 'W'" in cli
    assert "Print 'frame.png'; Exit;" in cli

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
    assert "gmsh_export_animation" in animation
    assert "PostProcessing.Link = 1" in animation
    assert "writes GIF/MP4 movies" in animation


def test_gmsh_knowledge_records_geo_companion_contract():
    policy = get_gmsh_documentation("policy")

    assert "case.geo.opt" in policy
    assert "case.msh.opt" in policy
    assert "A plain `case.opt` is not auto-loaded" in policy
    assert "gmsh_post_display_contract" in policy
    assert "cut-plane metadata" in policy


def test_gmsh_update_qualification_routes_to_canonical_sources():
    """Guards: the Gmsh 5 update workflow stays discoverable and pin-free."""
    text = get_gmsh_documentation(" update_qualification ")
    assert text in get_gmsh_documentation("all")
    assert "update_qualification" in get_gmsh_documentation("unknown")
    for contract in (
        "docs/operations/gmsh-development-runtime.md",
        ".agents/skills/gmsh-verify/SKILL.md",
        "radia_mcp.gmsh.runtime", "RADIA_GMSH_PROVENANCE",
        "gmsh_validate_msh", "gmsh_validate_geo",
        "not a native readback", "Headless", "GUI",
        "a Windows artifact does not qualify Linux",
        "runtime_provenance", "wheel_sha256", "msh42_roundtrip_probe",
        "it is not runtime provenance", "with no worker it reports none",
    ):
        assert contract in text, contract
    # The canonical guide owns URLs and hashes; this topic must not copy them.
    assert "https://" not in text
    import re
    assert not re.search(r"[0-9a-f]{64}", text)
