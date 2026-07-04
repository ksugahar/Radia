"""Response regression tests for mcp-server-gmsh knowledge tools."""

from radia_mcp.gmsh.gmsh_knowledge import get_gmsh_documentation


def test_gmsh_policy_documents_v41_msh_and_geo_launch_contract():
    body = get_gmsh_documentation("policy")

    assert "GmshPostExport.write()` and `vol2msh()` emit .msh v4.1" in body
    assert "Post-processing launch artifact: `case.geo`" in body
    assert "defaults to v2.2" not in body


def test_gmsh_overview_documents_python_wrapper_and_geo_association():
    body = get_gmsh_documentation("overview")

    assert "canonical launcher is `gmsh` on PATH" in body
    assert "Python wrapper" in body
    assert "Do not hard-code" in body
    assert r"`C:\gmsh.exe`" in body
    assert "`.geo` is the primary file association" in body
    assert "`.msh` association is optional raw mesh/data inspection" in body


def test_gmsh_opt_file_documents_exact_filename_autoload_contract():
    body = get_gmsh_documentation("opt_file")

    assert "Opening `case.msh` looks for `case.msh.opt`" in body
    assert "opening `case.geo` looks" in body
    assert "for `case.geo.opt`" in body
    assert "A sidecar named only `case.opt` is not an auto-load" in body
    assert "artifact writers should" in body
    assert "emit `case.geo.opt` next to `case.geo`" in body
    assert "prefer `case.geo` as the launch target" in body
    assert "`case.msh` as the raw mesh/field container" in body


def test_gmsh_geo_documents_double_click_safe_display_contract():
    body = get_gmsh_documentation("geo")

    assert "standard Radia post-processing launch artifact" in body
    assert "GmshPostExport.write()" in body
    assert "vol2msh()" in body
    assert "user-facing Open GMSH target is\n`case.geo`" in body
    assert "Do not rely on a plain `display.opt` being auto-loaded" in body
    assert "Put all critical display options directly in the `.geo`" in body
    assert "emit the exact auto-load twin `display.geo.opt`" in body
    assert "General.ClipOnlyVolume = 1" in body
    assert "View[0].Clip = 1" in body


def test_gmsh_workflow_documents_vol_to_existing_radia_msh_exporters():
    body = get_gmsh_documentation("workflow")

    assert "Do not plan a Radia workflow around GMSH directly opening Netgen `.vol`" in body
    assert "do not add a GMSH-side `.vol` reader/plugin" in body
    assert "Radia already owns the `.msh` output path" in body
    assert "GmshPostExport.write()" in body
    assert "vol2msh()" in body
    assert "cubit_mesh_export.export_Gmsh_ver4" in body
    assert "existing Radia .msh exporter -> case.geo" in body
    assert "GMSH v4.1 `.msh`" in body
    assert "Historical v2.2 snippets are legacy" in body


def test_gmsh_pitfalls_include_geo_opt_missing_settings_failure_mode():
    body = get_gmsh_documentation("pitfalls")

    assert ".geo Opens but the .opt Settings Are Missing" in body
    assert "`display.geo` can auto-load `display.geo.opt`, but not" in body
    assert "`display.opt`" in body
    assert "the `gmsh` command" in body
    assert "`View[0].Clip` clips the displayed" in body
    assert "Avoid relying on UserChoice" in body
    assert "Register\n`.geo` as the primary Radia post-processing launch association" in body
    assert "`.msh`\nassociation optional for raw mesh/data inspection only" in body
    assert "`case.geo` is the post-processing" in body
    assert "`case.msh` is the raw mesh/data" in body
    assert "`case.msh.opt`" in body
    assert "post-processing exporters should therefore emit\n`.geo` by default" in body
