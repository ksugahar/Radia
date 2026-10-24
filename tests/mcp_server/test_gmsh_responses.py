"""Response regression tests for mcp-server-gmsh knowledge tools."""

from radia_mcp.gmsh.gmsh_knowledge import get_gmsh_documentation


def test_gmsh_policy_documents_v42_post_only_and_geo_launch_contract():
    body = get_gmsh_documentation("policy")

    assert "GmshPostExport.write()` and `vol2msh()` emit .msh v4.2" in body
    assert "Post-processing launch artifact: `case.geo`" in body
    assert "visualization and post-processing only, NOT mesh generation" in body
    assert "**STEP -> Netgen**" in body
    assert "**STEP -> Cubit**" in body


def test_gmsh_overview_documents_qualified_runtime_and_solver_boundary():
    body = get_gmsh_documentation("overview")

    assert "qualified official Gmsh 5 with ASCII MSH 4.2" in body
    assert "checked Cubit/Netgen `.vol`; Gmsh is post-processing only" in body
    assert "gmsh-development-runtime.md" in body
    assert "Dependency requirements are declared only in the package's pyproject extras" in body
    assert "Native tools require `RADIA_GMSH_PROVENANCE`" in body
    assert "exact distribution/native version and wheel SHA256" in body
    assert "Pure-Python readers stay exact-4.2" in body


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


def test_gmsh_workflow_documents_inspected_results_and_checked_solver_inputs():
    body = get_gmsh_documentation("workflow")

    assert "Gmsh 5 / MSH 4.2 post-processing workflow" in body
    assert "existing result file exported by Radia or cubit-mesh-export" in body
    assert "Keep the original checked solver `.vol`, its label contract" in body
    assert "Gmsh does not create or import the solver mesh" in body
    assert body.index("gmsh_inspect_msh(") < body.index("gmsh_write_post_launch_artifact(")
    assert body.index("gmsh_validate_msh(") < body.index("gmsh_render(")
    assert "`field.geo`, `field.geo.opt`, `field.msh.opt`" in body
    assert "The existing MSH data are not modified" in body
    assert "They are not solver input" in body


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
