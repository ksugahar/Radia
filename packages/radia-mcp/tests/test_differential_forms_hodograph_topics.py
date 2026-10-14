from radia_mcp.differential_forms.mathematica_recipes_knowledge import (
    get_mathematica_recipes_documentation,
)


def test_hodograph_differential_geometry_runbook_is_exposed() -> None:
    doc = get_mathematica_recipes_documentation("differential_geometry")
    assert "weakform_hodge.wls" in doc
    assert "hodograph.wls" in doc
    assert "canonical.wls" in doc
    assert "surface_derham.wls" in doc
    assert "dtn_geometry.wls" in doc
    assert "cohomology.wls" in doc
    assert "curve_surface_basics.wls" in doc
    assert "curve_surface_reading_guide.md" in doc


def test_hodograph_symbolic_topics_are_exposed() -> None:
    assert "hodograph.wls" in get_mathematica_recipes_documentation("hodograph")
    assert "canonical.wls" in get_mathematica_recipes_documentation("canonical")
    assert "surface_derham.wls" in get_mathematica_recipes_documentation("hoibc")
    assert "dtn_geometry.wls" in get_mathematica_recipes_documentation("dtn")
    assert "weakform_hodge.wls" in get_mathematica_recipes_documentation(
        "weakform_hodge"
    )


def test_curve_surface_textbook_guide_is_exposed() -> None:
    doc = get_mathematica_recipes_documentation("curve_surface")
    assert "curve_surface_basics.wls" in doc
    assert "Mathematica 曲線と曲面の微分幾何.pdf" in doc
    assert "curve_surface_reading_guide.md" in doc
    assert "K = n x grad_Gamma psi" in doc
    assert "div_Gamma K = 0" in doc
    assert "hodograph.wls" in doc
    assert "surface_derham.wls" in doc
