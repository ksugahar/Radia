from radia_mcp.differential_forms.mathematica_recipes_knowledge import (
    get_mathematica_recipes_documentation,
)


def test_hodograph_differential_geometry_runbook_is_exposed() -> None:
    doc = get_mathematica_recipes_documentation("differential_geometry")
    assert "weakform_hodge.wls" in doc
    assert "hodograph.wls" in doc
    assert "canonical.wls" in doc
    assert "surface_derham.wls" in doc


def test_hodograph_symbolic_topics_are_exposed() -> None:
    assert "hodograph.wls" in get_mathematica_recipes_documentation("hodograph")
    assert "canonical.wls" in get_mathematica_recipes_documentation("canonical")
    assert "surface_derham.wls" in get_mathematica_recipes_documentation("hoibc")
    assert "weakform_hodge.wls" in get_mathematica_recipes_documentation(
        "weakform_hodge"
    )
