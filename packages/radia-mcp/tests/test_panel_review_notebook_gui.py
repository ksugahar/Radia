from radia_mcp.panel_review.panel_review_knowledge import (
    TOPICS,
    get_panel_review_documentation,
)
from radia_mcp.radia_ngsolve.knowledge.standalone_panels import (
    get_standalone_panels_documentation,
)


def test_panel_review_exposes_notebook_gui_builder_topic():
    assert "build_notebook_gui" in TOPICS
    doc = get_panel_review_documentation("build_notebook_gui")

    assert "DesignSpec" in doc
    assert "CommandWorkbench" in doc
    assert "panels/notebooks/radia_<app>.ipynb" in doc
    assert "panels/samples/<app>/..." in doc
    assert "NotebookFieldSpec" in doc
    assert "validation_test/panels/test_notebook_workbench.py" in doc
    assert "radia.notebook_panel_run.v2" in doc
    assert "IHDesignSpec" in doc
    assert "IHWorkbench" in doc
    assert "pointer events" in doc
    assert "PySide6/PyQt" in doc


def test_panel_review_exposes_cubit_panels_migration_route():
    assert "cubit_panels_migration" in TOPICS
    doc = get_panel_review_documentation("cubit_panels_migration")

    assert "validation_test/induction_heating/cubit_panels_legacy" in doc
    assert "examples/cubit_panels" in doc
    assert "panels/samples/em/c_type_dipole" in doc
    assert "src/radia" in doc
    assert "validation_test" in doc
    assert "result-saved docs notebooks" in doc
    assert "coil_dipole.py" in doc
    assert "verify_*.py" in doc
    assert "create_induction_model.py" in doc


def test_standalone_panels_redirects_to_notebook_gui_builder():
    doc = get_standalone_panels_documentation("build_notebook_gui")

    assert "panel_review(topic=\"build_notebook_gui\")" in doc
    assert "DesignSpec" in doc
    assert "CommandWorkbench" in doc
    assert "panels/calc_*.py" in doc
    assert "RADIA-IH.ipynb" in doc


def test_standalone_panels_redirects_to_cubit_panels_route():
    doc = get_standalone_panels_documentation("cubit_panels_migration")

    assert "examples/cubit_panels" in doc
    assert "panel_review(topic=\"cubit_panels_migration\")" in doc
    assert "validation_test/induction_heating" in doc
    assert "panels/samples/em/c_type_dipole" in doc
