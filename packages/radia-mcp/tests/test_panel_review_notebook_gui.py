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
    assert "masked Simulink block" in doc
    assert "radia.simulink.buildLibrary" in doc
    assert "radia.simulink.application_config.v1" in doc
    assert "radia.simulink.application_run.v1" in doc
    assert "MEX is optional" in doc
    assert "IH" in doc


def test_panel_review_exposes_cubit_panels_migration_route():
    assert "cubit_panels_migration" in TOPICS
    doc = get_panel_review_documentation("cubit_panels_migration")

    assert "validation_test/induction_heating/cubit_panels_legacy" in doc
    assert "examples/cubit_panels" in doc
    assert "src/radia/panels/samples/em/c_type_dipole" in doc
    assert "src/radia" in doc
    assert "validation_test" in doc
    assert "result-saved docs notebooks" in doc
    assert "coil_dipole.py" in doc
    assert "verify_*.py" in doc
    assert "create_induction_model.py" in doc


def test_standalone_panels_redirects_to_simulink_block_builder():
    doc = get_standalone_panels_documentation("build_notebook_gui")

    assert "DesignSpec" in doc
    assert "masked block" in doc
    assert "radia.simulink.buildLibrary" in doc
    assert "MEX/ROM" in doc
    assert "Python/headless CLI" in doc


def test_standalone_panels_redirects_to_cubit_panels_route():
    doc = get_standalone_panels_documentation("cubit_panels_migration")

    assert "examples/cubit_panels" in doc
    assert "validation_test" in doc
    assert "src/radia/panels/samples/" in doc
    assert "Simulink library" in doc


def test_standalone_panels_documents_ih_monitor_bus_boundary():
    doc = get_standalone_panels_documentation("ih_methods")

    assert "IHMonitorBusV1" in doc
    assert "RadiaMonitorHeaderV1" in doc
    assert "cell-weighted" in doc
    assert "raw heat-density `q`" in doc
