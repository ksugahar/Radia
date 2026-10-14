from radia_mcp.panel_review.panel_review_knowledge import (
    TOPICS,
    get_panel_review_documentation,
)


def test_panel_review_exposes_current_application_block_builder():
    assert set(TOPICS) == {
        "overview",
        "build_application_block",
        "cubit_boundary",
        "workflow",
    }
    doc = get_panel_review_documentation("build_application_block")

    assert "DesignSpec" in doc
    assert "masked Simulink block" in doc
    assert "radia.simulink.buildLibrary" in doc
    assert "radia.simulink.application_config.v1" in doc
    assert "radia.simulink.application_run.v1" in doc
    assert "MEX is optional" in doc
    assert "IH" in doc


def test_panel_review_keeps_cubit_and_qt_at_the_external_boundary():
    doc = get_panel_review_documentation("cubit_boundary")

    assert "SAT" in doc
    assert "STEP" in doc
    assert "check-vol" in doc
    assert "masked Simulink block" in doc
    assert "embedded PySide6" in doc
    assert "Radia itself must not acquire a Qt/PySide dependency" in doc
