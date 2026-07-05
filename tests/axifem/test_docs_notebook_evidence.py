import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs" / "axifem"


def test_axifem_element_evidence_json_covers_all_shipping_paths():
    evidence = json.loads(
        (DOCS / "axifem_element_evidence.json").read_text(encoding="utf-8")
    )

    assert "runtime_radia_version" in evidence
    assert "executed_at_utc" in evidence
    assert evidence["pytest"]["returncode"] == 0
    assert evidence["pytest"]["passed"] >= 34

    labels = {row["Element path"] for row in evidence["evidence_matrix"]}
    assert {
        "P1 triangle",
        "Q1 quad",
        "P2 triangle",
        "Q2 quad",
        "P2 curved triangle",
        "Q2 curved quad",
    } <= labels


def test_axifem_element_evidence_notebook_is_result_bearing():
    nb = json.loads(
        (DOCS / "AXIFEM_ELEMENT_EVIDENCE.ipynb").read_text(encoding="utf-8")
    )

    radia_meta = nb["metadata"]["radia"]
    assert radia_meta["artifact_type"] == "documentation-notebook"
    assert radia_meta["outputs_policy"] == "embedded-results-must-keep-version-stamp"
    assert radia_meta["result_json"] == "docs/axifem/axifem_element_evidence.json"

    text = json.dumps(nb)
    assert "runtime_radia_version" in text
    assert f"{radia_meta['version_stamp']['pytest_passed']} passed" in text
    assert ".vol" in text
    assert "du_rham_identity" in text
    for label in [
        "P1 triangle",
        "Q1 quad",
        "P2 triangle",
        "Q2 quad",
        "P2 curved triangle",
        "Q2 curved quad",
    ]:
        assert label in text

    code_cells = [cell for cell in nb["cells"] if cell["cell_type"] == "code"]
    assert code_cells
    assert all(cell.get("outputs") for cell in code_cells)
