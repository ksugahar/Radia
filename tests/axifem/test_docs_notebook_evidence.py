import ast
import json
import os
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs" / "axifem"
VALIDATION = ROOT / "validation_test" / "axifem"


def test_axifem_element_evidence_json_covers_all_shipping_paths():
    evidence = json.loads(
        (VALIDATION / "axifem_element_evidence.json").read_text(encoding="utf-8")
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
    assert (
        radia_meta["result_json"]
        == "validation_test/axifem/axifem_element_evidence.json"
    )

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


@pytest.mark.parametrize('source,names,expected', [
    ('axifem/research/validate_q2_codegen.py', ['PROTO'],
     {'PROTO': 'maglev/research_cln/axifem'}),
    ('axifem/research/verification/test_3way_cauer_cross_validation.py', ['BEM_DIR', 'BEM_REF'],
     {'BEM_REF': 'maglev/research_cln/ngsolve_validation/bem_disk_axisym_cauer_python_results.json'}),
    ('axifem/research/verification/test_hiruma_disk.py', ['BEM_TAU_REF_PATH'],
     {'BEM_TAU_REF_PATH': 'maglev/research_cln/ngsolve_validation/bem_disk_axisym_v3_refined.json'}),
    ('maglev/research_cln/ngsolve_validation/disk_bem_cauer.py', ['_VERIF', 'Q1_JSON', 'Q2_JSON'],
     {'Q1_JSON': 'axifem/research/verification/test_hiruma_disk_q1_results.json',
      'Q2_JSON': 'axifem/research/verification/test_hiruma_disk_q2_results.json'}),
])
def test_research_reference_paths_follow_a_relocated_checkout(tmp_path, monkeypatch, source, names, expected):
    validation = ROOT / 'validation_test'
    filename = validation / source
    tree = ast.parse(filename.read_text(encoding='utf-8'))
    # Evaluate only path declarations, never import or execute the numerical scripts.
    declarations = [node for node in tree.body if isinstance(node, ast.Assign)
                    and any(isinstance(t, ast.Name) and t.id in names for t in node.targets)]
    assert len(declarations) == len(names)
    relocated = tmp_path / 'relocated checkout' / 'validation_test'
    monkeypatch.chdir(tmp_path)
    scope = {'Path': Path, 'os': os, '__file__': str(relocated / source)}
    exec(compile(ast.Module(body=declarations, type_ignores=[]), str(filename), 'exec'), scope)
    for name, relative in expected.items():
        assert Path(scope[name]) == relocated / relative
        # JSON results can be regenerated; their owning directories must ship.
        target = validation / relative
        assert (target.parent if target.suffix else target).is_dir()
    assert (validation / 'maglev/research_cln/axifem/axifem_quad_q2.py').is_file()
    assert (validation / 'maglev/research_cln/ngsolve_validation/bem_disk_axisym_cauer.wls').is_file()
