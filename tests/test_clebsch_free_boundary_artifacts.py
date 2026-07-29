import ast
import json
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation_test" / "clebsch_legendre"
POLEFACE = VALIDATION / "results_poleface_design.json"
HORN = VALIDATION / "results_concentrator_horn.json"
SCRIPTS = (
    VALIDATION / "verify_poleface_design.py",
    VALIDATION / "verify_concentrator_horn.py",
)
README = VALIDATION / "README.md"
EM_KNOWLEDGE = (
    ROOT / "packages" / "radia-mcp" / "src" / "radia_mcp"
    / "electromagnet" / "em_knowledge.py"
)


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_free_boundary_artifacts_record_compute_runtime():
    expected_schemas = {
        POLEFACE: "radia.validation.clebsch-poleface-design.v1",
        HORN: "radia.validation.clebsch-concentrator-horn.v1",
    }
    for path, schema in expected_schemas.items():
        result = _load(path)
        assert result["schema"] == schema
        meta = result["meta"]
        assert meta["hostname"].lower() in {"mdx", "hibino"}
        assert meta["python_version"]
        assert meta["ngsolve_version"]
        assert meta["numpy_version"]
        assert datetime.fromisoformat(meta["generated_at_utc"]).tzinfo is not None


def test_poleface_artifact_keeps_design_and_fem_gates():
    result = _load(POLEFACE)
    checks = result["design_checks"]["B0_2"]
    verdict = result["verdict"]

    assert checks["J_dominant_negative"] is True
    assert checks["width_rel_err"] < 5.0e-4
    assert checks["face_Bn_rel_err_med"] < 5.0e-3
    assert 0.28 <= checks["sag_at_edge_mm"] <= 0.40
    assert verdict["improvement_factor"] > 2.5
    assert verdict["saturation_share_max"] < 1.0e-3


def test_concentrator_horn_artifact_keeps_cap_and_gain_gates():
    result = _load(HORN)
    design = result["design"]
    verdict = result["verdict"]

    assert design["J_single_sign"] is True
    assert 8.0 <= design["geometric_gain"] <= 9.6
    assert 0.02 <= verdict["gain_advantage"] <= 0.09
    assert verdict["iron_peak_at_rated_T"]["horn"] <= 1.03
    assert result["ladder"]["rel_shift"] < 0.01
    assert verdict["gain_at_3x"]["horn"] >= verdict["gain_at_3x"]["straight"]


def test_free_boundary_figure_sources_have_no_in_figure_titles():
    for path in SCRIPTS:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        title_calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr in {"set_title", "suptitle", "title"}
        ]
        assert not title_calls, f"in-figure title call remains in {path.name}"


def test_clebsch_research_anchors_use_validation_tree():
    readme = README.read_text(encoding="utf-8")
    knowledge = EM_KNOWLEDGE.read_text(encoding="utf-8")
    anchors = {
        "validation_test/feec/test_clebsch_legendre_3d.py": readme,
        "validation_test/feec/test_clebsch_hodograph_research.py": knowledge,
    }
    for relative, text in anchors.items():
        assert relative in text
        assert (ROOT / relative).is_file(), f"missing repository anchor: {relative}"
