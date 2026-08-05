import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_application_manifest_has_no_notebook_workbenches():
    manifest_path = (
        REPO_ROOT / "src" / "radia" / "panels"
        / "application_interface_manifest.json"
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    applications = {item["id"]: item for item in manifest["applications"]}
    spatial = manifest["library"]["spatial_artifact"]
    mesh_preflight = manifest["library"]["mesh_preflight"]
    assert spatial["format"] == "GMSH .msh v4.1"
    assert spatial["location"] == "application run directory"
    assert spatial["scalar_only"] == "not-applicable"
    assert mesh_preflight["checker"] == "check-vol"
    assert mesh_preflight["label_contract_schema"] == "radia.vol-label-contract.v1"
    assert mesh_preflight["report_schema"] == "cubit-mesh-export.vol-check.v1"
    assert "never inferred" in mesh_preflight["material_constants"]

    expected_blocks = {
        "radia-em",
        "radia-pcb",
        "radia-motor",
        "radia-streamfunction",
        "radia-ih",
        "radia-maglev",
    }
    for application_id in expected_blocks:
        application = applications[application_id]
        assert application["state"] == "active-simulink-block"
        assert application["block"].startswith("radia_simulink_library/Applications/")
        assert application["notebook"] is None
        assert "adapter" not in application

    ih = applications["radia-ih"]
    assert ih["backend"] == "matlab-level2+radia-mex-handles"
    assert "Level-2 MATLAB" in ih["backend_policy"]
    assert "radia_mex object handles" in ih["backend_policy"]
    assert "not a per-step fallback" in ih["backend_policy"]

    maglev = applications["radia-maglev"]
    assert maglev["backend"] == "matlab-level2-common-basis-cln"
    assert maglev["sample"] == "matlab/radia_maglev.slx"
    assert "never called per step" in maglev["backend_policy"]

    export_menu = applications["radia-export-menu"]
    assert export_menu["state"] == "active-cubit-toolbar"
    assert export_menu["notebook"] is None


def test_packaged_panel_notebook_directory_is_absent_or_empty():
    notebook_dir = REPO_ROOT / "src" / "radia" / "panels" / "notebooks"
    assert not notebook_dir.exists() or not list(notebook_dir.glob("*.ipynb"))
