import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_pyltspice_matlab_manifest_names_existing_classes():
    manifest = json.loads((ROOT / "matlab" / "pyltspice_api_compatibility.json").read_text())
    assert manifest["reference"] == {"PyLTSpice": "6.0.1", "spicelib": "1.6.3"}
    for family in manifest["families"]:
        class_name = family["matlab"].rsplit(".", 1)[-1]
        path = ROOT / "matlab" / "+radia" / "+ltspice" / f"{class_name}.m"
        assert path.is_file(), family


def test_complete_flag_cannot_hide_partial_families():
    manifest = json.loads((ROOT / "matlab" / "pyltspice_api_compatibility.json").read_text())
    if manifest["complete"]:
        assert all(family["status"] == "covered" for family in manifest["families"])
