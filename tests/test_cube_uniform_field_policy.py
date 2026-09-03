"""Fast placement guards for the current cube uniform-field validation."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "validation_test" / "cube_uniform_field"
HDIV = CORPUS / "hdiv"


def test_cube_uniform_field_keeps_only_the_current_hdiv_lane():
    offenders = [
        path.relative_to(CORPUS)
        for path in CORPUS.rglob("*")
        if path.is_file()
        and path.name != "README.md"
        and HDIV not in path.parents
    ]
    assert not offenders

    sources = sorted(HDIV.glob("*.py"))
    assert {path.name for path in sources} == {
        "bench_hdiv_cube.py",
        "bench_hdiv_tet_cube.py",
        "bench_hex_nonlattice.py",
    }
    for path in sources:
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def test_hdiv_cube_drivers_default_to_validation_json():
    for name in ("bench_hdiv_cube.py", "bench_hdiv_tet_cube.py"):
        source = (HDIV / name).read_text(encoding="utf-8")
        assert 'default=str(_HERE.parent / "results_hdiv_' in source
        assert "rad.Solve(" not in source


def test_hdiv_cube_evidence_manifest_pins_parseable_json():
    manifest = json.loads((HDIV / "evidence_manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema"] == "radia.validation.hdiv-cube-evidence.v1"
    assert manifest["policy"]["validation_json_required"] is True

    for artifact in manifest["files"]:
        path = HDIV / artifact["path"]
        payload = path.read_bytes()
        json.loads(payload.decode("utf-8-sig"))
        assert hashlib.sha256(payload).hexdigest() == artifact["sha256"], path
