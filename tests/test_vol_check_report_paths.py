"""Guards: check-vol reports record the mesh and CAD-sidecar paths with "/"
separators whatever the host OS, so committed reports do not depend on it.

Assumes: NGSolve is importable; the hdiv hex_pair fixture is a valid .vol.
"""
import os
import shutil
from pathlib import Path

import pytest

pytest.importorskip("ngsolve")
from radia.vol_check import check_consistency  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def test_report_paths_use_forward_slashes(tmp_path, monkeypatch):
    mesh_dir = tmp_path / "a" / "b"
    mesh_dir.mkdir(parents=True)
    shutil.copy(ROOT / "tests/fixtures/hdiv/hex_pair.vol", mesh_dir / "m.vol")
    (mesh_dir / "m.vol.json").write_text('{"volumes": {}, "areas": {}}', encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    given = os.path.join("a", "b", "m.vol")  # native separators, backslashes on Windows
    report = check_consistency(given)
    paths = {"vol_file": report["vol_file"], "json_file": report["json_file"],
             "cad_reference.file": report["cad_reference"]["file"]}
    expected = {"vol_file": "a/b/m.vol", "json_file": "a/b/m.vol.json",
                "cad_reference.file": "a/b/m.vol.json"}
    assert paths == expected, (paths, expected)
