"""Keep the public lifecycle fixture's quality options on the canonical API."""
import importlib.util
import inspect
import json
from pathlib import Path

import pytest


def test_projection_lifecycle_quality_options_match_checker():
    root = Path(__file__).resolve().parents[1]
    path = root / "packages/cubit-mesh-export/src/cubit_mesh_export/check.py"
    spec = importlib.util.spec_from_file_location("projection_label_checker", path)
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    contract = json.loads((root / "validation_test/hmatrix_cuda/affine_projection_lifecycle.contract.json").read_text())
    checker._load_label_contract(contract)
    signature = inspect.signature(checker.check_mesh_quality)
    signature.bind(object(), **contract["quality"])
    assert contract["quality"]["require_tetrahedra"] is True
    assert contract["metadata"]["length_unit"] == "m"
    with pytest.raises(TypeError, match="require_tet"):
        signature.bind(object(), require_tet=True)
