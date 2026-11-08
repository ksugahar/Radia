"""Keep the public lifecycle fixture's quality options on the canonical API."""
import inspect
import json
from pathlib import Path

import pytest


def test_projection_lifecycle_quality_options_match_checker():
    root = Path(__file__).resolve().parents[1]
    from cubit_mesh_export import check as checker
    contract = json.loads((root / "validation_test/hmatrix_cuda/affine_projection_lifecycle.contract.json").read_text())
    checker._load_label_contract(contract)
    signature = inspect.signature(checker.check_mesh_quality)
    signature.bind(object(), **contract["quality"])
    assert contract["quality"]["require_tetrahedra"] is True
    assert contract["metadata"]["length_unit"] == "m"
    with pytest.raises(TypeError, match="require_tet"):
        signature.bind(object(), require_tet=True)
