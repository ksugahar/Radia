"""Placement checks for validation families converted to script plus JSON."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_JSON_FAMILIES = (
    ROOT / "validation_test" / "eddy_current_analytical_validation",
    ROOT / "validation_test" / "force_validation",
    ROOT / "validation_test" / "ngsolve_user_meeting",
    ROOT / "validation_test" / "peec_bema_convergence",
)


def test_converted_validation_families_do_not_own_notebooks():
    offenders = [path for family in SCRIPT_JSON_FAMILIES for path in family.glob("*.ipynb")]
    assert not offenders


def test_converted_validation_families_keep_parseable_json_evidence():
    for family in SCRIPT_JSON_FAMILIES:
        records = sorted(family.glob("*.json"))
        assert records, family
        for path in records:
            json.loads(path.read_text(encoding="utf-8-sig"))
