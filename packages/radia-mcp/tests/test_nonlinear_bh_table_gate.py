from __future__ import annotations

import copy
import hashlib
import json

from radia_mcp.radia_ngsolve.nonlinear_bh_table_gate import (
    nonlinear_bh_canonical_table_gate,
)
from radia_mcp.radia_ngsolve.server import (
    nonlinear_bh_canonical_table_gate as mcp_gate,
)


def _contract() -> dict:
    rows = [
        [0.0, 0.0],
        [100.0, 0.5],
        [250.0, 1.0],
        [900.0, 1.5],
        [8000.0, 2.0],
    ]
    digest = hashlib.sha256(
        json.dumps(rows, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    return {
        "schema": "radia-nonlinear-bh-table/v1",
        "material_model": "single_valued_isotropic_soft_magnetic",
        "column_order": ["H", "B"],
        "units": {"H": "A/m", "B": "T"},
        "interpolation": "monotone_pchip",
        "extrapolation": "vacuum_slope",
        "rows": rows,
        "canonical_table_sha256": digest,
    }


def test_canonical_si_table_is_accepted_through_python_and_mcp() -> None:
    contract = _contract()
    result = nonlinear_bh_canonical_table_gate(contract)
    mcp_result = json.loads(mcp_gate(contract))

    assert result["status"] == "accepted"
    assert mcp_result["status"] == "accepted"
    assert result["checks"]["canonical_digest_is_present_and_matches"] is True


def test_missing_units_are_rejected() -> None:
    contract = _contract()
    del contract["units"]

    result = nonlinear_bh_canonical_table_gate(contract)

    assert result["status"] == "rejected"
    assert result["checks"]["units_are_explicit_si"] is False


def test_swapped_column_order_is_rejected_even_when_values_look_plausible() -> None:
    contract = _contract()
    contract["column_order"] = ["B", "H"]

    result = nonlinear_bh_canonical_table_gate(contract)

    assert result["status"] == "rejected"
    assert result["checks"]["column_order_is_h_then_b"] is False


def test_nonmonotone_h_and_stale_digest_are_rejected() -> None:
    contract = copy.deepcopy(_contract())
    contract["rows"][3][0] = 200.0

    result = nonlinear_bh_canonical_table_gate(contract)

    assert result["status"] == "rejected"
    assert result["checks"]["h_is_strictly_increasing"] is False
    assert result["checks"]["canonical_digest_is_present_and_matches"] is False


def test_single_valued_gate_rejects_sub_vacuum_soft_magnetic_rows() -> None:
    contract = copy.deepcopy(_contract())
    contract["rows"][-1] = [1.0e7, 0.1]
    contract["canonical_table_sha256"] = hashlib.sha256(
        json.dumps(contract["rows"], separators=(",", ":")).encode("ascii")
    ).hexdigest()

    result = nonlinear_bh_canonical_table_gate(contract)

    assert result["status"] == "rejected"
    assert result["checks"]["soft_magnetic_secant_mu_is_at_least_vacuum"] is False
