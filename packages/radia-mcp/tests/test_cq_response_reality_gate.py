from radia_mcp.radia_ngsolve.cq_urn import cq_response_reality_gate


def _row():
    return {"method": "BDF2", "coupling_form": "JohnsonNedelec", "num_time_steps": 16,
            "num_laplace_solves": 16, "min_real_laplace_parameter": 3.0,
            "max_relative_residual": 4e-18, "interior_imag_relative": 4e-14,
            "exterior_imag_relative": 2e-14, "double_layer_included": True}


def test_cq_response_reality_accepts_coupled_real_reconstruction():
    gate = cq_response_reality_gate(_row())
    assert gate["status"] == "ok" and all(gate["checks"].values())


def test_cq_response_reality_rejects_complex_leak_and_missing_double_layer():
    row = _row(); row["exterior_imag_relative"] = 1e-2; row["double_layer_included"] = False
    gate = cq_response_reality_gate(row)
    assert gate["status"] == "needs_attention"
