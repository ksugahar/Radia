from __future__ import annotations

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import nonlinear_inductance_sweep_gate
from test_cst_generalization_v41 import _summary_v41


_EMC = "emc_shielding_incident_transmitted_reflected_poynting_se_power_energy_generation_identity"
_CONNECTOR = "differential_connector_sdd_scc_modeconversion_passivity_loss_power_generation_identity"
_PROMOTED_CASE_IDS = (
    "v43_public_emc_shielding_incident_transmitted_reflected_poynting_se_power_energy_mismatch",
    "v43_public_differential_connector_sdd_scc_modeconversion_passivity_loss_power_mismatch",
)


def _summary_v43() -> dict:
    summary = _summary_v41()
    for index, run in enumerate(summary["runs"]):
        generation = f"emc-shield-843-{index}"
        values = {
            "incident_power_w": 100.0,
            "transmitted_power_w": 1.0,
            "reflected_power_w": 9.0,
            "absorbed_power_w": 90.0,
            "shielding_effectiveness_db": 20.0,
            "poynting_flux_orientation": "outward",
            "power_closure_residual_w": 0.0,
            "mesh_owner": f"mesh:{generation}",
        }
        run[_EMC] = {
            "emc_generation": generation,
            "generations": {name: generation for name in ("incident", "transmitted", "reflected", "poynting", "shielding", "power", "mesh", "result")},
            "values": values,
            "result_values": dict(values),
            "emc_result_sha256": "7" * 64,
            "accepted_emc_result_sha256": "7" * 64,
        }
        generation = f"connector-mixed-843-{index}"
        values = {
            "port_order": ["P1+", "P1-", "P2+", "P2-"],
            "mixed_mode_order": ["Sdd11", "Sdd21", "Scc11", "Scc21"],
            "sdd_ri": [[0.1, 0.0], [0.8, 0.0]],
            "scc_ri": [[0.2, 0.0], [0.7, 0.0]],
            "return_loss_db": [20.0, 1.938200260161128],
            "insertion_loss_db": [1.938200260161128, 3.098039199714863],
            "passivity_margin": 0.1,
            "dissipated_power_w": 0.2,
            "reference_impedance_ohm": 100.0,
            "mesh_owner": f"mesh:{generation}",
        }
        run[_CONNECTOR] = {
            "connector_generation": generation,
            "generations": {name: generation for name in ("port", "mixed_mode", "sdd", "scc", "loss", "passivity", "power", "impedance", "mesh", "result")},
            "values": values,
            "result_values": dict(values),
            "connector_result_sha256": "9" * 64,
            "accepted_connector_result_sha256": "9" * 64,
        }
    return summary


def test_v43_public_positive_emc_and_differential_connector() -> None:
    assert nonlinear_inductance_sweep_gate(_summary_v43())["status"] == "ok"


def test_v43_public_emc_rejects_power_closure_and_digest_mismatch() -> None:
    summary = _summary_v43()
    row = summary["runs"][0][_EMC]
    row["result_values"]["transmitted_power_w"] = 20.0
    row["result_values"]["power_closure_residual_w"] = -19.0
    row["accepted_emc_result_sha256"] = "d" * 64
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v43_public_connector_rejects_mode_order_and_passivity_mismatch() -> None:
    summary = _summary_v43()
    row = summary["runs"][0][_CONNECTOR]
    row["result_values"]["port_order"] = ["P1-", "P1+", "P2+", "P2-"]
    row["result_values"]["passivity_margin"] = -0.2
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v43_public_rejects_self_consistent_wrong_shielding_effectiveness() -> None:
    summary = _summary_v43()
    for run in summary["runs"]:
        row = run[_EMC]
        row["values"]["shielding_effectiveness_db"] = 10.0
        row["result_values"]["shielding_effectiveness_db"] = 10.0
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"
