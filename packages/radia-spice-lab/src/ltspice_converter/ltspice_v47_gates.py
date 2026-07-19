"""Replay identities for hierarchical currents and mixed-analysis rows."""

from __future__ import annotations

from collections.abc import Mapping


def _digest(value: object) -> bool:
    text = str(value or "")
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text.lower())


def _generation_closed(contract: Mapping[str, object], *names: str) -> bool:
    generation = str(contract.get("generation_id") or "")
    return bool(generation) and all(contract.get(name) == generation for name in names)


def _unique_text_list(value: object) -> bool:
    return isinstance(value, list) and bool(value) and all(isinstance(item, str) and item for item in value) and len(value) == len(set(value))


def _branch_ok(contract: Mapping[str, object]) -> bool:
    paths = contract.get("hierarchy_paths")
    currents = contract.get("branch_current_order")
    return (
        _generation_closed(
            contract,
            "hierarchy_generation_id",
            "current_generation_id",
            "kcl_generation_id",
            "result_generation_id",
        )
        and _unique_text_list(paths)
        and _unique_text_list(currents)
        and len(paths) == len(currents)
        and contract.get("result_hierarchy_paths") == paths
        and contract.get("result_branch_current_order") == currents
        and contract.get("current_sign_convention") == "positive_into_pin1"
        and contract.get("result_current_sign_convention") == contract.get("current_sign_convention")
        and str(contract.get("kcl_owner") or "").startswith("node:")
        and contract.get("result_kcl_owner") == contract.get("kcl_owner")
        and _digest(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _step_rows_ok(contract: Mapping[str, object]) -> bool:
    ac_steps = contract.get("ac_step_tuples")
    tran_steps = contract.get("tran_step_tuples")
    trace_rows = contract.get("trace_row_keys")
    measure_rows = contract.get("measure_row_keys")
    return (
        _generation_closed(
            contract,
            "ac_generation_id",
            "tran_generation_id",
            "trace_generation_id",
            "measure_generation_id",
            "result_generation_id",
        )
        and isinstance(ac_steps, list)
        and bool(ac_steps)
        and contract.get("result_ac_step_tuples") == ac_steps
        and isinstance(tran_steps, list)
        and bool(tran_steps)
        and contract.get("result_tran_step_tuples") == tran_steps
        and _unique_text_list(trace_rows)
        and trace_rows == contract.get("result_trace_row_keys")
        and _unique_text_list(measure_rows)
        and measure_rows == contract.get("result_measure_row_keys")
        and trace_rows == measure_rows
        and len(trace_rows) == len(ac_steps) + len(tran_steps)
        and str(contract.get("simulation_owner") or "").startswith("simulation:")
        and contract.get("result_simulation_owner") == contract.get("simulation_owner")
        and _digest(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def validate_ltspice_v47_identity(positive: Mapping[str, object]) -> bool:
    """Bind hierarchy/current/KCL and AC/TRAN row identities to one replay."""
    if not isinstance(positive, Mapping):
        return False
    branch = positive.get("hierarchical_branch_v47_current_sign_order_kcl_owner_identity")
    rows = positive.get("ac_tran_v47_step_tuple_trace_measure_row_identity")
    if branch is None and rows is None:
        return True
    return isinstance(branch, Mapping) and isinstance(rows, Mapping) and _branch_ok(branch) and _step_rows_ok(rows)
