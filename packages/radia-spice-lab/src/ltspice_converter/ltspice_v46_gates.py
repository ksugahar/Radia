"""Replay identities for partial LTspice transient and noise artifacts."""

from __future__ import annotations

import math
from collections.abc import Mapping


def _digest(value: object) -> bool:
    text = str(value or "")
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text.lower())


def _equal(contract: Mapping[str, object], *names: str) -> bool:
    return all(contract.get(f"result_{name}") == contract.get(name) for name in names)


def _finite_nonnegative(value: object) -> bool:
    try:
        return math.isfinite(float(value)) and float(value) >= 0.0
    except (TypeError, ValueError):
        return False


def _buck_ok(contract: Mapping[str, object]) -> bool:
    return (
        bool(str(contract.get("generation_id") or ""))
        and contract.get("result_generation_id") == contract.get("generation_id")
        and _equal(
            contract,
            "failed_timestep_index",
            "raw_truncated",
            "raw_point_count",
            "expected_raw_point_count",
            "measure_status",
            "nonfinite_count",
            "analysis_status",
            "waveform_owner",
            "raw_sha256",
        )
        and contract.get("raw_truncated") is False
        and contract.get("raw_point_count") == contract.get("expected_raw_point_count")
        and contract.get("measure_status") == "complete"
        and contract.get("analysis_status") == "completed"
        and _finite_nonnegative(contract.get("failed_timestep_index"))
        and _finite_nonnegative(contract.get("raw_point_count"))
        and contract.get("nonfinite_count") == 0
        and _digest(contract.get("raw_sha256"))
        and _digest(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _noise_ok(contract: Mapping[str, object]) -> bool:
    return (
        bool(str(contract.get("generation_id") or ""))
        and contract.get("result_generation_id") == contract.get("generation_id")
        and _equal(
            contract,
            "random_seed",
            "sample_count",
            "sample_filter",
            "psd_bin_hz",
            "partial_sweep",
            "sweep_status",
            "nonfinite_count",
            "measure_owner",
            "raw_sha256",
        )
        and isinstance(contract.get("random_seed"), int)
        and contract.get("sample_count") == 64
        and contract.get("sample_filter") == "discard_first_8"
        and _finite_nonnegative(contract.get("psd_bin_hz"))
        and contract.get("partial_sweep") is False
        and contract.get("sweep_status") == "complete"
        and contract.get("nonfinite_count") == 0
        and _digest(contract.get("raw_sha256"))
        and _digest(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def validate_ltspice_v46_identity(positive: Mapping[str, object]) -> bool:
    """Require complete RAW/measure identity for held-out transient/noise replays."""
    if not isinstance(positive, Mapping):
        return False
    buck = positive.get("buck_v46_failed_timestep_raw_truncation_measure_nan_inf_identity")
    noise = positive.get("noise_v46_monte_carlo_seed_sample_filter_psd_bin_partial_identity")
    if buck is None and noise is None:
        return True
    return isinstance(buck, Mapping) and isinstance(noise, Mapping) and _buck_ok(buck) and _noise_ok(noise)
