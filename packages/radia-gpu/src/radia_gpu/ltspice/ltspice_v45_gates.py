"""Release-bound replay identities for the frozen v45 LTSpice lane."""

from __future__ import annotations

from collections.abc import Mapping

from ._artifact_identity import (
    digest_is_sha256 as _digest,
)


def _release_owner_ok(value: Mapping[str, object], owner_key: str) -> bool:
    release = str(value.get("release_id") or "")
    owner = str(value.get(owner_key) or "")
    return (
        bool(release)
        and value.get("result_release_id") == release
        and bool(owner)
        and value.get(f"accepted_{owner_key}") == owner
        and value.get(f"result_{owner_key}") == owner
        and _digest(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def validate_ltspice_v45_identity(positive: Mapping[str, object]) -> bool:
    """Require release, owner, and digest closure for public replay identities."""
    if not isinstance(positive, Mapping):
        return False
    buck = positive.get(
        "buck_v45_startup_softstart_inductor_current_ripple_efficiency_energy_waveform_identity"
    )
    noise = positive.get(
        "noise_v45_ac_transfer_psd_sidedness_bandwidth_correlation_measure_identity"
    )
    if buck is None and noise is None:
        return True
    if not isinstance(buck, Mapping) or not isinstance(noise, Mapping):
        return False
    return _release_owner_ok(buck, "waveform_owner") and _release_owner_ok(
        noise, "measure_owner"
    )
