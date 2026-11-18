"""Shared solver acceptance roles; physical identities are private configuration.

Exact-wheel release acceptance runs on the validation role only. The CPU build
host prepares the release candidate and is not an acceptance target.
"""
import os

RELEASE_ACCEPTANCE_HOSTS = ("validation",)
PUBLIC_RELEASE_ACCEPTANCE_ROLES = {"validation": "validation"}
CURRENT_RELEASE_ACCEPTANCE_HOSTNAMES = {
    role: os.environ.get("RADIA_ACCEPTANCE_HOSTNAME_" + role.upper(), "").strip().lower()
    for role in RELEASE_ACCEPTANCE_HOSTS
}
RELEASE_ACCEPTANCE_HOSTNAMES = {
    role: frozenset({name}) if name else frozenset()
    for role, name in CURRENT_RELEASE_ACCEPTANCE_HOSTNAMES.items()
}

# Each release lane is built and validated by a fixed role on the same source.
RELEASE_LANES = {
    "cpu": ("cpu-build", "windows-validation"),
    "linux-python-gpu": ("linux-gpu", "linux-gpu"),
    "windows-python-gpu": ("windows-gpu", "windows-gpu"),
    "matlab": ("windows-gpu", "windows-gpu"),
}
REQUIRED_VALIDATION_CHECKS = {
    "linux-python-gpu": ("cpu_gpu_parity",),
    "windows-python-gpu": ("cpu_gpu_parity",),
    "matlab": ("simulink_acceptance", "ltspice_acceptance"),
}
# The published radia wheel is the Windows Python GPU lane artifact.
PUBLISHED_WHEEL_LANE = "windows-python-gpu"
RELEASE_LANES_SUMMARY = "release-lanes.json"
RELEASE_LANES_SUMMARY_SCHEMA = "radia.release-lanes.v1"


def verify_lane_receipts(data, source_sha, read):
    """Check the four-lane receipt manifest; ``read(relative_path)`` returns bytes.

    Returns the validated artifact hashes per lane. Raises ValueError on any
    missing receipt, role, source, check or hash mismatch.
    """
    import hashlib
    import re
    from pathlib import PurePosixPath

    def receipt(entry):
        rel = entry["path"]
        parts = PurePosixPath(rel).parts
        if not rel or "\\" in rel or rel.startswith("/") or ".." in parts or ":" in rel:
            raise ValueError(f"Receipt path must be relative inside the receipt directory: {rel}")
        digest = entry["sha256"]
        if not re.fullmatch(r"[0-9a-f]{64}", digest) or hashlib.sha256(read(rel)).hexdigest() != digest:
            raise ValueError(f"SHA-256 mismatch for {rel}")
        return digest

    if data["schema"] != 1 or data["source_sha"] != source_sha:
        raise ValueError("Release manifest schema or source SHA mismatch")
    native_ci = data["native_ci"]
    if native_ci["source_sha"] != source_sha or native_ci["passed"] is not True:
        raise ValueError("Exact-source native CI receipt must pass")
    receipt(native_ci["evidence"])
    lanes = data["lanes"]
    if set(lanes) != set(RELEASE_LANES):
        raise ValueError("A receipt is required for each of the four release lanes")
    validated = {}
    for name, (build_role, validation_role) in RELEASE_LANES.items():
        build, validation = lanes[name]["build"], lanes[name]["validation"]
        if build["role"] != build_role or validation["role"] != validation_role:
            raise ValueError(f"Incorrect build or validation role for {name}")
        if build["source_sha"] != source_sha or validation["source_sha"] != source_sha:
            raise ValueError(f"Source SHA mismatch for {name}")
        receipt(build["evidence"])
        if validation["passed"] is not True:
            raise ValueError(f"Validation did not pass for {name}")
        receipt(validation["evidence"])
        for check in REQUIRED_VALIDATION_CHECKS.get(name, ()):
            receipt(validation["checks"][check])
        artifacts = build["artifacts"]
        if not artifacts:
            raise ValueError(f"No artifacts for {name}")
        actual = {entry["path"]: receipt(entry) for entry in artifacts}
        if len(actual) != len(artifacts) or actual != validation["artifact_sha256"]:
            raise ValueError(f"Validated artifact set or hash mismatch for {name}")
        validated[name] = actual
    return validated


def lane_summary(data, validated):
    """Public digest of verified receipts: roles, checks and artifact hashes, no evidence."""
    from pathlib import PurePosixPath

    lanes = {}
    for name, (build_role, validation_role) in RELEASE_LANES.items():
        artifacts = {}
        for path, digest in validated[name].items():
            base = PurePosixPath(path).name
            if base in artifacts:
                raise ValueError(f"Duplicate artifact name in {name}: {base}")
            artifacts[base] = digest
        lanes[name] = {"build_role": build_role, "validation_role": validation_role, "passed": True,
                       "checks": sorted(data["lanes"][name]["validation"].get("checks", {})),
                       "artifact_sha256": artifacts}
    return {"schema": RELEASE_LANES_SUMMARY_SCHEMA, "source_sha": data["source_sha"],
            "native_ci_passed": True, "lanes": lanes}


def verify_lane_summary(summary, source_sha):
    """Check a committed public lane summary; returns the artifact hashes per lane."""
    import re

    if summary.get("schema") != RELEASE_LANES_SUMMARY_SCHEMA or summary.get("source_sha") != source_sha:
        raise ValueError("Release lane summary schema or source SHA mismatch")
    if summary.get("native_ci_passed") is not True:
        raise ValueError("Exact-source native CI did not pass")
    lanes = summary.get("lanes", {})
    if set(lanes) != set(RELEASE_LANES):
        raise ValueError("A summary is required for each of the four release lanes")
    result = {}
    for name, (build_role, validation_role) in RELEASE_LANES.items():
        lane = lanes[name]
        if lane.get("build_role") != build_role or lane.get("validation_role") != validation_role:
            raise ValueError(f"Incorrect build or validation role for {name}")
        if lane.get("passed") is not True:
            raise ValueError(f"Validation did not pass for {name}")
        missing = set(REQUIRED_VALIDATION_CHECKS.get(name, ())) - set(lane.get("checks", []))
        if missing:
            raise ValueError(f"Required checks missing for {name}: {sorted(missing)}")
        artifacts = lane.get("artifact_sha256", {})
        if not artifacts or not all(re.fullmatch(r"[0-9a-f]{64}", str(v)) for v in artifacts.values()):
            raise ValueError(f"No valid artifact hashes for {name}")
        result[name] = artifacts
    return result
