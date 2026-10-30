"""Validate independent CPU Python and GPU integration ownership."""

from __future__ import annotations

import argparse
import fnmatch
import json
from collections import Counter
from pathlib import Path


def audit(repo_root: Path) -> dict:
    repo_root = repo_root.resolve()
    manifest_path = repo_root / "packages/radia-simulink/matlab" / "python_api_parity_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    source_root = repo_root / manifest["source_root"]
    allowed = set(manifest["classifications"])
    rules = manifest["rules"]
    errors: list[str] = []
    assignments: list[dict] = []
    if manifest.get("schema") != "radia.cpu-gpu-ownership/v2":
        errors.append("expected the CPU/GPU ownership contract")
    if manifest.get("cpu_mex_supported") is not False or manifest.get("implicit_cpu_fallback") is not False:
        errors.append("CPU MEX and implicit CPU fallback must remain unsupported")
    expected_owners = {"cpu": "radia", "matlab": "radia-simulink", "simulink": "radia-simulink",
                       "ltspice": "radia-simulink", "mcp": "radia-mcp", "optuna": "radia-optuna"}
    if manifest.get("owners") != expected_owners:
        errors.append("integration distribution owners differ from policy")
    rule_ids = [rule["id"] for rule in rules]
    if len(rule_ids) != len(set(rule_ids)):
        errors.append("rule ids must be unique")
    for rule in rules:
        if rule.get("classification") not in allowed:
            errors.append(f"{rule['id']}: unknown classification")
        if rule.get("matlab"):
            errors.append(f"{rule['id']}: CPU classification cannot promise a MATLAB implementation")
    for feature in manifest.get("gpu_features", []):
        if feature.get("cpu_fallback") is not False:
            errors.append(f"{feature['id']}: implicit CPU fallback is forbidden")
        for field in ("source", "matlab", "validation"):
            relative = feature.get(field, "")
            resolved = (repo_root / relative).resolve()
            if not relative or not resolved.is_relative_to(repo_root) or not resolved.is_file():
                errors.append(f"{feature['id']}: missing or escaping {field}: {relative}")
            elif field == "source" and not relative.startswith(("packages/radia-simulink/native/", "src/core/gpu/")):
                errors.append(f"{feature['id']}: {field} is not GPU-owned")
            elif field == "matlab" and not relative.startswith("packages/radia-simulink/matlab/"):
                errors.append(f"{feature['id']}: {field} is not GPU-owned")
    for entry in manifest.get("binary_extensions", []):
        if entry.get("classification") != "cpu-python-only" or entry.get("matlab"):
            errors.append(f"{entry['python']}: CPU binary must not advertise CPU MEX parity")

    for path in sorted(source_root.rglob("*.py")):
        relative = path.relative_to(source_root).as_posix()
        matches = [
            rule
            for rule in rules
            if any(fnmatch.fnmatchcase(relative, pattern) for pattern in rule["patterns"])
        ]
        if not matches:
            errors.append(f"unclassified Python module: {relative}")
            continue
        rule = matches[0]
        assignments.append(
            {
                "python": relative,
                "rule": rule["id"],
                "classification": rule["classification"],
                "matlab": rule.get("matlab", []),
            }
        )

    counts = Counter(item["classification"] for item in assignments)
    return {
        "schema": manifest["schema"],
        "ok": not errors,
        "manifest": str(manifest_path),
        "python_file_count": len(list(source_root.rglob("*.py"))),
        "classified_file_count": len(assignments),
        "counts": dict(sorted(counts.items())),
        "gpu_features": manifest.get("gpu_features", []),
        "binary_extensions": manifest.get("binary_extensions", []),
        "assignments": assignments,
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--compact", action="store_true")
    args = parser.parse_args()
    result = audit(args.repo_root)
    print(json.dumps(result, indent=None if args.compact else 2, ensure_ascii=True))
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
