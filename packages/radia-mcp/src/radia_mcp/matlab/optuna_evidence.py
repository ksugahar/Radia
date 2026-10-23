"""Read package-owned Optuna evidence; never execute upstream producers."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath


def checked_evidence(package_root: Path) -> dict[str, Path]:
    """Resolve a source package or wheel and reject missing/altered payloads."""
    root = package_root.resolve()
    metadata = root if (root / "manifest.json").is_file() else root / "src/radia_optuna"
    distribution = json.loads((metadata / "manifest.json").read_text(encoding="utf-8"))
    name = distribution["compatibility_evidence"]
    if name != "evidence_manifest.json":
        raise ValueError("Unsupported compatibility evidence manifest location")
    manifest_path = metadata / name
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest.get("schema") != "radia-optuna.compatibility-evidence.v1"
            or manifest.get("version") != distribution["version"]
            or manifest.get("oracle_version") != distribution["oracle_version"]
            or manifest.get("hash_basis") != "utf8-lf-bytes"):
        raise ValueError("Invalid Optuna evidence identity/hash basis")
    records = manifest.get("files", {})
    required = {"compatibility", "api_coverage", "python_oracle", "public_api_inventory",
                "mcp_oracle", "test_manifest", "cuda_oracle", "cuda_density_oracle"}
    if set(records) != required:
        raise ValueError("Incomplete Optuna evidence inventory")
    paths = {}
    for key, record in records.items():
        relative = PurePosixPath(record["path"])
        if (len(relative.parts) != 2 or relative.parts[0] not in {"oracles", "matlab"}
                or relative.suffix != ".json" or ".." in relative.parts
                or "\\" in str(relative) or ":" in str(relative)):
            raise ValueError(f"Unsafe Optuna evidence path: {relative}")
        path = (root / str(relative)).resolve()
        if not path.is_relative_to(root):
            raise ValueError(f"Optuna evidence escapes its package: {relative}")
        actual = hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        if actual != record["sha256"]:
            raise ValueError(f"Optuna evidence hash mismatch: {key}")
        paths[key] = path
    paths["evidence_manifest"] = manifest_path
    return paths


def installed_package_root() -> Path | None:
    """Locate the installed package without importing Optuna or its backend."""
    spec = importlib.util.find_spec("radia_optuna")
    if spec is None or not spec.origin:
        return None
    metadata = Path(spec.origin).resolve().parent
    if metadata.parent.name == "src" and (metadata.parents[1] / "pyproject.toml").is_file():
        return None  # Source is handled as source, never as a wheel.
    return metadata
