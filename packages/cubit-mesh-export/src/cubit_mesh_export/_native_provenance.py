"""Content-addressed provenance for cubit-mesh-export native payloads."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


SCHEMA = "cubit-mesh-export.native-payloads.v2"
HASH_FORMAT = "sha256-path-nul-lf-content-nul-v1"
REQUIRED_PAYLOADS = ("cubit_mesh_export.ccm", "cubit_mesh_curver.pyd")
SOURCE_SUFFIXES = {
    ".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx",
    ".cmake", ".txt", ".ps1",
}
SKIP_DIRS = {"build", "build-ccm", "build-pyd", "__pycache__"}


def native_source_files(repo_root: Path) -> list[Path]:
    """Return the deterministic set of files that can affect native builds."""

    repo_root = repo_root.resolve()
    plugin_root = repo_root / "src" / "cubit_plugin"
    files = [
        path for path in plugin_root.rglob("*")
        if path.is_file()
        and not any(part in SKIP_DIRS for part in path.relative_to(plugin_root).parts)
        and path.suffix.lower() in SOURCE_SUFFIXES
    ]
    shared_discovery = repo_root / "tools" / "find_cubit.ps1"
    if shared_discovery.is_file():
        files.append(shared_discovery)
    return sorted(files, key=lambda path: path.relative_to(repo_root).as_posix())


def native_source_digest(repo_root: Path) -> tuple[str, int]:
    """Hash normalized paths and bytes, never filesystem timestamps."""

    repo_root = repo_root.resolve()
    files = native_source_files(repo_root)
    digest = hashlib.sha256()
    digest.update((HASH_FORMAT + "\0").encode("ascii"))
    for path in files:
        relative = path.relative_to(repo_root).as_posix().encode("utf-8")
        digest.update(relative)
        digest.update(b"\0")
        # Git may check text files out as LF or CRLF on different Windows
        # hosts. Line endings do not change the compiled program, so hash a
        # canonical LF representation to keep provenance portable.
        content = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
        digest.update(content)
        digest.update(b"\0")
    return digest.hexdigest(), len(files)


def _file_evidence(path: Path) -> tuple[str, int]:
    data = path.read_bytes()
    return hashlib.sha256(data).hexdigest(), len(data)


def _source_commit(repo_root: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo_root), "log", "-1", "--format=%H", "--",
         "src/cubit_plugin", "tools/find_cubit.ps1"],
        check=True, capture_output=True, text=True,
    )
    commit = result.stdout.strip()
    if len(commit) != 40:
        raise RuntimeError("could not resolve the native source commit")
    return commit


def verify_manifest(repo_root: Path | None, package_dir: Path) -> list[str]:
    """Verify payloads; omit source comparison only for source-free sdists."""

    manifest_path = package_dir / "native_payloads.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    errors: list[str] = []
    if manifest.get("schema") != SCHEMA:
        errors.append(f"manifest schema must be {SCHEMA}")

    source = manifest.get("source", {})
    if source.get("hash_format") != HASH_FORMAT:
        errors.append(f"source hash_format must be {HASH_FORMAT}")
    if repo_root is not None:
        actual_source_hash, actual_file_count = native_source_digest(repo_root)
        if source.get("tree_sha256") != actual_source_hash:
            errors.append(
                "native source content differs from manifest: "
                f"actual={actual_source_hash} recorded={source.get('tree_sha256')}"
            )
        if source.get("file_count") != actual_file_count:
            errors.append(
                "native source file count differs from manifest: "
                f"actual={actual_file_count} recorded={source.get('file_count')}"
            )
    if not source.get("commit"):
        errors.append("native source commit is not recorded")

    payloads = manifest.get("payloads", {})
    for name in REQUIRED_PAYLOADS:
        path = package_dir / name
        expected = payloads.get(name)
        if not isinstance(expected, dict):
            errors.append(f"manifest has no evidence for {name}")
            continue
        if not path.is_file():
            errors.append(f"required payload is missing: {path}")
            continue
        actual_hash, actual_size = _file_evidence(path)
        if expected.get("sha256") != actual_hash:
            errors.append(
                f"{name} SHA-256 differs: actual={actual_hash} "
                f"recorded={expected.get('sha256')}"
            )
        if expected.get("size") != actual_size:
            errors.append(
                f"{name} size differs: actual={actual_size} "
                f"recorded={expected.get('size')}"
            )
    return errors


def record_manifest(repo_root: Path, package_dir: Path) -> dict:
    """Record current native source and payload evidence after a real build."""

    manifest_path = package_dir / "native_payloads.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    source_hash, file_count = native_source_digest(repo_root)
    manifest["schema"] = SCHEMA
    manifest["source"] = {
        "hash_format": HASH_FORMAT,
        "tree_sha256": source_hash,
        "file_count": file_count,
        "commit": _source_commit(repo_root),
    }
    payloads = manifest.setdefault("payloads", {})
    for name in REQUIRED_PAYLOADS:
        path = package_dir / name
        if not path.is_file():
            raise FileNotFoundError(f"required payload is missing: {path}")
        sha256, size = _file_evidence(path)
        payload = payloads.setdefault(name, {})
        payload["sha256"] = sha256
        payload["size"] = size
        payload.setdefault("platform", "win_amd64")
        payload.setdefault("cubit_version", "2025.12")
        payload.setdefault("toolchain", "MSVC 19.50 / CMake / Ninja")
        if name.endswith(".pyd"):
            payload["asset_name"] = f"cubit_mesh_curver-{sha256}.pyd"
            payload.setdefault("python_abi", "cp312")
            payload.setdefault("netgen_version", "6.2.2606")
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("verify", "record"))
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--package-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.action == "record":
        manifest = record_manifest(args.repo_root, args.package_dir)
        print(
            "recorded native provenance: "
            f"source={manifest['source']['tree_sha256']}"
        )
        return 0
    errors = verify_manifest(args.repo_root, args.package_dir)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print("native source and payload provenance OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
