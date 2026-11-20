"""Hand the CPU natives from the CPU build job to the Windows GPU wheel job.

``record`` runs on the CPU build host. It copies the top-level native files of
``src/radia`` into a fresh directory and writes ``cpu-natives.manifest.json``
with the source commit and the size and SHA-256 of every file. ``stage`` runs on
the Windows GPU host before any file is copied into the checkout. It accepts the
artifact only when its file set equals the manifest exactly, every hash
matches, the required core members are present with a clean receipt for this
commit, and no CUDA module or private build record is included.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

MANIFEST = "cpu-natives.manifest.json"
SCHEMA = "radia.ci-cpu-natives.v1"
SUFFIXES = (".pyd", ".dll", ".pyd.build.json")
BACKUP_MARKERS = (".locked-old.", ".lockedold.", ".pre_", ".backup")
CORE = "_radia_pybind.pyd"
REQUIRED = frozenset({CORE, CORE + ".build.json", "radia_motor_rom.dll"})


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _reject_name(name: str) -> None:
    lowered = name.lower()
    if lowered.startswith("_") and "_cuda" in lowered:
        raise ValueError(f"CUDA module must not come from the CPU build host: {name}")
    if lowered.endswith(".private-build.json"):
        raise ValueError(f"Private build record must not be handed over: {name}")


def _check_core(files: dict[str, Path], sha: str, hashes: dict[str, str]) -> None:
    missing = sorted(REQUIRED - set(files))
    if missing:
        raise ValueError(f"Required CPU natives missing: {missing}")
    receipt = json.loads(files[CORE + ".build.json"].read_text(encoding="utf-8"))
    if receipt.get("source_commit") != sha or receipt.get("source_dirty") is not False:
        raise ValueError(f"Core receipt is not from clean commit {sha}: "
                         f"source_commit={receipt.get('source_commit')} dirty={receipt.get('source_dirty')}")
    if receipt.get("binary_sha256") != hashes[CORE]:
        raise ValueError(f"Core receipt hash {receipt.get('binary_sha256')} != binary {hashes[CORE]}")


def record(source_dir: Path, out_dir: Path, sha: str) -> dict:
    if out_dir.exists():
        raise ValueError(f"Output directory already exists: {out_dir}")
    files = {}
    for path in sorted(source_dir.iterdir()):
        if not path.is_file() or not path.name.endswith(SUFFIXES):
            continue
        if any(marker in path.name for marker in BACKUP_MARKERS):
            continue
        _reject_name(path.name)
        files[path.name] = path
    hashes = {name: _sha256(path) for name, path in files.items()}
    _check_core(files, sha, hashes)
    out_dir.mkdir(parents=True)
    for name, path in files.items():
        shutil.copy2(path, out_dir / name)
        if _sha256(out_dir / name) != hashes[name]:
            raise ValueError(f"Copy changed {name}")
    manifest = {"schema": SCHEMA, "source_commit": sha,
                "files": [{"name": name, "bytes": files[name].stat().st_size, "sha256": hashes[name]}
                          for name in sorted(files)]}
    (out_dir / MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def stage(artifact_dir: Path, dest: Path, sha: str, evidence: Path) -> dict:
    manifest = json.loads((artifact_dir / MANIFEST).read_text(encoding="utf-8"))
    if manifest.get("schema") != SCHEMA or manifest.get("source_commit") != sha:
        raise ValueError(f"Manifest is not for commit {sha}: schema={manifest.get('schema')} "
                         f"source_commit={manifest.get('source_commit')}")
    listed = {entry["name"]: entry for entry in manifest["files"]}
    if len(listed) != len(manifest["files"]):
        raise ValueError("Manifest lists a file twice")
    present = {}
    for path in artifact_dir.rglob("*"):
        if path.is_dir():
            raise ValueError(f"Unexpected directory in CPU natives: {path.relative_to(artifact_dir)}")
        if path.name != MANIFEST:
            present[path.relative_to(artifact_dir).as_posix()] = path
    if set(present) != set(listed):
        raise ValueError(f"CPU native set differs from its manifest: unlisted={sorted(set(present) - set(listed))} "
                         f"missing={sorted(set(listed) - set(present))}")
    hashes = {}
    for name, path in present.items():
        _reject_name(name)
        if not name.endswith(SUFFIXES) or "/" in name:
            raise ValueError(f"Unexpected file in CPU natives: {name}")
        hashes[name] = _sha256(path)
        if hashes[name] != listed[name]["sha256"] or path.stat().st_size != listed[name]["bytes"]:
            raise ValueError(f"Hash mismatch for {name}: {hashes[name]} != {listed[name]['sha256']}")
    _check_core(present, sha, hashes)
    clash = sorted(p.name for p in dest.iterdir() if p.is_file() and p.name.endswith(SUFFIXES))
    if clash:
        raise ValueError(f"Destination already holds native files: {clash}")
    for name, path in present.items():
        shutil.copy2(path, dest / name)
        if _sha256(dest / name) != hashes[name]:
            raise ValueError(f"Staging changed {name}")
    evidence.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    rec = sub.add_parser("record")
    rec.add_argument("--source-dir", type=Path, required=True)
    rec.add_argument("--out", type=Path, required=True)
    rec.add_argument("--sha", required=True)
    stg = sub.add_parser("stage")
    stg.add_argument("--artifact", type=Path, required=True)
    stg.add_argument("--dest", type=Path, required=True)
    stg.add_argument("--sha", required=True)
    stg.add_argument("--evidence", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "record":
        result = record(args.source_dir, args.out, args.sha)
    else:
        result = stage(args.artifact, args.dest, args.sha, args.evidence)
    print(json.dumps({"files": len(result["files"]), "source_commit": result["source_commit"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
