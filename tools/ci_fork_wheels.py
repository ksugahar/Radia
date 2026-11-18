"""Install and verify the locked Windows fork wheels in an isolated CI venv."""
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

LOCK = Path(__file__).with_name("ci_fork_wheels.json")


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def verify_installed(lock):
    site = Path(sys.prefix).resolve()
    for name, wheel in lock["wheels"].items():
        if importlib.metadata.version(name) != wheel["version"]:
            raise RuntimeError(f"Incorrect installed version: {name}")
    packages = {}
    for name in ("netgen", "ngsolve"):
        spec = importlib.util.find_spec(name)
        if spec is None or not spec.submodule_search_locations:
            raise RuntimeError(f"Missing package: {name}")
        directories = [Path(path).resolve() for path in spec.submodule_search_locations]
        if any(not path.is_relative_to(site) for path in directories):
            raise RuntimeError(f"{name} is outside the CI environment: {directories}")
        packages[name] = directories[0]
    for name, sha in lock["dlls"].items():
        if digest(packages["netgen"] / name) != sha:
            raise RuntimeError(f"Incorrect installed DLL: {name}")
    return lock["sources"]


def install(lock, directory):
    if sys.prefix == sys.base_prefix:
        raise RuntimeError("Install requires an isolated CI venv")
    with tempfile.TemporaryDirectory(prefix="radia-fork-wheels-") as tmp:
        paths = []
        for wheel in lock["wheels"].values():
            filename = wheel["filename"]
            path = Path(tmp) / filename
            shutil.copyfile(directory / filename, path)
            if digest(path) != wheel["sha256"]:
                raise RuntimeError(f"Incorrect wheel SHA-256: {filename}")
            paths.append(str(path))
        subprocess.run([sys.executable, "-I", "-m", "pip", "install", "--only-binary=:all:", *paths], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("install", "verify"))
    args = parser.parse_args()
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    if args.action == "install":
        directory = os.environ.get("RADIA_FORK_WHEEL_DIR")
        if not directory:
            raise RuntimeError("RADIA_FORK_WHEEL_DIR must point to the provisioned fork wheel directory")
        install(lock, Path(directory))
    print(json.dumps(verify_installed(lock)))


if __name__ == "__main__":
    main()
