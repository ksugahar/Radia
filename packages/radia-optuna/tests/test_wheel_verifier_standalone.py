"""Detached-package wheel verification contract.

Guards: Reintroducing monorepo-root reads into wheel verification, or accepting
a package whose upstream oracle pin drifts despite correct wheel metadata.
Assumes: A copy of the package source alone, with a synthetic PE header standing
in for the MEX. This checks packaging behavior, not native execution or ABI.
"""

import importlib.util
import json
from pathlib import Path
import shutil
import zipfile

import pytest


def test_detached_package_verifies_without_cpu_project(tmp_path):
    source = Path(__file__).resolve().parents[1]
    detached = tmp_path / "standalone"
    shutil.copytree(source, detached, ignore=shutil.ignore_patterns(
        "__pycache__", "*.mexw64", "build", "dist", "*.egg-info"))
    assert not (detached.parents[1] / "pyproject.toml").exists()
    spec = importlib.util.spec_from_file_location(
        "detached_optuna_verifier", detached / "verify_wheel.py")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    manifest = json.loads(
        (detached / "src/radia_optuna/manifest.json").read_text(encoding="utf-8"))
    pe = bytearray(1024)
    pe[:2] = b"MZ"
    pe[60:64] = (64).to_bytes(4, "little")
    pe[64:70] = b"PE\0\0\x64\x86"
    (detached / "matlab/optuna_mex.mexw64").write_bytes(pe)
    wheel = tmp_path / f"radia_optuna-{manifest['version']}-py3-none-win_amd64.whl"
    dist_info = f"radia_optuna-{manifest['version']}.dist-info"
    with zipfile.ZipFile(wheel, "w") as archive:
        for member, payload in verifier._source_payloads(manifest).items():
            archive.write(payload, member)
        archive.writestr(f"{dist_info}/METADATA", (
            "Metadata-Version: 2.1\nName: radia-optuna\n"
            f"Version: {manifest['version']}\n"
            'Requires-Dist: optuna==5.0.0; extra == "upstream"\n\n'))
        archive.writestr(f"{dist_info}/WHEEL", "Tag: py3-none-win_amd64\n")
    result = verifier.verify(wheel)
    assert result["ok"] is True, result
    assert result["source_fidelity_verified"] is True, result

    metadata = detached / "pyproject.toml"
    text = metadata.read_text(encoding="utf-8")
    assert '"optuna==5.0.0"' in text
    metadata.write_text(text.replace('"optuna==5.0.0"', '"optuna==4.9.0"'),
                        encoding="utf-8")
    with pytest.raises(SystemExit) as failure:
        verifier.verify(wheel)
    assert failure.value.code == 1, failure.value.code
