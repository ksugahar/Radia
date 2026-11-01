"""Detached-package wheel and coverage verification contract.

Guards: Reintroducing monorepo-root reads into wheel verification or API
coverage generation, changing frozen oracle inputs during verification, or accepting
a package whose upstream oracle pin drifts despite correct wheel metadata, or
missing/tampered installed evidence and accidentally bundled producers.
Assumes: A copy of the package source alone, with a synthetic PE header standing
in for the MEX, and stored upstream inputs. This checks packaging and coverage
metadata behavior, not native execution or ABI.
"""

import hashlib
import importlib.util
import os
import json
from pathlib import Path
import shutil
import subprocess
import sys
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
    tag = f"cp{sys.version_info.major}{sys.version_info.minor}"
    extension = f"_optuna_kernels.{tag}-win_amd64.pyd"
    build = {"archive_sha256":"a"*64,"objects":{"optuna_cuda.cu.obj":"b"*64},
             "python_extension_sha256s":{extension:hashlib.sha256(pe).hexdigest()},
             "mex_sha256":hashlib.sha256(pe).hexdigest()}
    (detached / "src/radia_optuna/native_build.json").write_text(json.dumps(build),encoding="utf-8")
    wheel = tmp_path / f"radia_optuna-{manifest['version']}-{tag}-{tag}-win_amd64.whl"
    dist_info = f"radia_optuna-{manifest['version']}.dist-info"
    with zipfile.ZipFile(wheel, "w") as archive:
        for member, payload in verifier._source_payloads(manifest).items():
            archive.write(payload, member)
        archive.writestr(f"{dist_info}/METADATA", (
            "Metadata-Version: 2.1\nName: radia-optuna\n"
            f"Version: {manifest['version']}\n"
            'Requires-Dist: optuna==5.0.0; extra == "upstream"\n\n'))
        archive.writestr(f"{dist_info}/WHEEL", f"Tag: {tag}-{tag}-win_amd64\n")
        archive.writestr(f"radia_optuna/{extension}", pe)
    result = verifier.verify(wheel)
    assert result["ok"] is True, result
    assert result["source_fidelity_verified"] is True, result

    assert result["compatibility_evidence_verified"] is True, result
    evidence = json.loads((detached / "src/radia_optuna/evidence_manifest.json").read_text(encoding="utf-8"))
    with zipfile.ZipFile(wheel) as archive:
        original = {name: archive.read(name) for name in archive.namelist()}
    target = "radia_optuna/" + evidence["files"]["python_oracle"]["path"]
    # Exercise the public verifier with damaged distributions, not mocked hashes.
    for damage in ("missing", "tampered", "producer"):
        members = dict(original)
        if damage == "missing":
            del members[target]
        elif damage == "tampered":
            members[target] += b" "
        else:
            members["radia_optuna/oracles/generate_optuna50_oracle.py"] = b"raise RuntimeError('must not ship')\n"
        with zipfile.ZipFile(wheel, "w") as archive:
            for name, payload in members.items():
                archive.writestr(name, payload)
        with pytest.raises(SystemExit) as failure:
            verifier.verify(wheel)
        assert failure.value.code == 1, damage
    with zipfile.ZipFile(wheel, "w") as archive:
        for name, payload in original.items():
            archive.writestr(name, payload)
    assert verifier.verify(wheel)["ok"] is True

    # Run the real package-local entry point without a monorepo or source-path
    # rescue. Only this detached copy's coverage ledger may be reproduced;
    # upstream producers and stored numerical expectations are read-only inputs.
    coverage = detached / "matlab/optuna50_api_coverage.json"
    expected_coverage = coverage.read_bytes().replace(b"\r\n", b"\n")
    inputs = {p: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (detached / "oracles").iterdir() if p.is_file()}
    environment = dict(os.environ)
    environment.pop("PYTHONPATH", None)
    environment.pop("RADIA_REPO_ROOT", None)
    completed = subprocess.run(
        [sys.executable, str(detached / "tools/generate_optuna50_api_coverage.py")],
        cwd=detached, env=environment, capture_output=True, text=True, timeout=60,
    )
    assert completed.returncode == 0, (completed.stdout, completed.stderr)
    assert coverage.read_bytes().replace(b"\r\n", b"\n") == expected_coverage
    assert {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs} == inputs

    metadata = detached / "pyproject.toml"
    text = metadata.read_text(encoding="utf-8")
    assert '"optuna==5.0.0"' in text
    metadata.write_text(text.replace('"optuna==5.0.0"', '"optuna==4.9.0"'),
                        encoding="utf-8")
    with pytest.raises(SystemExit) as failure:
        verifier.verify(wheel)
    assert failure.value.code == 1, failure.value.code
