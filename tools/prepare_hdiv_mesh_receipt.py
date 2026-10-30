"""Build-host producer using the authoritative installed check-vol distribution.

This runs only the existing strict-label/default-quality checker. It does not
approve its own output: an independent artifact manifest must pin the resulting
receipt digest and provenance before a consumer can use them.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def distribution_identity():
    dist = metadata.distribution("cubit-mesh-export")
    if dist.version.lower() in ("", "unknown") or json.loads(
            dist.read_text("direct_url.json") or "{}").get("dir_info", {}).get("editable", False):
        raise ValueError("Receipt producer requires an installed, non-editable checker distribution")
    records = [f for f in dist.files or () if str(f).endswith(".dist-info/RECORD")]
    if len(records) != 1:
        raise ValueError("Checker distribution RECORD is missing or ambiguous")
    for file in dist.files:
        if file.hash:
            data = Path(dist.locate_file(file)).read_bytes()
            actual = base64.urlsafe_b64encode(hashlib.new(file.hash.mode, data).digest()).decode().rstrip("=")
            if actual != file.hash.value:
                raise ValueError("Checker installed RECORD mismatch: " + str(file))
    return dict(distribution="cubit-mesh-export", version=dist.version,
                record_sha256=hashlib.sha256(Path(dist.locate_file(records[0])).read_bytes()).hexdigest())


def prepare(mesh_path, label_contract, output):
    mesh_path, label_contract, output = Path(mesh_path), Path(label_contract), Path(output)
    outputs = [output, Path(str(output) + ".report.json"), Path(str(output) + ".stderr.txt"),
               Path(str(output) + ".provenance.json")]
    if any(path.exists() for path in outputs):
        raise FileExistsError("Receipt producer never overwrites existing evidence")
    mesh_bytes, contract_bytes = mesh_path.read_bytes(), label_contract.read_bytes()
    suffix = ".vol.gz" if mesh_path.name.lower().endswith(".vol.gz") else ".vol"
    if not mesh_path.name.lower().endswith((".vol", ".vol.gz")):
        raise ValueError("Expected .vol or .vol.gz input")
    checker = distribution_identity()
    ng_version = metadata.version("ngsolve")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="radia-check-vol-") as directory:
        mesh = Path(directory) / ("mesh" + suffix)
        contract = Path(directory) / "labels.json"
        mesh.write_bytes(mesh_bytes); contract.write_bytes(contract_bytes)
        argv = [sys.executable, "-m", "cubit_mesh_export.check", str(mesh),
                "--strict-labels", "--contract", str(contract), "--format", "json"]
        env = dict(os.environ, PYTHONNOUSERSITE="1", PYTHONSAFEPATH="1",
                   PYTHONUTF8="1", PYTHONIOENCODING="utf-8")
        env.pop("PYTHONPATH", None); env.pop("PYTHONHOME", None)
        result = subprocess.run(argv, cwd=directory, env=env, capture_output=True, check=False)
        outputs[1].write_bytes(result.stdout); outputs[2].write_bytes(result.stderr)
        if result.returncode != 0:
            raise RuntimeError("Authoritative check-vol failed; raw report/stderr retained")
        report = json.loads(result.stdout)
        if report.get("passed") is not True:
            raise ValueError("Authoritative check-vol did not pass")
        if checker != distribution_identity() or ng_version != metadata.version("ngsolve"):
            raise ValueError("Checker distribution changed during validation")
        if mesh.read_bytes() != mesh_bytes or contract.read_bytes() != contract_bytes:
            raise ValueError("Check-vol inputs changed during execution")
        receipt = dict(schema="radia.checked-vol-receipt.v1", checker=checker,
                       ngsolve_version=ng_version, python_version=sys.version,
                       argv=argv, returncode=result.returncode,
                       mesh_sha256=hashlib.sha256(mesh_bytes).hexdigest(),
                       contract_sha256=hashlib.sha256(contract_bytes).hexdigest(),
                       report_sha256=hashlib.sha256(result.stdout).hexdigest(),
                       report_base64=base64.b64encode(result.stdout).decode())
    with output.open("x", encoding="utf-8") as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False)
    provenance = {key: receipt[key] for key in ("checker", "argv", "ngsolve_version")}
    with outputs[3].open("x", encoding="utf-8") as stream:
        json.dump(provenance, stream, indent=2, allow_nan=False)
    return dict(receipt_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                provenance=provenance, approval_status="requires independent artifact-manifest approval")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--label-contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.mesh, args.label_contract, args.output), indent=2))
