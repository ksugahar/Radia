"""Build/staging contracts for Radia-owned Python CUDA payloads, not runtime acceptance."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import shutil

MODULES = ("_gpu_csr", "_gpu_resident_amg", "_gpu_ams", "_hdiv_cuda")
WINDOWS_ARCHITECTURES = ("75", "80")
DIAGNOSTIC_FLAGS = ("fault_injection", "device_fault_injection", "lifecycle_diagnostic")


def verify_payload(binary_name, data, record, source_commit, architectures=WINDOWS_ARCHITECTURES):
    """Verify compiled artifact identity; a receipt does not prove a GPU solve."""
    if (record.get("schema") != "radia.native-build-provenance.v1"
            or record.get("source_commit") != source_commit
            or not re.fullmatch(r"[0-9a-f]{40}", source_commit)
            or record.get("source_dirty") is not False
            or record.get("binary_name") != binary_name
            or record.get("binary_sha256") != hashlib.sha256(data).hexdigest()
            or record.get("architectures") != list(architectures)
            or record.get("ptx_architectures") != []
            or any(record.get(k) is not False for k in DIAGNOSTIC_FLAGS)):
        raise ValueError("Python GPU artifact source/hash/architecture/build-mode mismatch: " + binary_name)
    for field in ("source_hashes", "pybind_header_hashes"):
        values = record.get(field)
        if (not isinstance(values, dict) or not values
                or any(not isinstance(k, str) or not isinstance(v, str)
                       or not re.fullmatch(r"[0-9a-f]{64}", v)
                       or Path(k).is_absolute() or ".." in Path(k).parts for k, v in values.items())):
            raise ValueError("Python GPU source/header closure is missing or malformed: " + binary_name)
    return record


def stage(output, package, source_commit, suffix):
    """Copy only verified owned binaries and public sidecars, never runtime DLLs."""
    output, package = Path(output), Path(package)
    verified = []
    for module in MODULES:
        binary = output / (module + suffix)
        receipt = Path(str(binary) + ".build.json")
        if binary.is_symlink() or receipt.is_symlink():
            raise ValueError("Python GPU staging refuses linked payloads")
        record = json.loads(receipt.read_text(encoding="utf-8"))
        verify_payload(binary.name, binary.read_bytes(), record, source_commit)
        verified.append((binary, receipt))
    # Fail before copying anything if a member or its provenance is invalid.
    for pair in verified:
        for path in pair:
            shutil.copyfile(path, package / path.name)
    return [p.name for pair in verified for p in pair]


if __name__ == "__main__":
    import argparse
    import sysconfig
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    args = parser.parse_args()
    print(json.dumps(stage(args.output, args.package, args.source_commit,
                           sysconfig.get_config_var("EXT_SUFFIX"))))
