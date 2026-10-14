"""Reject unregistered ESRF6 meshes before a candidate solve.

This identity gate does not certify discretization accuracy. A new mesh family
requires a reviewed contract update, not a filename change or a force option.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


DEFAULT_CONTRACT = Path(__file__).with_name("case6_mesh_contract.json")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_inputs(iron: Path, fem: Path, contract_path: Path = DEFAULT_CONTRACT) -> dict:
    contract_bytes = contract_path.read_bytes()
    contract = json.loads(contract_bytes)
    if contract.get("schema") != "radia.validation.case6-mesh-contract.v1" or contract.get("case") != 6:
        raise ValueError("not an ESRF6 mesh identity contract")
    verified = {}
    for role, path in (("iron", iron), ("fem", fem)):
        expected = contract[role]["sha256"]
        if not isinstance(expected, str) or len(expected) != 64 or any(c not in "0123456789abcdef" for c in expected):
            raise ValueError(f"invalid {role} SHA-256 in mesh contract")
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"{role} mesh identity mismatch: expected {expected}, got {actual}: {path}")
        verified[role] = {"path": str(path.resolve()), "sha256": actual}
    return {
        "schema": "radia.validation.case6-mesh-identity.v1",
        "identity_passed": True,
        "numerical_acceptance": "HOLD: identity is not numerical acceptance",
        "contract_sha256": hashlib.sha256(contract_bytes).hexdigest(),
        "meshes": verified,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iron", type=Path, required=True)
    parser.add_argument("--fem", type=Path, required=True)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify_inputs(args.iron, args.fem, args.contract)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
