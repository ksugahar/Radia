"""Generate the checked classification of every MATLAB Optuna test."""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TEST_DIRECTORY = ROOT / "tests" / "matlab"

CLASSIFICATION_BY_FILE = {
    "test_optuna_upstream_oracle.m": "upstream-python",
    "test_optuna_core_parity.m": "upstream-python",
    "test_optuna_pruner_parity.m": "upstream-python",
    "test_optuna_mcp_oracle.m": "upstream-mcp",
    "test_optuna_table.m": "matlab-integration",
    "test_optuna_simulink_block.m": "matlab-integration",
    "test_optuna_sampler_wrappers.m": "matlab-integration",
    "test_optuna_reliability.m": "matlab-integration",
    "test_optuna_nsgaii_joint.m": "matlab-integration",
}

INTEGRATION_SCOPE = {
    "test_optuna_table.m": (
        "MATLAB table storage, MAT persistence, and SimulinkRunner integration"
    ),
    "test_optuna_simulink_block.m": "MATLAB and Simulink block integration",
    "test_optuna_sampler_wrappers.m": (
        "MATLAB sampler-state persistence and AutoSampler policy integration"
    ),
    "test_optuna_reliability.m": (
        "MATLAB MAT-file transaction, recovery, and sampler-state persistence"
    ),
    "test_optuna_nsgaii_joint.m": (
        "MATLAB NSGA-II generation cache, parallel scheduling, and persistence"
    ),
}


def test_functions(path: Path) -> list[str]:
    source = path.read_text(encoding="utf-8")
    names = re.findall(
        r"(?m)^function\s+(?:\w+\s*=\s*)?(test\w+)\s*\(", source
    )
    return [name for name in names if name != path.stem]


def build_manifest() -> dict[str, object]:
    files = sorted(path.name for path in TEST_DIRECTORY.glob("test_optuna*.m"))
    missing = sorted(set(files) - CLASSIFICATION_BY_FILE.keys())
    stale = sorted(CLASSIFICATION_BY_FILE.keys() - set(files))
    if missing or stale:
        raise RuntimeError(f"Unclassified files: {missing}; stale files: {stale}")

    entries: list[dict[str, object]] = []
    for filename in files:
        classification = CLASSIFICATION_BY_FILE[filename]
        for test_name in test_functions(TEST_DIRECTORY / filename):
            effective_classification = classification
            if test_name == "testOraclePolicyManifestIsComplete":
                effective_classification = "matlab-integration"
            if effective_classification == "upstream-python":
                oracle = "optuna49_oracle.json"
                scope = "shared Optuna 4.9.0 behavior"
            elif effective_classification == "upstream-mcp":
                oracle = "optuna49_mcp_oracle.json"
                scope = "official Optuna MCP behavior"
            else:
                oracle = None
                scope = INTEGRATION_SCOPE.get(
                    filename, "MATLAB test-policy integration"
                )
            entries.append(
                {
                    "classification": effective_classification,
                    "file": filename,
                    "oracle": oracle,
                    "scope": scope,
                    "test": test_name,
                }
            )

    return {
        "schema": "radia.test.optuna-matlab-policy.v1",
        "upstream_version": "4.9.0",
        "entries": entries,
    }


def main() -> None:
    destination = Path(__file__).with_name("optuna_test_manifest.json")
    destination.write_text(
        json.dumps(build_manifest(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(destination)


if __name__ == "__main__":
    main()
