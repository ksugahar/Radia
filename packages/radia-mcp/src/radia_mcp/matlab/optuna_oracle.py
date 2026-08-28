"""Checked Optuna 4.9.0 compatibility evidence for the MATLAB domain layer.

This module deliberately does not import Optuna.  Upstream Optuna and the
official optuna-mcp server remain external oracle owners; radia-mcp reads and
audits the checked evidence produced by those runtimes.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
from typing import Any


_POLICY_HEADING = "### MATLAB Optuna Upstream Differential-Oracle Policy (2026-08-21)"

# Any local function declaration, including bracketed output lists.  The name
# filter below then applies MATLAB's own rule.  Matching only ``test``-prefixed
# names missed a shared-behavior test named ``...Test`` entirely: it appeared
# in neither the discovered set nor the unclassified list, so the manifest
# gate reported "pass" while the test went unclassified.
_FUNCTION_DECLARATION = re.compile(
    r"(?mi)^[ \t]*function\s+"
    r"(?:\[[^\]]*\]\s*=\s*|[A-Za-z]\w*\s*=\s*)?"
    r"([A-Za-z]\w*)\s*\("
)


def _is_matlab_test_name(name: str) -> bool:
    """Apply MATLAB's rule: a test begins or ends with "test", any case."""
    lowered = name.lower()
    return lowered.startswith("test") or lowered.endswith("test")


def _sha256(path: Path) -> str:
    content = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(content).hexdigest().upper()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _has_evidence(candidate: Path) -> bool:
    return (
        (candidate / "matlab" / "optuna_upstream_compatibility.json").is_file()
        and (candidate / "tests" / "matlab" / "fixtures").is_dir()
    )


def _repo_root(repo_root: str | os.PathLike[str] | None = None) -> Path:
    """Resolve the Radia repository that owns the checked Optuna evidence.

    An explicitly requested root is authoritative: if it does not carry the
    evidence the caller is told so, rather than silently auditing a different
    checkout that happens to sit above this module.
    """
    argument = str(repo_root or "").strip()
    environment = str(os.environ.get("RADIA_REPO_ROOT", "")).strip()
    configured = argument or environment
    if configured:
        candidate = Path(configured)
        if not _has_evidence(candidate):
            source = "repo_root argument" if argument else "RADIA_REPO_ROOT"
            raise FileNotFoundError(
                f"The {source} {configured!r} does not contain the checked "
                "Optuna evidence (matlab/optuna_upstream_compatibility.json "
                "and tests/matlab/fixtures/). Point it at the Radia repository "
                "root, or omit it to search upward from the installed package."
            )
        return candidate.resolve()
    for candidate in Path(__file__).resolve().parents:
        if _has_evidence(candidate):
            return candidate.resolve()
    raise FileNotFoundError(
        "Could not locate the Radia repository. Set RADIA_REPO_ROOT or pass repo_root."
    )


# The checked evidence the compatibility claim rests on, named explicitly.
# A key-suffix filter used to select these and silently dropped api_coverage
# and public_api_inventory -- precisely the artifacts that back the reported
# public-API closure counts -- so the provenance hashes omitted them.
_EVIDENCE_KEYS = (
    "compatibility",
    "api_coverage",
    "python_oracle",
    "public_api_inventory",
    "mcp_oracle",
    "test_manifest",
)


def _evidence_paths(root: Path) -> dict[str, Path]:
    fixture = root / "tests" / "matlab" / "fixtures"
    return {
        "compatibility": root / "matlab" / "optuna_upstream_compatibility.json",
        "api_coverage": root / "matlab" / "optuna49_api_coverage.json",
        "python_oracle": fixture / "optuna49_oracle.json",
        "public_api_inventory": fixture / "optuna49_public_api.json",
        "mcp_oracle": fixture / "optuna49_mcp_oracle.json",
        "test_manifest": fixture / "optuna_test_manifest.json",
        "python_generator": fixture / "generate_optuna49_oracle.py",
        "api_inventory_generator": fixture / "generate_optuna49_api_inventory.py",
        "api_coverage_generator": fixture / "generate_optuna49_api_coverage.py",
        "mcp_generator": fixture / "generate_optuna49_mcp_oracle.py",
        "manifest_generator": fixture / "generate_optuna_test_manifest.py",
    }


def matlab_optuna_compatibility_contract(
    repo_root: str | os.PathLike[str] | None = None,
) -> dict[str, Any]:
    """Return the checked supported/integration/unsupported Optuna boundary."""
    root = _repo_root(repo_root)
    paths = _evidence_paths(root)
    compatibility = _json(paths["compatibility"])
    api_coverage = _json(paths["api_coverage"])
    python_oracle = _json(paths["python_oracle"])
    public_api = _json(paths["public_api_inventory"])
    mcp_oracle = _json(paths["mcp_oracle"])
    manifest = _json(paths["test_manifest"])
    counts = Counter(
        str(entry.get("classification", "")) for entry in manifest.get("entries", [])
    )
    errors: list[str] = []
    if compatibility.get("schema") != "radia.optuna-upstream-compatibility.v1":
        errors.append("unexpected compatibility schema")
    claim = str(compatibility.get("claim", ""))
    if not any(
        boundary in claim
        for boundary in ("not a drop-in replacement", "not a Python binary drop-in")
    ):
        errors.append("compatibility contract must reject binary drop-in parity")
    if python_oracle.get("optuna_version") != "4.9.0":
        errors.append("direct Python oracle is not optuna==4.9.0")
    if mcp_oracle.get("optuna_version") != "4.9.0":
        errors.append("official MCP oracle is not backed by optuna==4.9.0")
    if mcp_oracle.get("transport") != "stdio":
        errors.append("official MCP oracle was not captured over stdio")
    if manifest.get("upstream_version") != "4.9.0":
        errors.append("test manifest does not pin Optuna 4.9.0")
    if public_api.get("schema") != "radia.test.optuna49-public-api.v1":
        errors.append("unexpected Optuna public API inventory schema")
    if public_api.get("optuna_version") != "4.9.0":
        errors.append("public API inventory is not pinned to Optuna 4.9.0")
    if api_coverage.get("schema") != "radia.optuna49-api-coverage.v1":
        errors.append("unexpected MATLAB Optuna API coverage schema")
    if api_coverage.get("upstream_version") != "4.9.0":
        errors.append("MATLAB API coverage is not pinned to Optuna 4.9.0")
    if api_coverage.get("upstream_inventory_sha256") != _sha256(
        paths["public_api_inventory"]
    ):
        errors.append("MATLAB API coverage was not generated from the checked inventory")
    if bool(api_coverage.get("full_compatibility_complete")):
        if (
            int(api_coverage.get("surface_missing_count", -1)) != 0
            or int(api_coverage.get("oracle_unmapped_count", -1)) != 0
            or int(api_coverage.get("oracle_partial_count", -1)) != 0
        ):
            errors.append("full API compatibility was claimed before closure")
    return {
        "schema": "radia-mcp.matlab-optuna-compatibility/v1",
        "status": "ready" if not errors else "error",
        "ok": not errors,
        "runtime_owner": "MathWorks MATLAB MCP Server",
        "oracle_owner": "optuna==4.9.0 and official optuna/optuna-mcp",
        "claim": claim,
        "repo_root": str(root),
        "oracle_versions": {
            "optuna": python_oracle.get("optuna_version"),
            "python": python_oracle.get("python_version"),
            "numpy": python_oracle.get("numpy_version"),
            "scipy": python_oracle.get("scipy_version"),
            "torch": python_oracle.get("torch_version"),
            "cmaes": python_oracle.get("cmaes_version"),
            "optuna_mcp": mcp_oracle.get("optuna_mcp_version"),
        },
        "transport": {
            "seeded_numerical_oracle": "direct-python",
            "public_mcp_contract": mcp_oracle.get("transport"),
            "mcp_sampler_seed_supported": bool(
                mcp_oracle.get("sampler_seed_supported", False)
            ),
        },
        "test_counts": {
            "upstream_python": counts["upstream-python"],
            "upstream_mcp": counts["upstream-mcp"],
            "matlab_integration": counts["matlab-integration"],
            "total": sum(counts.values()),
        },
        "shared_behavior_verified": compatibility.get(
            "shared_behavior_verified", []
        ),
        "matlab_integration_only": compatibility.get(
            "matlab_integration_only", []
        ),
        "unsupported_or_not_yet_oracled": compatibility.get(
            "unsupported_or_not_yet_oracled", []
        ),
        "public_api_closure": {
            key: api_coverage.get(key)
            for key in (
                "surface_entry_count",
                "surface_present_count",
                "surface_missing_count",
                "oracle_verified_count",
                "oracle_partial_count",
                "oracle_unmapped_count",
                "full_compatibility_complete",
            )
        },
        "evidence_sha256": {
            name: _sha256(paths[name]) if paths[name].is_file() else None
            for name in _EVIDENCE_KEYS
        },
        "errors": errors,
    }


def _policy_block(path: Path) -> str:
    source = path.read_text(encoding="utf-8")
    start = source.find(_POLICY_HEADING)
    if start < 0:
        return ""
    end = source.find("\n---\n", start)
    return source[start:] if end < 0 else source[start:end]


def matlab_optuna_oracle_audit(
    repo_root: str | os.PathLike[str] | None = None,
) -> dict[str, Any]:
    """Audit fixture provenance, policy sync, and test-manifest completeness."""
    root = _repo_root(repo_root)
    contract = matlab_optuna_compatibility_contract(root)
    paths = _evidence_paths(root)
    manifest = _json(paths["test_manifest"])
    expected = {
        (str(entry.get("file")), str(entry.get("test")))
        for entry in manifest.get("entries", [])
    }
    discovered: set[tuple[str, str]] = set()
    for path in sorted((root / "tests" / "matlab").glob("test_optuna*.m")):
        for name in _FUNCTION_DECLARATION.findall(path.read_text(encoding="utf-8")):
            if name != path.stem and _is_matlab_test_name(name):
                discovered.add((path.name, name))

    agent_policy = _policy_block(root / "AGENTS.md")
    claude_policy = _policy_block(root / "CLAUDE.md")
    missing = sorted(discovered - expected)
    stale = sorted(expected - discovered)
    errors = list(contract["errors"])
    if missing:
        errors.append(f"unclassified MATLAB tests: {missing}")
    if stale:
        errors.append(f"stale manifest entries: {stale}")
    if not agent_policy or agent_policy != claude_policy:
        errors.append("AGENTS.md and CLAUDE.md Optuna policies are not synchronized")
    missing_files = [str(path) for path in paths.values() if not path.is_file()]
    if missing_files:
        errors.append(f"missing evidence files: {missing_files}")

    return {
        "schema": "radia-mcp.matlab-optuna-oracle-audit/v1",
        "status": "pass" if not errors else "fail",
        "ok": not errors,
        "repo_root": str(root),
        "upstream_version": "4.9.0",
        "test_function_count": len(discovered),
        "manifest_entry_count": len(expected),
        "missing_manifest_entries": [list(item) for item in missing],
        "stale_manifest_entries": [list(item) for item in stale],
        "policy_identical": bool(agent_policy and agent_policy == claude_policy),
        "policy_sha256": (
            hashlib.sha256(agent_policy.encode("utf-8")).hexdigest().upper()
            if agent_policy
            else None
        ),
        "evidence_sha256": {
            name: _sha256(path) for name, path in paths.items() if path.is_file()
        },
        "regeneration_commands": [
            "python tests/matlab/fixtures/generate_optuna49_oracle.py",
            "python tests/matlab/fixtures/generate_optuna49_api_inventory.py",
            "python tests/matlab/fixtures/generate_optuna49_api_coverage.py",
            "python tests/matlab/fixtures/generate_optuna49_mcp_oracle.py",
            "python tests/matlab/fixtures/generate_optuna_test_manifest.py",
        ],
        "matlab_verification": {
            "execute_with": "official MATLAB MCP run_matlab_test_file",
            "test_glob": "tests/matlab/test_optuna*.m",
        },
        "errors": errors,
    }


__all__ = [
    "matlab_optuna_compatibility_contract",
    "matlab_optuna_oracle_audit",
]
