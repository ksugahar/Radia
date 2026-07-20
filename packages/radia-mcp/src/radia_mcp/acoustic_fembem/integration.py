"""Radia-owned integration for MathWorks' MATLAB MCP Server."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path
from typing import Any, Mapping


_HERE = Path(__file__).resolve().parent
_EXTENSION_FILE = _HERE / "extensions" / "radia-acoustic-fembem-tools.json"
_PROJECT_ROOT_ENV = "ACOUSTIC_FEMBEM_ROOT"
_SERVER_ENV = "RADIA_MATLAB_MCP_SERVER"
_PROFILES = {
    "existing": ("existing", "desktop"),
    "auto_nodesktop": ("auto", "nodesktop"),
    "new_nodesktop": ("new", "nodesktop"),
}


def acoustic_fembem_extension_path() -> Path:
    """Return the installed custom-tool extension file."""
    if not _EXTENSION_FILE.is_file():
        raise FileNotFoundError(f"MATLAB MCP extension is missing: {_EXTENSION_FILE}")
    return _EXTENSION_FILE


def acoustic_fembem_extension_contract() -> dict[str, Any]:
    """Load and validate the Radia-owned MATLAB custom-tool contract."""
    path = acoustic_fembem_extension_path()
    raw = path.read_bytes()
    manifest = json.loads(raw.decode("utf-8"))
    errors = _extension_errors(manifest)
    tools = manifest.get("tools", [])
    signatures = manifest.get("signatures", {})
    return {
        "schema": "radia-mcp.acoustic-fembem-extension/v1",
        "status": "ok" if not errors else "error",
        "ok": not errors,
        "runtime_owner": "MathWorks MATLAB MCP Server",
        "domain_owner": "radia_mcp.acoustic_fembem",
        "extension_file": str(path),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "tool_count": len(tools) if isinstance(tools, list) else 0,
        "tool_names": [tool.get("name", "") for tool in tools]
        if isinstance(tools, list)
        else [],
        "signature_count": len(signatures) if isinstance(signatures, dict) else 0,
        "errors": errors,
    }


def acoustic_fembem_server_config(
    project_root: str | os.PathLike[str] | None = None,
    profile: str = "existing",
) -> dict[str, Any]:
    """Build official-server arguments for the acoustic FEM-BEM extension.

    ``project_root`` is the companion MATLAB solver checkout. Radia owns the
    MCP contract; the GPL teaching solver remains outside the BSD wheel.
    """
    if profile not in _PROFILES:
        raise ValueError(
            f"Unknown profile {profile!r}; expected one of {sorted(_PROFILES)}"
        )

    contract = acoustic_fembem_extension_contract()
    if not contract["ok"]:
        raise RuntimeError(f"Invalid MATLAB MCP extension: {contract['errors']}")

    root = _resolve_project_root(project_root)
    entrypoints: list[str] = []
    if root is not None:
        entrypoints = _companion_entrypoints()
        missing = [name for name in entrypoints if not _matlab_function_path(root, name).is_file()]
        if missing:
            raise ValueError(
                "Companion MATLAB project does not implement the packaged MCP contract: "
                + ", ".join(missing)
            )
    extension_path = acoustic_fembem_extension_path()
    setup_code = _matlab_setup_code(root) if root is not None else ""
    session_mode, display_mode = _PROFILES[profile]
    args = [
        f"--matlab-session-mode={session_mode}",
        f"--extension-file={extension_path}",
    ]
    if display_mode == "nodesktop":
        args.append("--matlab-display-mode=nodesktop")
    if root is not None:
        args.append(f"--initial-working-folder={root}")

    preflight = [
        "Install or update the official MathWorks MATLAB MCP Server.",
        "Make the companion acoustic_fembem MATLAB package visible before using custom tools.",
        "Keep this readable education solver separate from the radia.acoustics production API.",
    ]
    if profile == "existing":
        preflight.extend(
            [
                "Run shareMATLABSession() in the MATLAB session to reuse.",
                "Run matlab_setup_code in that shared session before calling a domain tool.",
            ]
        )
    else:
        preflight.append(
            "Run matlab_setup_code with evaluate_matlab_code before the first domain tool call."
        )

    return {
        "schema": "radia-mcp.acoustic-fembem-server-config/v1",
        "status": "ok",
        "runtime_owner": "MathWorks MATLAB MCP Server",
        "domain_owner": "radia_mcp.acoustic_fembem",
        "solver_owner": "companion MATLAB education project",
        "profile": profile,
        "command_id": "matlab-mcp-server",
        "command": _official_server_command(),
        "command_env": _SERVER_ENV,
        "args": args,
        "extension_file": str(extension_path),
        "project_root": str(root) if root is not None else "",
        "project_root_env": _PROJECT_ROOT_ENV,
        "matlab_setup_code": setup_code,
        "matlab_entrypoints": entrypoints,
        "preflight": preflight,
        "tool_count": contract["tool_count"],
        "tool_names": contract["tool_names"],
    }


def _extension_errors(manifest: Any) -> list[str]:
    if not isinstance(manifest, Mapping):
        return ["root must be an object"]

    tools = manifest.get("tools")
    signatures = manifest.get("signatures")
    errors: list[str] = []
    if not isinstance(tools, list) or not tools:
        errors.append("tools must be a non-empty array")
        tools = []
    if not isinstance(signatures, Mapping):
        errors.append("signatures must be an object")
        signatures = {}

    names: list[str] = []
    for index, tool in enumerate(tools):
        if not isinstance(tool, Mapping):
            errors.append(f"tools[{index}] must be an object")
            continue
        name = tool.get("name")
        if not isinstance(name, str) or not name:
            errors.append(f"tools[{index}].name must be a non-empty string")
            continue
        names.append(name)
        schema = tool.get("inputSchema")
        if not isinstance(schema, Mapping) or schema.get("type") != "object":
            errors.append(f"{name}: inputSchema must be an object schema")
        annotations = tool.get("annotations")
        if not isinstance(annotations, Mapping) or not annotations.get("readOnlyHint"):
            errors.append(f"{name}: readOnlyHint must be true")

    if len(names) != len(set(names)):
        errors.append("tool names must be unique")
    if set(names) != set(signatures):
        errors.append("tool names and signature names must match")

    for name, signature in signatures.items():
        if not isinstance(signature, Mapping):
            errors.append(f"{name}: signature must be an object")
            continue
        function = signature.get("function")
        if not isinstance(function, str) or not function.startswith("acoustic_fembem."):
            errors.append(f"{name}: function must use the acoustic_fembem package")
        order = signature.get("input", {}).get("order")
        if not isinstance(order, list) or not all(isinstance(item, str) for item in order):
            errors.append(f"{name}: input.order must be a string array")
    return errors


def _resolve_project_root(
    project_root: str | os.PathLike[str] | None,
) -> Path | None:
    raw = str(project_root or "").strip() or os.environ.get(_PROJECT_ROOT_ENV, "").strip()
    if not raw:
        return None
    root = Path(raw).expanduser().resolve()
    required = [root / "+acoustic_fembem", root / "matlab_api"]
    missing = [str(path.name) for path in required if not path.is_dir()]
    if missing:
        raise ValueError(
            f"Not an acoustic FEM-BEM MATLAB project: {root}; missing {', '.join(missing)}"
        )
    return root


def _companion_entrypoints() -> list[str]:
    manifest = json.loads(acoustic_fembem_extension_path().read_text(encoding="utf-8"))
    return [
        signature["function"]
        for signature in manifest["signatures"].values()
    ]


def _matlab_function_path(root: Path, qualified_name: str) -> Path:
    package_name, function_name = qualified_name.split(".", 1)
    return root / f"+{package_name}" / f"{function_name}.m"


def _matlab_setup_code(root: Path) -> str:
    quoted = str(root).replace("'", "''")
    return (
        f"root = '{quoted}'; "
        "addpath(root); "
        "addpath(genpath(fullfile(root, 'matlab_api'))); "
        "addpath(fullfile(root, 'validation')); "
        f"setenv('{_PROJECT_ROOT_ENV}', root);"
    )


def _official_server_command() -> str:
    configured = os.environ.get(_SERVER_ENV, "").strip()
    if configured:
        return configured

    toolkit_candidate = (
        Path.home()
        / ".matlab"
        / "agentic-toolkits"
        / "bin"
        / "matlab-mcp-server-windows-x64.exe"
    )
    if toolkit_candidate.is_file():
        return str(toolkit_candidate)

    candidates = [
        "matlab-mcp-server",
        "matlab-mcp-server-windows-x64.exe",
        "matlab-mcp-core-server",
        "matlab-mcp-core-server-win64.exe",
    ]
    for candidate in candidates:
        resolved = shutil.which(candidate)
        if resolved:
            return resolved

    return "matlab-mcp-server"


def main(argv: list[str] | None = None) -> int:
    """Print a machine-readable official MATLAB MCP server configuration."""
    parser = argparse.ArgumentParser(
        description="Print Radia's MATLAB acoustic FEM-BEM MCP configuration."
    )
    parser.add_argument("--project-root", default="")
    parser.add_argument("--profile", choices=sorted(_PROFILES), default="existing")
    parser.add_argument(
        "--contract-only",
        action="store_true",
        help="Print the packaged extension contract instead of server arguments.",
    )
    args = parser.parse_args(argv)
    payload = (
        acoustic_fembem_extension_contract()
        if args.contract_only
        else acoustic_fembem_server_config(args.project_root, args.profile)
    )
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


__all__ = [
    "main",
    "acoustic_fembem_extension_contract",
    "acoustic_fembem_extension_path",
    "acoustic_fembem_server_config",
]
