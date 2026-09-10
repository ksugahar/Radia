"""Explicit, non-networked maintenance; serving never emits diagnostic JSON."""
from __future__ import annotations

import argparse
import importlib
import importlib.metadata as metadata
import hashlib
import json
import os
import re
from pathlib import Path
import sys
import subprocess
import uuid
from urllib.parse import urlsplit
from urllib.request import url2pathname

BASELINE = {
    "radia-meta": "meta",
    "radia-build123d": "build123d",
    "radia-cubit": "cubit",
    "radia-gmsh": "gmsh",
    "radia-analysis": "radia-analysis",
    "radia-motion": "radia-motion",
    "radia-publication": "radia-publication",
}


def server_module(key: str) -> str:
    from .meta.catalog import CATALOG
    if key not in CATALOG:
        raise ValueError(f"Unknown catalog server: {key}")
    return CATALOG[key]["subpackage"] + ".server"


def standard_entry(key: str, python: str) -> dict:
    server_module(key)  # Validate without importing domain/optional dependencies.
    executable = Path(python)
    if not executable.is_absolute() or not executable.is_file():
        raise ValueError("Python must be an existing absolute executable path")
    return {"command": str(executable), "args": [
        "-s", "-m", "radia_mcp.maintenance", "serve", key,
    ]}


def _migratable(entry: dict, key: str, expected: dict) -> bool:
    # Retain env/cwd/approvals/timeouts. Unknown args, wrappers and remote
    # transports require review, rather than silently widening tool access.
    if entry.get("url") or entry.get("type", "stdio") != "stdio":
        return False
    command = str(entry.get("command", "")).replace("\\", "/").split("/")[-1].lower()
    if not re.fullmatch(r"python(?:\d+(?:\.\d+)*)?(?:\.exe)?", command):
        return False
    args = entry.get("args", [])
    module = server_module(key)
    return args in (expected["args"], ["-m", module], ["-s", "-m", module])


def _unique_object(pairs: list[tuple]) -> dict:
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("Duplicate JSON key")
        result[key] = value
    return result


def plan_config(path: Path, python: str) -> tuple[bytes | None, bytes, dict]:
    """Prepare an all-or-nothing change; never expose client secrets in report."""
    original = path.read_bytes() if path.exists() else None
    text = (original or b"").decode("utf-8-sig")
    toml = path.suffix.lower() == ".toml"
    if toml:
        try:
            import tomlkit
        except ImportError as exc:
            raise ValueError("Install radia-mcp[maintenance] for TOML configuration") from exc
        document = tomlkit.parse(text)
        section = "mcp_servers"
    else:
        document = json.loads(text, object_pairs_hook=_unique_object) if text.strip() else {}
        section = "mcpServers"
    if not isinstance(document, dict):
        raise ValueError("Client configuration must be an object/table")
    if section not in document:
        document[section] = {}
    entries = document[section]
    if not hasattr(entries, "items"):
        raise ValueError(f"{section} must be an object/table")
    actions, conflicts = [], []
    for name, key in BASELINE.items():
        expected = standard_entry(key, python)
        if name not in entries:
            entries[name] = expected
            actions.append({"server": name, "action": "add"})
        else:
            entry = entries[name]
            if not hasattr(entry, "get") or not _migratable(entry, key, expected):
                conflicts.append(name)
                continue
            if all(entry.get(k) == v for k, v in expected.items()):
                continue
            for field, value in expected.items():
                entry[field] = value
            actions.append({"server": name, "action": "normalize-launcher"})
    if toml:
        rendered = tomlkit.dumps(document)
    else:
        rendered = json.dumps(document, ensure_ascii=False, indent=2) + "\n"
    newline = "\r\n" if original and b"\r\n" in original else "\n"
    rendered = rendered.replace("\r\n", "\n").replace("\n", newline)
    candidate = rendered.encode("utf-8")
    if original and original.startswith(b"\xef\xbb\xbf"):
        candidate = b"\xef\xbb\xbf" + candidate
    # True no-op: preserve exact bytes, whitespace and comments on repeat runs.
    if not actions or conflicts:
        candidate = original or b""
    return original, candidate, {"path": str(path), "actions": actions,
        "conflicts": conflicts, "changed": candidate != (original or b""),
        "status": "conflict" if conflicts else "ready"}


def apply_config(path: Path, original: bytes | None, candidate: bytes, *,
                 owner: str = "", reason: str = "", clients_idle: bool = False,
                 targets: str = "") -> str | None:
    """Apply only within a recorded, exclusive maintenance operation."""
    if candidate == (original or b""):
        return _write_config(path, original, candidate)
    from . import _maintenance_guard as guard
    if not targets.strip():
        raise ValueError("Affected client targets must be named")
    with guard.change(owner=owner, reason=reason, clients_idle=clients_idle,
                      scope={"operation": "config", "path": str(path.resolve()),
                             "targets": targets},
                      expected={"sha256": _digest(original)},
                      proposed={"sha256": _digest(candidate)}) as record:
        backup = _write_config(path, original, candidate)
        record["observed"] = {"sha256": _digest(path.read_bytes()), "backup": backup}
    return backup


def _digest(value: bytes | None) -> str | None:
    return hashlib.sha256(value).hexdigest() if value is not None else None


def _write_config(path: Path, original: bytes | None, candidate: bytes) -> str | None:
    """Back up exact bytes, detect intervening edits, preserve existing ACLs.

    Stop clients editing this file before applying. This is not a distributed
    transaction or a lock against non-cooperating client writers.
    """
    current = path.read_bytes() if path.exists() else None
    if current != original:
        raise ValueError("Configuration changed since planning; re-run the plan")
    if candidate == (original or b""):
        return None
    backup = None
    if original is not None:
        backup = path.with_name(path.name + ".before-radia-" + uuid.uuid4().hex)
        # Keep backups in the same protected directory. Never print their content.
        with backup.open("xb") as stream:
            stream.write(original)
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = "r+b" if original is not None else "xb"
    with path.open(mode) as stream:
        if original is not None and stream.read() != original:
            raise ValueError("Configuration changed before write; re-run the plan")
        stream.seek(0)
        stream.write(candidate)
        stream.truncate()
        stream.flush()
        os.fsync(stream.fileno())
    if path.read_bytes() != candidate:
        raise OSError(f"Configuration verification failed; backup: {backup}")
    return str(backup) if backup else None


def doctor(expected_root: Path | None = None, expected_version: str | None = None,
           expected_commit: str | None = None) -> dict:
    import radia_mcp
    distribution = metadata.distribution("radia-mcp")
    direct = json.loads(distribution.read_text("direct_url.json") or "{}")
    loaded = Path(radia_mcp.__file__).resolve().parent
    issues = []
    if not direct.get("dir_info", {}).get("editable"):
        issues.append("not-editable")
    source_url = urlsplit(direct.get("url", ""))
    editable_source = None
    if source_url.scheme == "file":
        local_url_path = ("//" + source_url.netloc if source_url.netloc else "") + source_url.path
        editable_source = Path(url2pathname(local_url_path)).resolve()
        if loaded not in (editable_source / "src" / "radia_mcp", editable_source / "radia_mcp"):
            issues.append("editable-metadata-source-mismatch")
    if not sys.flags.no_user_site:
        issues.append("user-site-not-isolated: run Python with -s")
    if expected_root is not None and loaded != expected_root.resolve():
        issues.append("loaded-source-mismatch")
    if expected_version is not None and distribution.version != expected_version:
        issues.append("installed-version-mismatch")
    revision = None
    repository = next((p for p in loaded.parents if (p / ".git").exists()), None)
    if repository is not None:
        try:
            result = subprocess.run(["git", "-c", f"safe.directory={repository.as_posix()}",
                "-C", str(repository), "rev-parse", "HEAD"],
                capture_output=True, text=True, timeout=10, check=True)
            revision = result.stdout.strip()
        except (OSError, subprocess.SubprocessError):
            pass
    if expected_commit is not None and revision != expected_commit:
        issues.append("loaded-commit-mismatch-or-unavailable")
    return {"schema": "radia-mcp.maintenance.v1", "status": "fail" if issues else "pass",
        "python": sys.executable, "version": distribution.version,
        "loaded_package": str(loaded), "editable_source": str(editable_source) if editable_source else None,
        "loaded_commit": revision,
        "issues": issues, "scope": "this new process only; running clients require reconnect"}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=metadata.version("radia-mcp"))
    commands = parser.add_subparsers(dest="mode", required=True)
    serve = commands.add_parser("serve", help="Run a catalog server without launcher wrappers")
    serve.add_argument("server")
    serve.add_argument("server_args", nargs=argparse.REMAINDER)
    check = commands.add_parser("doctor", help="Check this interpreter, not other live sessions")
    check.add_argument("--expected-root", type=Path, help="Exact src/radia_mcp directory")
    check.add_argument("--expected-version")
    check.add_argument("--expected-commit", help="Approved full Git SHA, not a moving branch")
    config = commands.add_parser("config", help="Plan by default; never remove existing servers")
    config.add_argument("path", type=Path)
    config.add_argument("--python", default=sys.executable)
    config.add_argument("--apply", action="store_true")
    activate = commands.add_parser("activate-editable", help="Explicit radia-mcp-only activation")
    activate.add_argument("source", type=Path, help="Approved radia-mcp package directory")
    activate.add_argument("--commit", required=True, help="Exact approved full Git SHA")
    activate.add_argument("--expected-source", type=Path, required=True)
    activate.add_argument("--expected-commit", required=True)
    recovery = commands.add_parser("reconcile", help="Record audited recovery; never retry a write")
    recovery.add_argument("change_id")
    for mutation in (config, activate, recovery):
        mutation.add_argument("--owner", default="")
        mutation.add_argument("--reason", default="")
        if mutation is not recovery:
            mutation.add_argument("--targets", default="", help="Affected host/user/client/server scope")
        mutation.add_argument("--clients-idle", action="store_true",
                              help="Operator confirms all affected consumers are idle")
    args = parser.parse_args(argv)
    try:
        if args.mode == "serve":
            module = server_module(args.server)
            old_argv = sys.argv
            try:
                sys.argv = [module, *args.server_args]
                result = importlib.import_module(module).main()
                return result if isinstance(result, int) else 0
            finally:
                sys.argv = old_argv
        if args.mode == "doctor":
            report = doctor(args.expected_root, args.expected_version, args.expected_commit)
        elif args.mode == "activate-editable":
            from ._editable_activation import activate_editable
            report = activate_editable(args.source, args.commit, args.expected_source,
                                       args.expected_commit, owner=args.owner, reason=args.reason,
                                       clients_idle=args.clients_idle, targets=args.targets)
        elif args.mode == "reconcile":
            from . import _maintenance_guard as guard
            from ._editable_activation import observe_scope
            report = guard.reconcile(args.change_id, owner=args.owner, reason=args.reason,
                                     clients_idle=args.clients_idle, observe=observe_scope)
        else:
            original, candidate, report = plan_config(args.path, args.python)
            if args.apply and not report["conflicts"]:
                report["backup"] = apply_config(
                    args.path, original, candidate, owner=args.owner, reason=args.reason,
                    clients_idle=args.clients_idle, targets=args.targets,
                )
                report["applied"] = report["changed"]
        print(json.dumps(report, ensure_ascii=True, indent=2))
        return 1 if report["status"] in {"fail", "conflict"} else 0
    except (ValueError, OSError, subprocess.SubprocessError, metadata.PackageNotFoundError) as exc:
        # Parser/IO error strings may contain configuration values: keep logs clean.
        print(f"Maintenance failed ({type(exc).__name__}); inspect the input locally.", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
