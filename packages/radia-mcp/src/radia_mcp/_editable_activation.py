"""Opt-in radia-mcp-only activation; never restart or uninstall live services."""

from __future__ import annotations

import ast
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import url2pathname

from . import _maintenance_guard as guard


def _git(source: Path, *args: str) -> str:
    repository = next((p for p in (source, *source.parents) if (p / ".git").exists()), None)
    if repository is None:
        raise ValueError("Source must belong to an approved Git snapshot")
    return subprocess.run(
        ["git", "-c", f"safe.directory={repository.as_posix()}", "-C", str(repository), *args],
        capture_output=True,
        text=True,
        encoding="utf-8",
        check=True,
        timeout=30,
    ).stdout.strip()


def _snapshot(source: Path, commit: str, *, clean: bool) -> dict:
    source = source.resolve(strict=True)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("An exact full Git SHA is required")
    if not (source / "src/radia_mcp/__init__.py").is_file():
        raise ValueError("Source is not a radia-mcp package")
    actual = _git(source, "rev-parse", "HEAD")
    if actual != commit:
        raise ValueError("Snapshot commit differs from the approved commit")
    if clean and _git(source, "status", "--porcelain", "--untracked-files=normal"):
        raise ValueError("Candidate snapshot has uncommitted changes")
    tree = ast.parse((source / "src/radia_mcp/__init__.py").read_text(encoding="utf-8"))
    versions = [
        node.value.value
        for node in tree.body
        if isinstance(node, ast.Assign)
        and isinstance(node.value, ast.Constant)
        and isinstance(node.value.value, str)
        and any(
            isinstance(target, ast.Name) and target.id == "__version__" for target in node.targets
        )
    ]
    if len(versions) != 1:
        raise ValueError("Candidate package version must be explicitly declared")
    try:
        import tomllib

        manifest = tomllib.loads((source / "pyproject.toml").read_text(encoding="utf-8"))
    except ImportError:
        try:
            import tomlkit
        except ImportError as exc:
            raise ValueError("Python 3.10 activation requires radia-mcp[maintenance]") from exc
        manifest = tomlkit.parse((source / "pyproject.toml").read_text(encoding="utf-8"))
    project = manifest.get("project", {})
    if project.get("name") != "radia-mcp" or project.get("version") != versions[0]:
        raise ValueError("Only a consistent radia-mcp distribution may be activated")
    return {"source": guard.identity(source), "commit": actual, "version": versions[0]}


def installed() -> dict:
    # -I excludes task-local PYTHONPATH/user-site. This is intentionally only
    # fresh-process evidence, not a claim about any already running client.
    script = (
        "import importlib.metadata as m,importlib.util as u,json;"
        "d=m.distribution('radia-mcp');s=u.find_spec('radia_mcp');"
        "print(json.dumps({'version':d.version,'direct':json.loads(d.read_text('direct_url.json') or '{}'),"
        "'module':s.origin if s else None}))"
    )
    result = subprocess.run(
        [sys.executable, "-I", "-c", script],
        capture_output=True,
        text=True,
        encoding="utf-8",
        check=True,
        timeout=30,
    )
    data = json.loads(result.stdout)
    direct = data["direct"]
    parsed = urlsplit(direct.get("url", ""))
    if not direct.get("dir_info", {}).get("editable") or parsed.scheme != "file":
        raise ValueError("Existing installation must be editable with known provenance")
    source = Path(url2pathname(("//" + parsed.netloc if parsed.netloc else "") + parsed.path))
    if not data["module"] or guard.identity(Path(data["module"])) != guard.identity(
        source / "src/radia_mcp/__init__.py"
    ):
        raise ValueError("Fresh-process resolution differs from editable registration")
    return {
        "source": guard.identity(source),
        "commit": _git(source, "rev-parse", "HEAD"),
        "version": data["version"],
        "scope": "fresh-process-only",
    }


def _pip(source: Path) -> None:
    # No separate uninstall window, dependency upgrades, arbitrary package or
    # unbounded retry. Child errors are not copied into the shared receipt.
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    subprocess.run(
        [
            sys.executable,
            "-I",
            "-m",
            "pip",
            "--isolated",
            "install",
            "--no-deps",
            "--no-build-isolation",
            "--no-cache-dir",
            "--force-reinstall",
            "--prefix",
            sys.prefix,
            "-e",
            str(source),
        ],
        env=env,
        capture_output=True,
        check=True,
        timeout=300,
    )


def activate_editable(
    source: Path,
    commit: str,
    expected_source: Path,
    expected_commit: str,
    *,
    owner: str,
    reason: str,
    clients_idle: bool,
    targets: str,
) -> dict:
    if not targets.strip():
        raise ValueError("All affected user/client/server targets must be named")
    before = _snapshot(expected_source, expected_commit, clean=False)
    proposed = _snapshot(source, commit, clean=True)
    with guard.change(
        owner=owner,
        reason=reason,
        clients_idle=clients_idle,
        scope={"operation": "radia-mcp-editable", "targets": targets},
        expected=before,
        proposed=proposed,
    ) as record:
        observed = installed()
        record["observed_before"] = observed
        if any(observed[k] != before[k] for k in ("source", "commit", "version")):
            raise ValueError("Installed source changed since approval; re-plan")
        # Recheck inside exclusion, immediately before the potentially mutating call.
        _snapshot(source, commit, clean=True)
        if proposed != before:
            _pip(source.resolve(strict=True))
        after = installed()
        record["observed"] = after
        if any(after[k] != proposed[k] for k in ("source", "commit", "version")):
            raise ValueError("Editable activation verification failed")
        _snapshot(source, commit, clean=True)
    return {
        "status": "completed",
        "change_id": record["change_id"],
        "receipt": str(guard.STATE_ROOT / "receipts" / (record["change_id"] + ".json")),
        "observed": after,
        "live_status": "unverified",
        "client_reconnect": "not-performed",
    }


def observe_scope(scope: dict) -> dict:
    if scope["operation"] == "radia-mcp-editable":
        return installed()
    if scope["operation"] == "config":
        from .maintenance import _digest

        path = Path(scope["path"])
        return {
            "sha256": _digest(path.read_bytes() if path.exists() else None),
            "scope": "on-disk-only",
        }
    raise ValueError("Unknown maintenance operation; manual investigation required")
