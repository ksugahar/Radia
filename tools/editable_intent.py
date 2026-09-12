"""Recorded editable-source intent for the LAB / 100号機 development tier.

The record answers one question per (host, interpreter, package): which
source tree is *meant* to be installed editable right now. Verification
compares the installed registration and a fresh-process import against that
record. A package with no record is reported as UNVERIFIED, not as drift, and
no repair toward any historical or "canonical" path is ever suggested: the only
way to change the intent is an explicit ``repoint`` with a reason.

Layout of the record file (``RADIA_EDITABLE_INTENT_FILE`` or
``%ProgramData%\\Radia\\editable-intent.json``)::

    {"schema": "radia.editable-intent.v1", "host": "LAB",
     "interpreters": {"c:/program files/python312/python.exe": {"packages": {
         "radia-mcp": {"source": "S:/Radia/release-quad/x/packages/radia-mcp",
                       "commit": "<40 hex>", "tracked_clean": true,
                       "recorded_at": "...Z", "recorded_by": "user@host",
                       "recorded_via": "repoint", "reason": "...",
                       "pushed_refs": [...], "previous": {...} | null}}}}}

Every change appends one JSON line to ``editable-intent.log.jsonl`` next to
the record, so the previous pointer is never lost.

``repoint`` performs ``pip install -e <source>`` for the named packages only.
It records the previous pointer before pip runs, verifies registration and
import afterwards, and leaves the record untouched when either step fails so
``repoint --rollback`` can reinstall the previous source. It never stops
processes and never uninstalls: a locked entry point is reported as a blocked
update for the operator to resolve at a quiet boundary.

The file is dependency-free so it can be piped to a remote interpreter as
``python - --argv-b64 <base64 json argv>``.
"""
from __future__ import annotations

import argparse
import base64
import getpass
import json
import os
import platform
import shutil
import subprocess
import sys
import sysconfig
import tempfile
from datetime import UTC, datetime
from importlib import metadata as importlib_metadata
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import url2pathname

SCHEMA = "radia.editable-intent.v1"
INTENT_FILE_ENV = "RADIA_EDITABLE_INTENT_FILE"
DEFAULT_PACKAGES = ("radia", "cubit-mesh-export", "radia-mcp")
IMPORT_MODULES = {
    "radia": "radia",
    "cubit-mesh-export": "cubit_mesh_export",
    "radia-mcp": "radia_mcp",
    "mcp-server-document": "mcp_server_document",
}

# Both UNC spellings of the lab NAS, LAB's S: drive and 100号機's W: drive
# all name the same Radia tree.
_PATH_ALIASES = (
    ("//192.168.11.100/work/00_cae/radia/", "s:/radia/"),
    ("//192.168.121.100/work/00_cae/radia/", "s:/radia/"),
    ("w:/00_cae/radia/", "s:/radia/"),
)

EXIT_OK = 0
EXIT_DRIFT = 1
EXIT_PRECONDITION = 2
EXIT_ACTION = 3
EXIT_VERIFY = 4
EXIT_UNVERIFIED = 5


# ------------------------------------------------------------------
# record store
# ------------------------------------------------------------------

def intent_file_path() -> Path:
    explicit = os.environ.get(INTENT_FILE_ENV, "").strip()
    if explicit:
        return Path(explicit)
    if os.name == "nt":
        base = Path(os.environ.get("ProgramData", r"C:\ProgramData"))
        return base / "Radia" / "editable-intent.json"
    return Path.home() / ".radia" / "editable-intent.json"


def log_file_path(path: Path | None = None) -> Path:
    path = path or intent_file_path()
    return path.with_name(path.stem + ".log.jsonl")


def interpreter_key(executable: str | None = None) -> str:
    return norm_path(executable or sys.executable)


def load_intent(path: Path | None = None) -> dict:
    path = path or intent_file_path()
    if not path.exists():
        return {"schema": SCHEMA, "host": platform.node(), "interpreters": {}}
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema") != SCHEMA:
        raise ValueError(f"{path}: unexpected schema {data.get('schema')!r}")
    data.setdefault("interpreters", {})
    return data


def save_intent(data: dict, path: Path | None = None) -> Path:
    path = path or intent_file_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    data["schema"] = SCHEMA
    data.setdefault("host", platform.node())
    text = json.dumps(data, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    handle, tmp_name = tempfile.mkstemp(dir=str(path.parent), prefix=".intent-", suffix=".tmp")
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as tmp:
            tmp.write(text)
        os.replace(tmp_name, path)
    finally:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
    return path


def append_log(event: dict, path: Path | None = None) -> Path:
    log_path = log_file_path(path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    record = {"at": _now(), "host": platform.node(), "by": _actor(), **event}
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True, ensure_ascii=False) + "\n")
    return log_path


def recorded_entry(data: dict, package: str, executable: str | None = None) -> dict | None:
    interp = data.get("interpreters", {}).get(interpreter_key(executable), {})
    return interp.get("packages", {}).get(package)


def recorded_packages(data: dict, executable: str | None = None) -> list[str]:
    interp = data.get("interpreters", {}).get(interpreter_key(executable), {})
    return sorted(interp.get("packages", {}))


def set_entry(data: dict, package: str, entry: dict, executable: str | None = None) -> None:
    interp = data.setdefault("interpreters", {}).setdefault(interpreter_key(executable), {})
    interp.setdefault("executable", executable or sys.executable)
    interp.setdefault("packages", {})[package] = entry


# ------------------------------------------------------------------
# observation
# ------------------------------------------------------------------

def norm_path(p: str | os.PathLike | None) -> str:
    text = str(p or "").replace("\\", "/").rstrip("/").lower()
    for alias, canonical in _PATH_ALIASES:
        text = text.replace(alias, canonical)
    return text


def _editable_location(dist) -> str | None:
    raw = dist.read_text("direct_url.json")
    if not raw:
        return None
    try:
        data = json.loads(raw)
    except json.JSONDecodeError:
        return None
    if not data.get("dir_info", {}).get("editable") or not data.get("url"):
        return None
    parsed = urlparse(data["url"])
    if parsed.scheme != "file":
        return data["url"]
    path = url2pathname(parsed.path)
    if parsed.netloc:
        return f"//{parsed.netloc}{path}"
    if os.name == "nt" and path.startswith("/") and len(path) > 2 and path[2] == ":":
        path = path[1:]
    return path


def _pth_files(package: str) -> list[dict]:
    """List ``__editable__*.pth`` files belonging to the distribution."""
    stem = package.lower().replace("-", "_")
    found = []
    purelib = sysconfig.get_paths().get("purelib")
    if not purelib:
        return found
    try:
        entries = sorted(Path(purelib).glob("__editable__*.pth"))
    except OSError:
        return found
    for entry in entries:
        name = entry.name.lower().replace("-", "_")
        if not name.startswith(f"__editable__.{stem}"):
            continue
        try:
            content = entry.read_text(encoding="utf-8", errors="replace").strip()
        except OSError:
            content = "<unreadable>"
        found.append({"file": str(entry), "content": content})
    return found


def observe_registration(package: str) -> dict:
    """What this interpreter has registered for ``package`` (no import)."""
    try:
        dist = importlib_metadata.distribution(package)
    except importlib_metadata.PackageNotFoundError:
        return {"installed": False, "version": None, "editable": False,
                "source": None, "pth": []}
    location = _editable_location(dist)
    return {
        "installed": True,
        "version": dist.version,
        "editable": location is not None,
        "source": location,
        "pth": _pth_files(package),
    }


def fresh_import_origin(package: str, executable: str | None = None) -> str | None:
    """``__file__`` of the module a *new* process imports, or None if unknown."""
    module = IMPORT_MODULES.get(package)
    if not module:
        return None
    probe = (
        "import importlib; "
        f"m = importlib.import_module({module!r}); "
        "print(m.__file__ or '')"
    )
    try:
        result = subprocess.run(
            [executable or sys.executable, "-s", "-c", probe],
            capture_output=True, text=True, timeout=300, check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if result.returncode != 0:
        return None
    return (result.stdout or "").strip() or None


def _repo_root(path: Path) -> Path | None:
    """Nearest ancestor holding ``.git`` (a directory, or a worktree's file)."""
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    return None


def _git(root: Path, *args: str) -> subprocess.CompletedProcess | None:
    """Run git in ``root``; safe.directory must name the repository top level."""
    git = shutil.which("git")
    if git is None:
        return None
    try:
        return subprocess.run(
            [git, "-c", f"safe.directory={root.as_posix()}", "-C", str(root), *args],
            capture_output=True, text=True, timeout=60, check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None


def describe_source(source: str | os.PathLike) -> dict:
    """Existence, project marker, commit and tracked-clean state of a source tree."""
    path = Path(source)
    info = {
        "source": str(source),
        "exists": path.is_dir(),
        "project": (path / "pyproject.toml").is_file() or (path / "setup.py").is_file(),
        "repository": None,
        "commit": None,
        "tracked_clean": None,
    }
    if not info["exists"]:
        return info
    root = _repo_root(path)
    if root is None:
        return info
    info["repository"] = str(root)
    head = _git(root, "rev-parse", "HEAD")
    if head is not None and head.returncode == 0:
        info["commit"] = head.stdout.strip().lower() or None
        status = _git(root, "status", "--porcelain", "--untracked-files=no")
        if status is not None and status.returncode == 0:
            info["tracked_clean"] = not status.stdout.strip()
    return info


def pushed_refs(source: str | os.PathLike, commit: str | None) -> list[str] | None:
    """Remote branches and tags that contain ``commit``; None when unknown."""
    if not commit:
        return None
    root = _repo_root(Path(source))
    if root is None:
        return None
    result = _git(root, "for-each-ref", "--format=%(refname:short)",
                  f"--contains={commit}", "refs/remotes", "refs/tags")
    if result is None or result.returncode != 0:
        return None
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def evaluate(package: str, expected: str | None, executable: str | None = None) -> dict:
    """One verification row: registration + fresh import against ``expected``."""
    observed = observe_registration(package)
    origin = fresh_import_origin(package, executable) if observed["installed"] else None
    row = {
        "package": package,
        "expected": expected,
        "version": observed["version"],
        "registration": observed["source"],
        "import": origin,
        "pth": observed["pth"],
    }
    if not observed["installed"]:
        row["status"] = "not_installed"
    elif not observed["editable"]:
        row["status"] = "not_editable"
    elif expected is None:
        row["status"] = "unverified"
    else:
        want = norm_path(expected)
        registration_ok = norm_path(observed["source"]) == want
        import_known = origin is not None
        import_ok = import_known and norm_path(origin).startswith(want + "/")
        if registration_ok and (import_ok or not import_known):
            row["status"] = "match" if import_known else "match_registration_only"
        else:
            row["status"] = "drift"
    return row


def verify(packages: list[str], expectations: dict[str, str | None] | None = None,
           data: dict | None = None, executable: str | None = None) -> dict:
    """Verify ``packages``; expectation = explicit override, else the record."""
    data = data if data is not None else load_intent()
    expectations = expectations or {}
    rows = []
    for package in packages:
        entry = recorded_entry(data, package, executable)
        if package in expectations:
            expected, basis = expectations[package], "explicit"
        elif entry:
            expected, basis = entry["source"], "record"
        else:
            expected, basis = None, "none"
        row = evaluate(package, expected, executable)
        row["expectation_basis"] = basis
        row["recorded"] = entry
        rows.append(row)
    counts = {"match": 0, "drift": 0, "unverified": 0}
    for row in rows:
        if row["status"] in ("match", "match_registration_only"):
            counts["match"] += 1
        elif row["status"] == "unverified":
            counts["unverified"] += 1
        else:
            counts["drift"] += 1
    return {"schema": SCHEMA, "host": platform.node(), "interpreter": executable or sys.executable,
            "record_file": str(intent_file_path()), "rows": rows, "counts": counts,
            "exit_code": _verify_exit(counts)}


def _verify_exit(counts: dict) -> int:
    if counts["drift"]:
        return EXIT_DRIFT
    if counts["unverified"]:
        return EXIT_UNVERIFIED
    return EXIT_OK


def format_verify(report: dict) -> str:
    lines = [f"editable intent record: {report['record_file']}",
             f"interpreter: {report['interpreter']}"]
    for row in report["rows"]:
        status = row["status"].upper()
        lines.append(f"  [{status:<24}] {row['package']} v{row['version'] or '?'}")
        lines.append(f"      expected ({row['expectation_basis']}): {row['expected'] or '<no record>'}")
        lines.append(f"      registration:        {row['registration'] or '<none>'}")
        lines.append(f"      fresh-process import: {row['import'] or '<not probed>'}")
        recorded = row.get("recorded")
        if recorded:
            lines.append(f"      recorded {recorded.get('recorded_at')} via {recorded.get('recorded_via')}"
                         f" by {recorded.get('recorded_by')}: {recorded.get('reason') or ''}")
        if row["status"] == "unverified":
            lines.append("      no recorded intent: not counted as drift. If this pointer is the "
                         "intended one, run `repoint --record-current --reason ...`; "
                         "otherwise repoint explicitly. Nothing is repaired automatically.")
        elif row["status"] == "drift":
            lines.append("      differs from the recorded intent. Either the record or the "
                         "installation is wrong: decide which, then `repoint --package "
                         f"{row['package']} --source <intended> --reason ...` or re-record. "
                         "No default target exists.")
    counts = report["counts"]
    lines.append(f"  match={counts['match']} drift={counts['drift']} unverified={counts['unverified']}")
    lines.append("  Running server processes keep the code they loaded; verify each live client "
                 "separately (status tools report runtime_provenance).")
    return "\n".join(lines)


# ------------------------------------------------------------------
# repoint / record / rollback
# ------------------------------------------------------------------

def _now() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _actor() -> str:
    try:
        user = getpass.getuser()
    except (KeyError, OSError, ImportError):  # pragma: no cover - platform specific
        user = "unknown"
    return f"{user}@{platform.node()}"


def _entry_from_source(source: str, reason: str, via: str, previous: dict | None,
                       require_pushed: bool = False) -> tuple[dict | None, str | None]:
    info = describe_source(source)
    if not info["exists"]:
        return None, f"source does not exist: {source}"
    if not info["project"]:
        return None, f"source has no pyproject.toml/setup.py: {source}"
    refs = pushed_refs(source, info["commit"])
    if require_pushed and not refs:
        return None, (f"commit {info['commit'] or '<unknown>'} of {source} is not contained in "
                      "any remote branch or tag; --require-pushed refuses it")
    entry = {
        "source": str(source),
        "commit": info["commit"],
        "tracked_clean": info["tracked_clean"],
        "recorded_at": _now(),
        "recorded_by": _actor(),
        "recorded_via": via,
        "reason": reason,
        "pushed_refs": refs,
        "previous": previous,
    }
    return entry, None


def _previous_of(data: dict, package: str, executable: str | None = None) -> dict | None:
    """Snapshot of what is installed now, to keep under ``previous``."""
    observed = observe_registration(package)
    recorded = recorded_entry(data, package, executable)
    if not observed["installed"] or not observed["editable"]:
        snapshot = {"source": None, "installed": observed["installed"],
                    "editable": observed["editable"], "version": observed["version"]}
    else:
        snapshot = {"source": observed["source"], "version": observed["version"]}
        info = describe_source(observed["source"])
        snapshot["commit"] = info["commit"]
    if recorded:
        snapshot["recorded"] = {k: v for k, v in recorded.items() if k != "previous"}
    snapshot["observed_at"] = _now()
    return snapshot


def pip_install_editable(source: str, executable: str | None = None,
                         runner=subprocess.run) -> subprocess.CompletedProcess:
    command = [executable or sys.executable, "-m", "pip", "install", "--no-deps",
               "--no-cache-dir", "--no-build-isolation", "-e", str(source)]
    return runner(command, capture_output=True, text=True)


def repoint(targets: list[tuple[str, str]], reason: str, *, via: str = "repoint",
            dry_run: bool = False, require_pushed: bool = False,
            executable: str | None = None, runner=subprocess.run,
            data: dict | None = None, path: Path | None = None) -> dict:
    """Point each (package, source) at its new editable source, recording both ends."""
    path = path or intent_file_path()
    data = data if data is not None else load_intent(path)
    results = []
    exit_code = EXIT_OK
    for package, source in targets:
        previous = _previous_of(data, package, executable)
        entry, error = _entry_from_source(source, reason, via, previous, require_pushed)
        if error:
            results.append({"package": package, "source": source, "status": "refused", "detail": error})
            exit_code = EXIT_PRECONDITION
            break
        if dry_run:
            results.append({"package": package, "source": source, "status": "planned",
                            "previous": previous, "entry": entry,
                            "command": [executable or sys.executable, "-m", "pip", "install",
                                        "--no-deps", "--no-cache-dir", "--no-build-isolation",
                                        "-e", str(source)]})
            continue
        append_log({"event": "repoint-pre", "package": package, "target": source,
                    "previous": previous, "reason": reason, "via": via}, path)
        proc = pip_install_editable(source, executable, runner)
        tail = "\n".join((proc.stdout or "").splitlines()[-8:] + (proc.stderr or "").splitlines()[-8:])
        if proc.returncode != 0:
            append_log({"event": "repoint-blocked", "package": package, "target": source,
                        "returncode": proc.returncode, "tail": tail}, path)
            results.append({"package": package, "source": source, "status": "blocked",
                            "returncode": proc.returncode, "detail": tail,
                            "previous": previous,
                            "hint": ("pip did not complete; the record was not changed. A locked "
                                     "entry point means an affected server is still running: stop "
                                     "only that client at a quiet boundary and retry, or run "
                                     "`repoint --rollback --package " + package + "`.")})
            exit_code = EXIT_ACTION
            break
        row = evaluate(package, source, executable)
        if row["status"] not in ("match", "match_registration_only"):
            append_log({"event": "repoint-verify-failed", "package": package, "target": source,
                        "row": {k: row[k] for k in ("status", "registration", "import")}}, path)
            results.append({"package": package, "source": source, "status": "verify_failed",
                            "row": row, "previous": previous,
                            "hint": ("pip reported success but the interpreter does not resolve "
                                     "the new source; the record was not changed. Inspect the "
                                     ".pth files listed above, then retry or "
                                     "`repoint --rollback --package " + package + "`.")})
            exit_code = EXIT_VERIFY
            break
        set_entry(data, package, entry, executable)
        save_intent(data, path)
        append_log({"event": "repoint-post", "package": package, "target": source,
                    "commit": entry["commit"], "tracked_clean": entry["tracked_clean"],
                    "reason": reason, "via": via}, path)
        results.append({"package": package, "source": source, "status": "repointed",
                        "row": row, "entry": entry, "previous": previous})
    return {"schema": SCHEMA, "host": platform.node(), "interpreter": executable or sys.executable,
            "record_file": str(path), "dry_run": dry_run, "results": results, "exit_code": exit_code,
            "note": ("Existing server processes keep the previous source until they are "
                     "reconnected; verify each live client separately.")}


def record_current(packages: list[str], reason: str, *, via: str = "record-current",
                   require_pushed: bool = False, executable: str | None = None,
                   data: dict | None = None, path: Path | None = None) -> dict:
    """Adopt the currently installed editable pointers as the recorded intent."""
    path = path or intent_file_path()
    data = data if data is not None else load_intent(path)
    results = []
    exit_code = EXIT_OK
    for package in packages:
        observed = observe_registration(package)
        if not observed["installed"] or not observed["editable"]:
            results.append({"package": package, "status": "refused",
                            "detail": "not installed editable; nothing to adopt as intent"})
            exit_code = EXIT_PRECONDITION
            continue
        previous = _previous_of(data, package, executable)
        entry, error = _entry_from_source(observed["source"], reason, via, previous, require_pushed)
        if error:
            results.append({"package": package, "status": "refused", "detail": error})
            exit_code = EXIT_PRECONDITION
            continue
        set_entry(data, package, entry, executable)
        save_intent(data, path)
        append_log({"event": "record-current", "package": package, "target": entry["source"],
                    "commit": entry["commit"], "reason": reason, "via": via}, path)
        results.append({"package": package, "source": entry["source"], "status": "recorded",
                        "entry": entry})
    return {"schema": SCHEMA, "host": platform.node(), "interpreter": executable or sys.executable,
            "record_file": str(path), "results": results, "exit_code": exit_code}


def rollback(packages: list[str], reason: str, *, executable: str | None = None,
             runner=subprocess.run, data: dict | None = None, path: Path | None = None) -> dict:
    """Reinstall each package from the ``previous`` pointer kept in its record."""
    path = path or intent_file_path()
    data = data if data is not None else load_intent(path)
    targets = []
    refused = []
    for package in packages:
        entry = recorded_entry(data, package, executable)
        previous = (entry or {}).get("previous") or {}
        source = previous.get("source")
        if not source:
            refused.append({"package": package, "status": "refused",
                            "detail": "no previous editable pointer recorded for this package"})
            continue
        targets.append((package, source))
    if refused:
        return {"schema": SCHEMA, "record_file": str(path), "results": refused,
                "exit_code": EXIT_PRECONDITION}
    return repoint(targets, reason, via="rollback", executable=executable, runner=runner,
                   data=data, path=path)


def format_action(report: dict) -> str:
    lines = [f"editable intent record: {report['record_file']}"]
    if report.get("dry_run"):
        lines.append("  dry run: nothing was installed or recorded")
    for result in report["results"]:
        lines.append(f"  [{result['status'].upper():<14}] {result['package']} -> {result.get('source') or ''}")
        if result.get("detail"):
            lines.append("      " + str(result["detail"]).replace("\n", "\n      "))
        if result.get("hint"):
            lines.append("      " + result["hint"])
        entry = result.get("entry")
        if entry:
            lines.append(f"      commit {entry['commit'] or '<not a git tree>'}"
                         f" tracked_clean={entry['tracked_clean']} pushed_refs={entry['pushed_refs']}")
        previous = result.get("previous")
        if previous:
            lines.append(f"      previous pointer: {previous.get('source') or '<not editable>'}")
    if report.get("note"):
        lines.append("  " + report["note"])
    return "\n".join(lines)


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------

def _parse_expectations(values: list[str]) -> dict[str, str]:
    expectations = {}
    for value in values or []:
        if "=" not in value:
            raise SystemExit(f"--expect needs PACKAGE=PATH, got {value!r}")
        package, path = value.split("=", 1)
        expectations[package.strip()] = path.strip()
    return expectations


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="editable_intent", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json", action="store_true", help="print the JSON report only")
    parser.add_argument("--argv-b64", help=argparse.SUPPRESS)
    sub = parser.add_subparsers(dest="cmd", required=True)

    show = sub.add_parser("show", help="print the recorded intent for this interpreter")
    show.set_defaults(func=_cmd_show)

    ver = sub.add_parser("verify", help="compare installed pointers with the recorded intent")
    ver.add_argument("--package", action="append", default=[],
                     help="package to verify (default: the three repo packages plus any recorded)")
    ver.add_argument("--expect", action="append", default=[],
                     help="PACKAGE=PATH explicit expectation overriding the record")
    ver.set_defaults(func=_cmd_verify)

    rep = sub.add_parser("repoint", help="install a package editable from a named source and record it")
    rep.add_argument("--package", action="append", default=[], help="package name (repeatable)")
    rep.add_argument("--source", action="append", default=[],
                     help="source tree for the package in the same position (repeatable)")
    rep.add_argument("--reason", default="", help="why the pointer moves (recorded)")
    rep.add_argument("--via", default="repoint", help="recorded_via label")
    rep.add_argument("--record-current", action="store_true",
                     help="adopt the currently installed pointers as intent; no pip")
    rep.add_argument("--rollback", action="store_true",
                     help="reinstall from the previous pointer kept in the record")
    rep.add_argument("--dry-run", action="store_true", help="plan only")
    rep.add_argument("--require-pushed", action="store_true",
                     help="refuse a source commit that no remote branch or tag contains "
                          "(formal handoff / completion evidence)")
    rep.set_defaults(func=_cmd_repoint)
    return parser


def _emit(args, report: dict, text: str) -> int:
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False))
    else:
        print(text)
    return int(report.get("exit_code", EXIT_OK))


def _cmd_show(args) -> int:
    data = load_intent()
    report = {"schema": SCHEMA, "record_file": str(intent_file_path()),
              "interpreter": sys.executable,
              "packages": {p: recorded_entry(data, p) for p in recorded_packages(data)},
              "exit_code": EXIT_OK}
    text = json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False)
    return _emit(args, report, text)


def _cmd_verify(args) -> int:
    data = load_intent()
    packages = list(args.package) or sorted(set(DEFAULT_PACKAGES) | set(recorded_packages(data)),
                                            key=lambda p: (p not in DEFAULT_PACKAGES, p))
    report = verify(packages, _parse_expectations(args.expect), data=data)
    return _emit(args, report, format_verify(report))


def _cmd_repoint(args) -> int:
    if args.record_current and args.rollback:
        raise SystemExit("--record-current and --rollback are exclusive")
    if args.rollback:
        packages = list(args.package) or list(DEFAULT_PACKAGES)
        report = rollback(packages, args.reason or "rollback to the previous recorded pointer")
        return _emit(args, report, format_action(report))
    if args.record_current:
        if args.source:
            raise SystemExit("--record-current adopts installed pointers; do not pass --source")
        if not args.reason:
            raise SystemExit("--reason is required so the record says why this pointer is intended")
        packages = list(args.package) or list(DEFAULT_PACKAGES)
        report = record_current(packages, args.reason, require_pushed=args.require_pushed)
        return _emit(args, report, format_action(report))
    if not args.package or len(args.package) != len(args.source):
        raise SystemExit("repoint needs matching --package/--source pairs")
    if not args.reason:
        raise SystemExit("--reason is required so the record says why this pointer moves")
    report = repoint(list(zip(args.package, args.source)), args.reason, via=args.via,
                     dry_run=args.dry_run, require_pushed=args.require_pushed)
    return _emit(args, report, format_action(report))


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) == 2 and argv[0] == "--argv-b64":
        argv = json.loads(base64.b64decode(argv[1]).decode("utf-8"))
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
