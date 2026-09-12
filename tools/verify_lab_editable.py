"""Read-only development editable-source verification.

Defaults to the canonical LAB checkout. Use --source-root explicitly for an
approved clean monorepo worktree, or --mcp-source for a separately approved
MCP runtime. Never infer the expectation from the installed package: that
would accept drift as its own source of truth.

Pointing at the right path is not the same as running current code. The LAB
checkout can sit on a backup branch, so an install can satisfy the path check
and still import something older than origin/main. Measured 2026-09-10:
radia-mcp resolved to the canonical tree carrying 1.4.39 while origin/main
carried 1.4.53, and the path check called that clean. So also compare the
version that actually imports against origin/main. origin/main is the
reference; the installed package still never is.
"""
from __future__ import annotations

import argparse
import importlib.util
import pathlib
import re
import subprocess

import release_quad

# distribution name -> (import name, repo-relative file holding __version__).
# mcp-server-document lives in another repository and has no origin/main here.
VERSIONED_PACKAGES = {
    "radia": ("radia", "src/radia/__init__.py"),
    "cubit-mesh-export": (
        "cubit_mesh_export",
        "packages/cubit-mesh-export/src/cubit_mesh_export/__init__.py"),
    "radia-mcp": (
        "radia_mcp", "packages/radia-mcp/src/radia_mcp/__init__.py"),
}

_VERSION = re.compile(r"""^__version__\s*=\s*["']([^"']+)["']""", re.M)
_REPO = pathlib.Path(__file__).resolve().parents[1]


def _declared_version(text):
    found = _VERSION.search(text or "")
    return found.group(1) if found else None


def _running_version(import_name):
    """Version of the file that would actually import, read without importing.

    Importing radia loads its native extension, which is far too much work for
    a check that only needs one assignment out of __init__.py.
    """
    try:
        spec = importlib.util.find_spec(import_name)
    except (ImportError, ValueError):
        return None, None
    if spec is None or not spec.origin:
        return None, None
    origin = pathlib.Path(spec.origin)
    try:
        return _declared_version(
            origin.read_text(encoding="utf-8", errors="replace")), origin
    except OSError:
        return None, origin


def _origin_main_version(relative_path):
    # Decode as UTF-8 explicitly: these files carry Japanese comments, and the
    # Windows default (cp932) raises UnicodeDecodeError on them, which would
    # silently turn a real version mismatch into a skipped comparison.
    try:
        shown = subprocess.run(
            ["git", "-c", f"safe.directory={_REPO.as_posix()}",
             "-C", str(_REPO), "show", f"origin/main:{relative_path}"],
            capture_output=True, text=True, encoding="utf-8", errors="replace")
    except OSError:
        return None
    if shown.returncode != 0:
        return None
    return _declared_version(shown.stdout)


def _release_order(version):
    return tuple(int(part) for part in re.findall(r"\d+", version))


def verify_against_origin_main(packages):
    """Count stale packages and packages whose comparison cannot be verified.

    A package that is ahead of origin/main is normal mid-development and is
    reported but not counted. Unknown versions fail the check: a missing answer
    must not read as a passing one. Explicit --skip-origin-check remains the
    opt-out for an intentionally path-only, offline check.
    """
    stale = 0
    for name, _expected_path in packages:
        if name not in VERSIONED_PACKAGES:
            continue
        import_name, relative_path = VERSIONED_PACKAGES[name]
        running, origin = _running_version(import_name)
        reference = _origin_main_version(relative_path)
        if running is None or reference is None:
            stale += 1
            release_quad.fail(
                f"{name:<26} version comparison UNVERIFIED "
                f"(running={running or '?'}, origin/main={reference or '?'})\n"
                "        Check the import source and local origin/main ref. "
                "No installation was changed.")
            continue
        if running == reference:
            release_quad.ok(f"{name:<26} {running} matches origin/main")
        elif _release_order(running) < _release_order(reference):
            stale += 1
            release_quad.fail(
                f"{name:<26} runs {running}, origin/main carries {reference}\n"
                f"        importing {origin}\n"
                f"        The path check passes, so this is not a stray\n"
                f"        pointer: the tree it points at is behind. Fetch and\n"
                f"        align that checkout, or aim the install at a tree\n"
                f"        that carries origin/main.")
        else:
            release_quad.ok(
                f"{name:<26} {running} is ahead of origin/main ({reference})")
    return stale


def expected_packages(mcp_source=None, source_root=None):
    packages = release_quad._canonical_lab_editable_packages()
    if source_root:
        root = pathlib.Path(source_root)
        monorepo_paths = {
            "radia": str(root),
            "cubit-mesh-export": str(root / "packages" / "cubit-mesh-export"),
            "radia-mcp": str(root / "packages" / "radia-mcp"),
        }
        packages = [(name, monorepo_paths.get(name, path))
                    for name, path in packages]
    if mcp_source:
        packages = [(name, mcp_source if name == "radia-mcp" else path)
                    for name, path in packages]
    return packages


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-root",
        help="Approved clean monorepo worktree used by release-quad")
    parser.add_argument("--mcp-source", help="Approved editable MCP project root")
    parser.add_argument(
        "--skip-origin-check", action="store_true",
        help="skip the origin/main comparison (offline, or origin not fetched)")
    args = parser.parse_args(argv)
    packages = expected_packages(
        mcp_source=args.mcp_source, source_root=args.source_root)
    drift = release_quad._verify_lab_editable(packages)
    if not args.skip_origin_check:
        drift += verify_against_origin_main(packages)
    return 4 if drift else 0


if __name__ == "__main__":
    raise SystemExit(main())
