#!/usr/bin/env python
"""Fail a Ninja/MSVC build whose header dependencies Ninja cannot record.

cl.exe writes its localized ``/showIncludes`` notes in the console output
code page, while CMake decodes its prefix probe with the console input code
page.  When the two differ, ``rules.ninja`` receives a garbled
``msvc_deps_prefix``; Ninja then records no header dependencies and a
header-only change rebuilds nothing while the build still succeeds.

``prefix`` checks, before compiling, that the compiler's current
``/showIncludes`` line starts with the configured prefix bytes, which is
the comparison Ninja makes.  ``deps`` checks, after a build, that every
object of a target carries a non-empty header dependency record that Ninja
marks VALID for the current object.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

_PREFIX = re.compile(rb"^[ \t]*msvc_deps_prefix[ \t]*=[ \t]*(.*?)\r?$", re.MULTILINE)
_RECORD = re.compile(r"^(?P<path>\S.*): #deps (?P<count>\d+), deps mtime \d+ \((?P<state>\w+)\)$")
_REMEDY = ("Run Build.ps1 -Rebuild; Build.ps1 aligns the console input code page with the output "
           "code page so CMake decodes cl.exe output as cl.exe writes it.")


class DepsError(RuntimeError):
    """Ninja cannot or did not record header dependencies."""


def _unescape(value: bytes) -> bytes:
    return re.sub(rb"\$([$ :])", rb"\1", value)


def configured_prefix(build_dir: Path) -> bytes:
    """Return the single ``msvc_deps_prefix`` value of a CMake Ninja tree."""
    rules = build_dir / "CMakeFiles" / "rules.ninja"
    if not rules.is_file():
        raise DepsError(f"{rules} is missing; configure the Ninja build tree first")
    values = {_unescape(v) for v in _PREFIX.findall(rules.read_bytes())}
    if len(values) != 1 or not next(iter(values)).strip():
        raise DepsError(f"{rules} must define one non-empty msvc_deps_prefix; found {sorted(values)!r}. "
                        + _REMEDY)
    return values.pop()


def compiler_include_line(compiler: str) -> bytes:
    """Compile a one-header probe and return cl.exe's raw ``/showIncludes`` line."""
    with tempfile.TemporaryDirectory(prefix="ninja-msvc-deps-") as scratch:
        root = Path(scratch)
        (root / "probe_header.h").write_bytes(b"\n")
        (root / "probe.c").write_bytes(b'#include "probe_header.h"\nint main(void){return 0;}\n')
        result = subprocess.run([compiler, "/nologo", "/showIncludes", "/c", "probe.c"],
                                cwd=root, capture_output=True, timeout=120)
    lines = [line.rstrip(b"\r") for line in result.stdout.splitlines() if b"probe_header.h" in line]
    if result.returncode != 0 or len(lines) != 1:
        raise DepsError(f"{compiler} /showIncludes probe failed (exit {result.returncode}): "
                        f"{result.stdout[-400:]!r} {result.stderr[-400:]!r}")
    return lines[0]


def check_prefix(prefix: bytes, include_line: bytes) -> None:
    """Fail unless Ninja's prefix match accepts the compiler's include line."""
    if not include_line.startswith(prefix):
        raise DepsError(
            "msvc_deps_prefix does not match the compiler's /showIncludes output, so Ninja would "
            "record no header dependencies.\n"
            f"  configured: {prefix.hex(' ')}  {prefix.decode('utf-8', 'replace')!a}\n"
            f"  compiler:   {include_line[:48].hex(' ')}  "
            f"{include_line.decode('utf-8', 'replace')!a}\n" + _REMEDY)


def recorded_deps(ninja_deps_output: str) -> dict[str, tuple[int, str]]:
    """Parse ``ninja -t deps`` into {output path: (dependency count, record state)}."""
    records = {}
    for line in ninja_deps_output.splitlines():
        match = _RECORD.match(line)
        if match:
            records[match["path"]] = (int(match["count"]), match["state"])
    return records


def check_target_deps(build_dir: Path, target: str, records: dict[str, tuple[int, str]]) -> int:
    """Fail unless every object of ``target`` has a non-empty VALID dependency record."""
    target_dir = build_dir / "CMakeFiles" / f"{target}.dir"
    objects = sorted(p.relative_to(build_dir).as_posix() for p in target_dir.rglob("*.obj"))
    if not objects:
        raise DepsError(f"{target_dir} holds no compiled objects")
    empty = [obj for obj in objects if records.get(obj, (0, ""))[0] == 0]
    if empty:
        raise DepsError(f"Ninja recorded no header dependencies for {len(empty)}/{len(objects)} "
                        f"{target} objects (first: {empty[:3]}); a header-only change would not "
                        "rebuild them. " + _REMEDY)
    invalid = [(obj, records[obj][1]) for obj in objects if records[obj][1] != "VALID"]
    if invalid:
        raise DepsError(f"Ninja's header dependency records are not VALID for {len(invalid)}/{len(objects)} "
                        f"{target} objects (first: {invalid[:3]}); they do not describe the current "
                        "objects. Run Build.ps1 -Rebuild so Ninja records them again.")
    return len(objects)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="phase", required=True)
    prefix = sub.add_parser("prefix", help="after configure: compare the prefix with cl.exe output")
    prefix.add_argument("--build-dir", type=Path, required=True)
    prefix.add_argument("--compiler", required=True)
    deps = sub.add_parser("deps", help="after a build: require recorded header dependencies")
    deps.add_argument("--build-dir", type=Path, required=True)
    deps.add_argument("--ninja", required=True)
    deps.add_argument("--target", required=True)
    args = parser.parse_args(argv)
    try:
        if args.phase == "prefix":
            check_prefix(configured_prefix(args.build_dir), compiler_include_line(args.compiler))
            print(f"Ninja msvc_deps_prefix matches {args.compiler} /showIncludes output")
        else:
            result = subprocess.run([args.ninja, "-C", str(args.build_dir), "-t", "deps"],
                                    capture_output=True, timeout=300)
            if result.returncode != 0:
                raise DepsError(f"ninja -t deps failed (exit {result.returncode}): {result.stderr[-400:]!r}")
            count = check_target_deps(args.build_dir, args.target,
                                      recorded_deps(result.stdout.decode("utf-8", "replace")))
            print(f"Ninja recorded VALID header dependencies for all {count} {args.target} objects")
    except DepsError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
