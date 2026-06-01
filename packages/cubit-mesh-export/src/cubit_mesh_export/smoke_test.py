"""
End-to-end smoke test for a deployed cubit-mesh-export plugin.

Runs Cubit in ``-batch -nographics`` mode on the canonical
``ih_bem_sample.jou`` from the radia package, exports a high-order .vol
via ``export netgen``, then asserts that the .vol contains the
expected boundary labels (``source``, ``sink``, ``sibc``). Used as the
last gate after ``cubit-plugin-install`` to prove the full round-trip
is healthy.

Usage::

    cubit-smoke-test                  # uses ih_bem_sample.jou
    cubit-smoke-test --order 3        # mesh order (default 2)
    cubit-smoke-test --keep           # keep temp workdir for inspection
    cubit-smoke-test --jou X.jou      # override sample file
    cubit-smoke-test --expect src sink sibc   # override required labels

Exit codes:
    0 -- every expected boundary label is present in the exported .vol
    1 -- any step failed (Cubit not found, .jou missing, export failed,
        labels missing)
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
import textwrap
import time
from pathlib import Path


def _find_cubit_exe():
    """Locate coreform_cubit.exe (Windows) or cubit (POSIX)."""
    env_path = os.environ.get("CUBIT_EXE")
    if env_path and Path(env_path).is_file():
        return Path(env_path)
    if sys.platform == "win32":
        import glob
        for base in [os.environ.get("ProgramFiles", "")]:
            if not base:
                continue
            for c in sorted(
                    glob.glob(os.path.join(base, "Coreform Cubit *",
                                           "bin", "coreform_cubit.exe")),
                    reverse=True):
                if Path(c).is_file():
                    return Path(c)
    else:
        for name in ("coreform_cubit", "cubit"):
            found = shutil_which_fallback(name)
            if found:
                return Path(found)
    return None


def shutil_which_fallback(name):
    """``shutil.which`` without importing shutil lazily -- portable."""
    import shutil
    return shutil.which(name)


def _find_sample_jou(override: str = "") -> Path:
    """Locate the canonical IH BEM sample .jou."""
    if override:
        p = Path(override)
        if not p.is_file():
            raise RuntimeError(f"--jou argument does not exist: {p}")
        return p

    # From installed radia package
    try:
        import radia
        rad_dir = Path(radia.__file__).resolve().parent
        candidate = rad_dir / "panels" / "samples" / "ih_bem_sample.jou"
        if candidate.is_file():
            return candidate
    except ImportError:
        pass

    raise RuntimeError(
        "Cannot locate ih_bem_sample.jou. Either install the radia "
        "package (`pip install radia`) or pass --jou <path>.")


def _read_vol_named_section(vol_path: Path, section: str):
    """Read a named-items section (e.g. ``bcnames``, ``materials``) of a
    Netgen .vol text file. Both sections share the layout::

        <section>
        <N>
        <id> "<name>"
        ...

    Returns a list of names in declaration order.
    """
    names = []
    with open(vol_path, "r", encoding="utf-8", errors="replace") as f:
        inside = False
        count_line_seen = False
        remaining = 0
        for line in f:
            s = line.strip()
            if not inside:
                if s == section:
                    inside = True
                    count_line_seen = False
                continue
            if not count_line_seen:
                try:
                    remaining = int(s.split()[0])
                except Exception:
                    remaining = 0
                count_line_seen = True
                continue
            if remaining <= 0 or not s:
                break
            parts = s.split(None, 1)
            if len(parts) == 2:
                names.append(parts[1].strip('"'))
            elif len(parts) == 1:
                names.append(parts[0].strip('"'))
            remaining -= 1
    return names


def _read_vol_bcnames(vol_path: Path):
    return _read_vol_named_section(vol_path, "bcnames")


def _read_vol_materials(vol_path: Path):
    return _read_vol_named_section(vol_path, "materials")


def run_smoke_test(*, jou: str = "", order: int = 2,
                    expect: list[str] | None = None,
                    expect_materials: list[str] | None = None,
                    keep: bool = False, timeout: float = 600.0) -> int:
    """Run the round-trip. Returns 0 on success, 1 on any failure.

    expect: boundary (sideset) labels required in the .vol bcnames section.
    expect_materials: block (volume / material) names required in the .vol
        materials section. Block names live in a different section than
        sidesets -- do not mix them up (2026-04-14 false-FAIL).
    """
    if expect is None:
        expect = ["source", "sink", "sibc"]
    if expect_materials is None:
        expect_materials = ["coil", "workpiece", "air"]

    cubit_exe = _find_cubit_exe()
    if cubit_exe is None:
        print("[FAIL] coreform_cubit.exe not found. "
              "Install Cubit or set CUBIT_EXE.")
        return 1
    print(f"  Cubit:  {cubit_exe}")

    try:
        sample = _find_sample_jou(jou)
    except RuntimeError as e:
        print(f"[FAIL] {e}")
        return 1
    print(f"  Sample: {sample}")

    work = Path(tempfile.mkdtemp(prefix="cubit-smoke-"))
    vol_path = work / "smoke.vol"
    driver = work / "driver.jou"

    # Cubit-side driver: play the sample, export .vol via the plugin,
    # then exit cleanly.  The expected labels (expect_materials +
    # expect) are validated AFTER the run by parsing the .vol file's
    # `materials` and `bcnames` sections directly -- the .vol is the
    # source of truth.
    #
    # Trailing `exit 0` matters on slower boxes (100号機 2026-04-22):
    # Cubit's headless teardown sometimes access-violates before the
    # mesh database destructor flushes the .vol writer; an explicit
    # ``exit 0`` forces the Python/Qt shutdown path through the
    # normal exit handler so the .vol is flushed first.
    driver.write_text(textwrap.dedent(f"""\
        play "{sample.as_posix()}"
        export netgen "{vol_path.as_posix()}" order {order} overwrite
        exit 0
    """), encoding="utf-8")
    print(f"  Work:   {work}")
    print(f"  Vol:    {vol_path}")
    print()

    t0 = time.time()
    cmd = [str(cubit_exe), "-batch", "-nographics",
           "-nojournal", str(driver)]
    print(f"  Running: {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=timeout, cwd=str(work),
                           encoding="utf-8", errors="replace")
    dt = time.time() - t0
    print(f"  Cubit exit={proc.returncode} ({dt:.1f}s)")

    log_path = work / "cubit.log"
    log_path.write_text(
        "=== STDOUT ===\n" + (proc.stdout or "") +
        "\n=== STDERR ===\n" + (proc.stderr or ""),
        encoding="utf-8", errors="replace")
    print(f"  Log:    {log_path}")

    # The .vol is the source of truth.  The verify_launcher command was
    # removed in cubit-mesh-export 0.8.0 along with the unified panel
    # launcher (target-centric architecture, 2026-05-05).

    # Cubit's headless mode is flaky: it often segfaults in the mesh-cleanup
    # stage AFTER export has written the .vol. We therefore trust the
    # .vol as the source of truth -- its presence + valid bcnames means the
    # plugin round-trip succeeded, regardless of Cubit's exit code.
    if not vol_path.is_file():
        print(f"[FAIL] export did not produce {vol_path}")
        if proc.returncode != 0:
            print(f"[DIAG] Cubit exited with {proc.returncode} "
                  "(0xC0000005 = access violation is common on exit).")
        print("----- last 30 lines of stderr -----")
        for line in (proc.stderr or "").splitlines()[-30:]:
            print(f"    {line}")
        return 1

    if proc.returncode != 0:
        print(f"[WARN] Cubit exited {proc.returncode} after exporting the "
              ".vol. Cubit's headless teardown can segfault after a "
              "successful export; continuing because the .vol looks valid.")
    print(f"  .vol size: {vol_path.stat().st_size} bytes")

    # Inspect boundary labels (sidesets)
    try:
        bcnames = _read_vol_bcnames(vol_path)
    except Exception as e:
        print(f"[FAIL] could not parse .vol bcnames: {e}")
        return 1
    print(f"  bcnames  ({len(bcnames)}): {bcnames}")

    # Inspect block / material labels
    try:
        materials = _read_vol_materials(vol_path)
    except Exception as e:
        print(f"[FAIL] could not parse .vol materials: {e}")
        return 1
    print(f"  materials ({len(materials)}): {materials}")

    missing_bnd = [name for name in expect if name not in bcnames]
    missing_mat = [name for name in expect_materials if name not in materials]
    if missing_bnd or missing_mat:
        if missing_bnd:
            print(f"[FAIL] missing expected boundary label(s): {missing_bnd}")
            print(f"       bcnames present: {bcnames}")
        if missing_mat:
            print(f"[FAIL] missing expected material/block name(s): {missing_mat}")
            print(f"       materials present: {materials}")
        return 1

    print()
    print(f"[OK] round-trip healthy. "
          f"Boundaries {expect} + materials {expect_materials} all present.")

    if not keep:
        import shutil as _sh
        try:
            _sh.rmtree(work)
        except Exception:
            pass
    else:
        print(f"  (--keep: workdir retained at {work})")
    return 0


def main():
    parser = argparse.ArgumentParser(
        prog="cubit-smoke-test",
        description="End-to-end smoke test: Cubit -batch -> export "
                    "netgen -> boundary-label check on the exported .vol.")
    parser.add_argument("--jou", default="",
                        help="override the source .jou (default: "
                             "radia/panels/samples/ih_bem_sample.jou)")
    parser.add_argument("--order", type=int, default=2,
                        help="mesh curving order passed to export "
                             "netgen (default 2)")
    parser.add_argument("--expect", nargs="+",
                        default=["source", "sink", "sibc"],
                        help="boundary (sideset) labels required in the "
                             ".vol bcnames section "
                             "(default: source sink sibc)")
    parser.add_argument("--expect-materials", nargs="+",
                        default=["coil", "workpiece", "air"],
                        help="block / material names required in the .vol "
                             "materials section "
                             "(default: coil workpiece air)")
    parser.add_argument("--keep", action="store_true",
                        help="keep the temp work directory for inspection")
    parser.add_argument("--timeout", type=float, default=600.0,
                        help="seconds to wait for Cubit (default 600)")
    args = parser.parse_args()

    print("=" * 60)
    print("  cubit-mesh-export: smoke test")
    print("=" * 60)
    print()
    rc = run_smoke_test(jou=args.jou, order=args.order,
                         expect=args.expect,
                         expect_materials=args.expect_materials,
                         keep=args.keep,
                         timeout=args.timeout)
    raise SystemExit(rc)


if __name__ == "__main__":
    main()
