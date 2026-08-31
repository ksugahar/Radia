"""
End-to-end smoke test for a deployed cubit-mesh-export plugin.

Runs Cubit in ``-batch -nographics`` mode on the canonical
``ih_bem_sample.jou`` from the radia package and exports a high-order .vol
via ``export netgen``. The result must then pass the same production
``check-vol`` gate used before solver initialization: strict labels, complete
DomainIn/DomainOut ownership, NGSolve reload, curved-map Jacobian sampling,
tetrahedral topology, CAD-sidecar metadata, and positive required material /
boundary measures. Used as the last gate after ``cubit-plugin-install`` to
prove the full round-trip is solver-ready rather than merely file-producing.

Usage::

    cubit-smoke-test                  # uses ih_bem_sample.jou
    cubit-smoke-test --order 3        # mesh order (default 2)
    cubit-smoke-test --keep           # keep temp workdir for inspection
    cubit-smoke-test --jou X.jou      # override sample file
    cubit-smoke-test --expect src sink sibc   # override required labels

Exit codes:
    0 -- the exported .vol passes the complete solver-ready gate
    1 -- any step failed (Cubit not found, .jou missing, export failed,
        or the solver-ready gate failed)
"""

from __future__ import annotations

import argparse
import json
import math
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
    lines = Path(vol_path).read_text(
        encoding="utf-8", errors="replace").splitlines()
    try:
        cursor = next(
            index for index, line in enumerate(lines)
            if line.strip() == section
        ) + 1
    except StopIteration as exc:
        raise ValueError(f".vol section {section!r} is missing") from exc

    while cursor < len(lines) and not lines[cursor].strip():
        cursor += 1
    if cursor >= len(lines):
        raise ValueError(f".vol section {section!r} has no item count")
    try:
        expected_count = int(lines[cursor].strip())
    except ValueError as exc:
        raise ValueError(
            f".vol section {section!r} has invalid item count: "
            f"{lines[cursor]!r}"
        ) from exc
    if expected_count < 0:
        raise ValueError(f".vol section {section!r} has negative item count")

    names = []
    identifiers = set()
    cursor += 1
    while len(names) < expected_count:
        if cursor >= len(lines):
            raise ValueError(
                f".vol section {section!r} is truncated: expected "
                f"{expected_count} rows, found {len(names)}"
            )
        row = lines[cursor].strip()
        cursor += 1
        if not row:
            continue
        parts = row.split(None, 1)
        if len(parts) != 2:
            raise ValueError(
                f".vol section {section!r} is truncated or has a malformed "
                f"row: {row!r}"
            )
        try:
            identifier = int(parts[0])
        except ValueError as exc:
            raise ValueError(
                f".vol section {section!r} has non-integer id: {row!r}"
            ) from exc
        name = parts[1].strip().strip('"')
        if identifier in identifiers:
            raise ValueError(
                f".vol section {section!r} repeats id {identifier}"
            )
        if not name:
            raise ValueError(
                f".vol section {section!r} has an empty name for id "
                f"{identifier}"
            )
        identifiers.add(identifier)
        names.append(name)
    return names


def _read_vol_bcnames(vol_path: Path):
    return _read_vol_named_section(vol_path, "bcnames")


def _read_vol_materials(vol_path: Path):
    return _read_vol_named_section(vol_path, "materials")


def _positive_named_measures(entries, required_names, value_key, category):
    """Return diagnostics for absent, non-finite, or non-positive measures."""
    by_name = {str(entry.get("name")): entry.get(value_key) for entry in entries}
    issues = []
    for name in required_names:
        if name not in by_name:
            issues.append(f"{category} {name!r} has no measured entry")
            continue
        try:
            value = float(by_name[name])
        except (TypeError, ValueError):
            issues.append(f"{category} {name!r} has invalid {value_key}")
            continue
        if not math.isfinite(value) or value <= 0.0:
            issues.append(
                f"{category} {name!r} has non-positive {value_key}={value!r}"
            )
    return issues


def _validate_exported_vol(
    vol_path: Path,
    *,
    order: int,
    expect: list[str],
    expect_materials: list[str],
    threshold: float,
):
    """Run the canonical production gate and return a structured summary."""
    from .check import REPORT_SCHEMA, check_consistency

    sidecar_path = Path(str(vol_path) + ".json")
    report_path = vol_path.parent / "vol-check.json"
    if not sidecar_path.is_file():
        raise RuntimeError(
            f"export did not produce required CAD sidecar {sidecar_path}"
        )
    try:
        sidecar = json.loads(sidecar_path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"could not read CAD sidecar {sidecar_path}: {exc}") from exc
    if not isinstance(sidecar, dict):
        raise RuntimeError("CAD sidecar root must be a JSON object")
    missing_sidecar_keys = sorted(
        {"materials", "boundaries", "edges", "n_elements", "n_points",
         "order", "export_time_s"} - set(sidecar)
    )
    if missing_sidecar_keys:
        raise RuntimeError(
            "CAD sidecar is missing required keys: "
            + ", ".join(missing_sidecar_keys)
        )
    try:
        export_time_s = float(sidecar["export_time_s"])
    except (TypeError, ValueError) as exc:
        raise RuntimeError("CAD sidecar export_time_s is not numeric") from exc
    if not math.isfinite(export_time_s) or export_time_s < 0.0:
        raise RuntimeError(
            f"CAD sidecar export_time_s is invalid: {export_time_s!r}"
        )

    bcnames = _read_vol_bcnames(vol_path)
    materials = _read_vol_materials(vol_path)
    issues = []
    issues.extend(
        f"missing expected boundary label: {name}"
        for name in expect if name not in bcnames
    )
    issues.extend(
        f"missing expected material label: {name}"
        for name in expect_materials if name not in materials
    )

    quality_options = {
        "min_curve_order": order,
        "require_tetrahedra": True,
    }
    # The canonical IH sample marks the entire workpiece/air interface SIBC.
    # Keep terminals out of this check by classifying only the workpiece as a
    # conductor; source/sink belong to the separate coil material.
    if {
        "workpiece", "air"
    }.issubset(expect_materials) and "sibc" in expect:
        quality_options.update({
            "conductive_materials": ("workpiece",),
            "air_materials": ("air",),
            "sibc_boundaries": ("sibc",),
            "require_all_sibc_labeled": True,
        })

    report = check_consistency(
        vol_path,
        threshold=threshold,
        strict_labels=True,
        required_materials=tuple(expect_materials),
        required_boundaries=tuple(expect),
        **quality_options,
    )
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    if report.get("schema") != REPORT_SCHEMA:
        issues.append(
            f"unexpected check-vol report schema: {report.get('schema')!r}"
        )
    if report.get("passed") is not True:
        issues.extend(str(value) for value in report.get("warnings", []))

    mesh = report.get("mesh", {})
    quality = report.get("quality") or {}
    ownership = report.get("boundary_domain_ownership") or {}
    cad = report.get("cad_reference") or {}
    if cad.get("available") is not True:
        issues.append("check-vol did not auto-discover the CAD sidecar")
    metadata = cad.get("metadata_checks") or {}
    for key in ("n_elements", "n_points", "order"):
        check = metadata.get(key)
        if not isinstance(check, dict) or check.get("passed") is not True:
            issues.append(f"CAD sidecar metadata check failed or missing: {key}")
    if mesh.get("curve_order") != order:
        issues.append(
            f"loaded .vol curve order {mesh.get('curve_order')!r} != requested {order}"
        )
    if int(mesh.get("n_elements", 0)) <= 0 or int(mesh.get("n_points", 0)) <= 0:
        issues.append("loaded .vol has no volume elements or points")
    if quality.get("passed") is not True:
        issues.append("curved-map quality gate did not pass")
    if quality.get("tetrahedron_count") != quality.get("volume_element_count"):
        issues.append("exported volume mesh contains non-tetrahedral elements")
    if int(quality.get("mapping_sample_count", 0)) <= 0:
        issues.append("curved-map quality gate sampled no Jacobians")
    if order >= 2 and int(quality.get("mapping_sample_count", 0)) <= int(
        quality.get("volume_element_count", 0)
    ):
        issues.append("high-order mapping was not sampled beyond one point per element")
    if quality.get("invalid_jacobian_sample_count") != 0:
        issues.append("invalid or orientation-flipping Jacobian samples were found")
    if ownership.get("passed") is not True:
        issues.append("DomainIn/DomainOut boundary ownership gate did not pass")
    if ownership.get("unreferenced_volume_domains"):
        issues.append(
            "volume domains absent from boundary ownership: "
            f"{ownership.get('unreferenced_volume_domains')}"
        )
    if ownership.get("duplicate_surface_connectivity_rows"):
        issues.append(
            "duplicate exported surface connectivity rows: "
            f"{ownership.get('duplicate_surface_connectivity_rows')}"
        )

    issues.extend(_positive_named_measures(
        report.get("materials", []), expect_materials, "ng_volume", "material"
    ))
    issues.extend(_positive_named_measures(
        report.get("boundaries", []), expect, "ng_area", "boundary"
    ))
    issues = list(dict.fromkeys(issues))
    return {
        "passed": not issues,
        "issues": issues,
        "bcnames": bcnames,
        "materials": materials,
        "sidecar_path": sidecar_path,
        "report_path": report_path,
        "report": report,
    }


def run_smoke_test(*, jou: str = "", order: int = 2,
                    expect: list[str] | None = None,
                    expect_materials: list[str] | None = None,
                    keep: bool = False, timeout: float = 600.0,
                    threshold: float = 1.0) -> int:
    """Run the round-trip. Returns 0 on success, 1 on any failure.

    expect: boundary (sideset) labels required in the .vol bcnames section.
    expect_materials: block (volume / material) names required in the .vol
        materials section. Block names live in a different section than
        sidesets -- do not mix them up (2026-04-14 false-FAIL).
    """
    if expect is None:
        expect = ["source", "sink", "sibc", "coil_surface", "outer"]
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

    temp_root = Path(r"C:\temp") if os.name == "nt" else Path(
        tempfile.gettempdir()
    )
    temp_root.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="cubit-smoke-", dir=temp_root))
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

    try:
        validation = _validate_exported_vol(
            vol_path,
            order=order,
            expect=expect,
            expect_materials=expect_materials,
            threshold=threshold,
        )
    except Exception as e:
        error_path = work / "vol-check-error.json"
        error_path.write_text(
            json.dumps({
                "passed": False,
                "error_type": type(e).__name__,
                "error": str(e),
                "vol_file": str(vol_path),
            }, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        print(f"[FAIL] solver-ready .vol validation raised: "
              f"{type(e).__name__}: {e}")
        print(f"       error report: {error_path}")
        print(f"       workdir retained for diagnosis: {work}")
        return 1

    report = validation["report"]
    quality = report["quality"]
    ownership = report["boundary_domain_ownership"]
    print(f"  bcnames  ({len(validation['bcnames'])}): "
          f"{validation['bcnames']}")
    print(f"  materials ({len(validation['materials'])}): "
          f"{validation['materials']}")
    print(f"  sidecar: {validation['sidecar_path']}")
    print(f"  check-vol report: {validation['report_path']}")
    print(f"  loaded mesh: order={report['mesh']['curve_order']}, "
          f"points={report['mesh']['n_points']}, "
          f"elements={report['mesh']['n_elements']}")
    print(f"  Jacobians: samples={quality['mapping_sample_count']}, "
          f"min_abs={quality['minimum_absolute_jacobian']!r}, "
          f"min_scaled={quality['minimum_scaled_jacobian']!r}, "
          f"invalid={quality['invalid_jacobian_sample_count']}")
    print(f"  boundary ownership: exterior="
          f"{ownership['exterior_surface_element_count']}, "
          f"interfaces={ownership['internal_interface_element_count']}, "
          f"unreferenced_domains={ownership['unreferenced_volume_domains']}")

    if not validation["passed"]:
        print("[FAIL] exported .vol is not solver-ready:")
        for issue in validation["issues"]:
            print(f"       - {issue}")
        print(f"       workdir retained for diagnosis: {work}")
        return 1

    print()
    print("[OK] round-trip solver-ready. "
          f"Boundaries {expect} + materials {expect_materials}; "
          "strict labels, CAD metadata, topology ownership, NGSolve reload, "
          "and curved-map quality all passed.")

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
                    "netgen -> complete check-vol solver-ready gate.")
    parser.add_argument("--jou", default="",
                        help="override the source .jou (default: "
                             "radia/panels/samples/ih_bem_sample.jou)")
    parser.add_argument("--order", type=int, default=2,
                        help="mesh curving order passed to export "
                             "netgen (default 2)")
    parser.add_argument("--expect", nargs="+",
                        default=["source", "sink", "sibc",
                                 "coil_surface", "outer"],
                        help="boundary (sideset) labels required in the "
                             ".vol bcnames section "
                             "(default: source sink sibc coil_surface outer)")
    parser.add_argument("--expect-materials", nargs="+",
                        default=["coil", "workpiece", "air"],
                        help="block / material names required in the .vol "
                             "materials section "
                             "(default: coil workpiece air)")
    parser.add_argument("--keep", action="store_true",
                        help="keep the temp work directory for inspection")
    parser.add_argument("--timeout", type=float, default=600.0,
                        help="seconds to wait for Cubit (default 600)")
    parser.add_argument("--threshold", type=float, default=1.0,
                        help="maximum CAD/mesh relative error in percent "
                             "(default 1.0)")
    args = parser.parse_args()

    print("=" * 60)
    print("  cubit-mesh-export: smoke test")
    print("=" * 60)
    print()
    rc = run_smoke_test(jou=args.jou, order=args.order,
                         expect=args.expect,
                         expect_materials=args.expect_materials,
                         keep=args.keep,
                         timeout=args.timeout,
                         threshold=args.threshold)
    raise SystemExit(rc)


if __name__ == "__main__":
    main()
