"""Sphere p-sweep: export netgen at orders 1-5 and check-vol every .vol.

Plays sphere_p_convergence.jou through Cubit's Python API (batch,
nographics) with the Radia command plugin loaded the same way as the
production export worker, in a child process, then checks each exported .vol with radia.vol_check (labels,
curved-map quality, finite measures). Exit 1 on any missing file or failed
check.

Usage:
    python run_test.py --command-plugin-dir <dir with radia_cubit_export.dll> [--cubit-bin <Cubit bin dir>]
"""

import argparse
import json
import math
import os
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
SPHERE_VOLUME = 4.0 / 3.0 * math.pi * 0.05 ** 3
# Flat facets at h/R = 0.4 lose about 6 % of the sphere volume, so order 1 is
# checked against a 10 % CAD gate; curved orders must meet the default 1 %.
LINEAR_THRESHOLD_PCT = 10.0
CURVED_THRESHOLD_PCT = 1.0


def find_cubit_bin():
    from radia.cubit_locator import find_cubit_bin as locate

    bin_dir = locate()
    if bin_dir is None:
        raise SystemExit("Cubit not found: set CUBIT_BIN_DIR or pass --cubit-bin")
    return Path(bin_dir)


def play_journal(cubit_bin, plugin_dir, journal):
    """Child-process body: load the fork nglib, then Cubit, then play the journal."""
    import ngsolve  # noqa: F401  -- the plugin binds to the already loaded nglib
    import netgen

    os.add_dll_directory(os.path.dirname(netgen.__file__))
    sys.path.append(str(cubit_bin))
    os.add_dll_directory(str(cubit_bin))
    os.environ["CUBIT_PLUGIN_DIR"] = str(plugin_dir)
    import cubit

    cubit.init(["cubit", "-nographics", "-nojournal", "-batch", "-noinitfile"])
    if not cubit.cmd(f'playback "{Path(journal).as_posix()}"'):
        raise SystemExit(f"Cubit journal playback failed: {journal}")
    return 0


def exact_volume_error(vol, order):
    """Relative sphere-volume error with quadrature exact for the order-p map.

    The geometry Jacobian of an order-p tetrahedron has degree 3(p-1); vol_check's
    CAD gate integrates at NGSolve's default order, which is not exact for p >= 3.
    """
    from ngsolve import CF, Integrate, Mesh

    volume = Integrate(CF(1), Mesh(str(vol)), order=3 * order)
    return volume / SPHERE_VOLUME - 1.0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cubit-bin", type=Path)
    parser.add_argument("--command-plugin-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=SCRIPT_DIR / "output")
    parser.add_argument("--play", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    cubit_bin = (args.cubit_bin or find_cubit_bin()).resolve()
    plugin_dir = args.command_plugin_dir.resolve()
    if args.play:
        return play_journal(cubit_bin, plugin_dir, args.play)
    args.out.mkdir(parents=True, exist_ok=True)
    journal = (SCRIPT_DIR / "sphere_p_convergence.jou").read_text(encoding="utf-8")
    run_jou = args.out / "_sphere_p_convergence.jou"
    run_jou.write_text(journal.replace("{DIR}", args.out.as_posix()), encoding="utf-8")
    for order in range(1, 6):
        (args.out / f"sphere_o{order}.vol").unlink(missing_ok=True)
    run = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), "--cubit-bin", str(cubit_bin),
         "--command-plugin-dir", str(plugin_dir), "--out", str(args.out), "--play", str(run_jou)],
        capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=600)
    (args.out / "cubit.log").write_text(run.stdout + run.stderr, encoding="utf-8")
    rows = []
    for order in range(1, 6):
        vol = args.out / f"sphere_o{order}.vol"
        if not vol.is_file():
            rows.append({"order": order, "passed": False, "error": "not created"})
            continue
        check = subprocess.run(
            [sys.executable, "-m", "radia.vol_check", str(vol), "--strict-labels",
             "--required-materials", "sphere", "--required-boundaries", "surface",
             "--integration-order", str(max(4, 2 * order + 2)),
             "--threshold", str(LINEAR_THRESHOLD_PCT if order == 1 else CURVED_THRESHOLD_PCT),
             "--format", "json"],
            capture_output=True, text=True)
        report = json.loads(check.stdout) if check.stdout.strip().startswith("{") else {}
        material = next((m for m in report.get("materials", []) if m.get("name") == "sphere"), {})
        rows.append({"order": order, "passed": check.returncode == 0 and report.get("passed") is True,
                     "n_elements": report.get("n_elements"), "ng_volume": material.get("ng_volume"),
                     "volume_error": exact_volume_error(vol, order)})
    for row in rows:
        print(json.dumps(row))
    failed = [row["order"] for row in rows if not row["passed"]]
    errors = [abs(row["volume_error"]) for row in rows if row.get("volume_error") is not None]
    if len(errors) != 5 or any(later >= earlier for earlier, later in zip(errors, errors[1:])):
        failed.append("volume error not strictly decreasing with order")
    print(f"cubit exit {run.returncode}; failed orders: {failed or 'none'}")
    return 1 if failed or run.returncode else 0


if __name__ == "__main__":
    sys.exit(main())
