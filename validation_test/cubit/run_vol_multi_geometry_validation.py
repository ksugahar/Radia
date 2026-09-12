"""Run the Cubit multi-geometry validation and attest its native evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
TEST_SCRIPT = HERE / "test_vol_multi_geometry.py"
DEFAULT_RESULTS = HERE / "export_netgen_test" / "multi_geometry_results.json"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--plugin-dir",
        type=Path,
        default=Path(os.environ.get("CUBIT_PLUGIN_DIR", REPO_ROOT / "src" / "cubit_plugin" / "build-ccm")),
    )
    parser.add_argument("--log", type=Path, default=Path(r"C:\temp\cubit-multi-geometry-validation.log"))
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    args = parser.parse_args()

    plugin_dir = args.plugin_dir.resolve()
    ccm_path = plugin_dir / "cubit_mesh_export.ccm"
    if not ccm_path.is_file():
        parser.error(f"Cubit plugin not found: {ccm_path}")

    args.log.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["CUBIT_PLUGIN_DIR"] = str(plugin_dir)
    env["CUBIT_VALIDATION_LOG"] = str(args.log.resolve())
    completed = subprocess.run(
        [sys.executable, str(TEST_SCRIPT)],
        cwd=REPO_ROOT,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    args.log.write_text(completed.stdout, encoding="utf-8", newline="\n")
    console_encoding = sys.stdout.encoding or "utf-8"
    sys.stdout.write(
        completed.stdout.encode(console_encoding, errors="replace").decode(console_encoding)
    )

    evidence = json.loads(args.results.read_text(encoding="utf-8"))
    rejection = evidence["labelled_same_material_internal_surface"]
    expected_diagnostic = rejection["expected_diagnostic"]
    diagnostic_matched = expected_diagnostic in completed.stdout

    manifest_ccm = evidence["native_provenance"]["payloads"]["cubit_mesh_export.ccm"]
    actual_ccm = {
        "path": str(ccm_path),
        "sha256": _sha256(ccm_path),
        "size": ccm_path.stat().st_size,
    }
    expected_load_line = f"Loading Plugin: '{ccm_path}'"
    actual_ccm["load_line"] = expected_load_line
    actual_ccm["observed_in_stdout"] = expected_load_line in completed.stdout
    actual_ccm["matches_manifest"] = (
        actual_ccm["sha256"] == manifest_ccm["sha256"]
        and actual_ccm["size"] == manifest_ccm["size"]
    )

    rejection["diagnostic_matched_in_stdout"] = diagnostic_matched
    evidence["subprocess_returncode"] = completed.returncode
    evidence["loaded_ccm"] = actual_ccm
    evidence["all_passed"] = bool(
        evidence["all_passed"]
        and completed.returncode == 0
        and diagnostic_matched
        and actual_ccm["observed_in_stdout"]
        and actual_ccm["matches_manifest"]
    )
    args.results.write_text(
        json.dumps(evidence, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )

    if not evidence["all_passed"]:
        print(
            "VALIDATION EVIDENCE FAILED: "
            f"returncode={completed.returncode}, diagnostic={diagnostic_matched}, "
            f"loaded_ccm={actual_ccm['observed_in_stdout']}, "
            f"manifest_match={actual_ccm['matches_manifest']}",
            file=sys.stderr,
        )
        return 1
    print(f"VALIDATION EVIDENCE: {args.results} (PASS)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
