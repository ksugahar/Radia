"""Console entry points for locating and checking the MATLAB distribution."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from . import __version__, matlab_path, mex_path


def _manifest() -> dict[str, object]:
    path = Path(__file__).with_name("manifest.json")
    return json.loads(path.read_text(encoding="utf-8"))


def _doctor_payload() -> dict[str, object]:
    root = matlab_path()
    mex = mex_path()
    notice = root / "THIRD_PARTY_NOTICES.md"
    notice_text = notice.read_text(encoding="utf-8") if notice.is_file() else ""
    required_notices = (
        "Copyright (c) 2018 Preferred Networks, Inc.",
        "Copyright (c) 2025 Preferred Networks, Inc.",
        (
            "Optuna, the Optuna logo and any related marks are trademarks of "
            "Preferred Networks, Inc."
        ),
        "independent, unofficial project",
    )
    notices_complete = all(token in notice_text for token in required_notices)
    matlab_files = sorted((root / "+radia" / "+optuna").rglob("*.m"))
    manifest = _manifest()
    expected_files = int(manifest["matlab_file_count"])
    simulink_entries = [root / entry for entry in manifest["simulink_entry_points"]]
    ok = (
        root.is_dir()
        and mex.is_file()
        and len(matlab_files) == expected_files
        and all(path.is_file() for path in simulink_entries)
        and notices_complete
    )
    return {
        "schema": "radia-optuna.doctor.v1",
        "ok": ok,
        "version": __version__,
        "matlab_path": str(root),
        "mex_path": str(mex),
        "mex_size_bytes": mex.stat().st_size if mex.is_file() else None,
        "mex_sha256": (
            hashlib.sha256(mex.read_bytes()).hexdigest() if mex.is_file() else None
        ),
        "matlab_file_count": len(matlab_files),
        "expected_matlab_file_count": expected_files,
        "native_command_count": manifest["native_command_count"],
        "simulink_standalone": manifest["simulink_standalone"],
        "simulink_entry_points": manifest["simulink_entry_points"],
        "simulink_entry_count": sum(path.is_file() for path in simulink_entries),
        "oracle_version": manifest["oracle_version"],
        "radia_integration_adapters": manifest["radia_integration_adapters"],
        "third_party_notice": str(notice),
        "upstream_notices_complete": notices_complete,
    }


def path_main() -> None:
    print(matlab_path())


def doctor_main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", dest="as_json")
    args = parser.parse_args()
    payload = _doctor_payload()
    if args.as_json:
        print(json.dumps(payload, indent=2))
    else:
        status = "PASS" if payload["ok"] else "FAIL"
        print(f"radia-optuna {payload['version']}: {status}")
        print(f"MATLAB path: {payload['matlab_path']}")
        print(f"MEX: {payload['mex_path']}")
        print(
            f"MATLAB files: {payload['matlab_file_count']}/"
            f"{payload['expected_matlab_file_count']}"
        )
    if not payload["ok"]:
        raise SystemExit(1)
