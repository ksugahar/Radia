"""Read-only verification after the owner removes obsolete PyPI releases."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import urllib.request


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--keep", required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    report = {"schema": "eqnedit64.pypi-latest-only.v1", "keep": args.keep,
              "read_only": True, "status": "UNVERIFIED",
              "checked_at": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    code = 2
    try:
        request = urllib.request.Request(
            "https://pypi.org/pypi/eqnedit64/json", headers={"Cache-Control": "no-cache"})
        with urllib.request.urlopen(request, timeout=30) as response:
            raw = response.read()
        data = json.loads(raw)
        versions = sorted(data["releases"])
        report.update(json_sha256=hashlib.sha256(raw).hexdigest(), versions=versions,
                      retained_file_count=len(data["releases"].get(args.keep, [])))
        process = subprocess.run(
            [sys.executable, "-m", "pip", "--isolated", "--no-cache-dir", "index",
             "versions", "eqnedit64", "--index-url", "https://pypi.org/simple"],
            capture_output=True, text=True, timeout=60)
        report.update(pip_exit_code=process.returncode, pip_stdout=process.stdout,
                      pip_stderr=process.stderr)
        match = re.search(r"^Available versions:\s*(.+)$", process.stdout, re.MULTILINE)
        if process.returncode or not match:
            raise RuntimeError("pip index failed or its version list could not be read")
        pip_versions = [value.strip() for value in match.group(1).split(",")]
        report["pip_versions"] = pip_versions
        passed = (versions == [args.keep] and pip_versions == [args.keep]
                  and report["retained_file_count"] > 0
                  and data["info"]["version"] == args.keep)
        report["status"] = "PASS" if passed else "FAIL"
        code = 0 if passed else 1
    except Exception as error:
        report["error"] = str(error)
    finally:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"{report['status']}: read-only PyPI latest-only check; report={args.report}")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
