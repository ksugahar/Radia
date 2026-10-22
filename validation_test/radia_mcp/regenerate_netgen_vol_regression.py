"""Produce a separate regression-reference proposal; never rewrite approved records."""
import argparse
import difflib
import hashlib
import json
from pathlib import Path
import subprocess
import sys

# Run from the repository root with its selected MCP source on PYTHONPATH.
from validation_test.radia_mcp.netgen_vol_examples import CASES, ROOT, example_module

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reason", required=True)
    args = parser.parse_args()
    if not args.reason.strip():
        parser.error("A reference-change reason is required")
    args.out.mkdir(parents=True, exist_ok=False)
    approved = ROOT / "tests/mcp_integration/fixtures/netgen_vol_analytic"
    manifest = json.loads((approved / "manifest.json").read_text())
    rows = []
    for name in CASES:
        with example_module(name) as module:
            data = json.loads(json.dumps(module.build_summary()))
        filename = manifest["cases"][name]["file"]
        original = (approved / filename).read_bytes()
        previous = json.loads(original)
        changed = data != previous
        # Preserve approved bytes for unchanged JSON values, including signed
        # zero representation and whitespace. Never normalize the reference.
        encoded = ((json.dumps(data, indent=2, sort_keys=True) + "\n").encode()
                   if changed else original)
        (args.out / filename).write_bytes(encoded)
        delta = list(difflib.unified_diff(
            original.decode().splitlines(),
            encoded.decode().splitlines(), fromfile="approved/" + filename,
            tofile="proposed/" + filename, lineterm=""))
        rows.append({"diff": delta, "case": name, "sha256": hashlib.sha256(encoded).hexdigest(),
                     "changed": changed, "byte_changed": encoded != original,
                     "formatting_only": encoded != original and not changed})
    report = {"reason": args.reason, "source_sha": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT)),
        "python": sys.version, "conditions": "Built-in synthetic inputs, no external VOL",
        "cases": rows, "approval": "Proposal only; user approval and review required"}
    (args.out / "proposal.json").write_text(json.dumps(report, indent=2) + "\n")

if __name__ == "__main__":
    main()
