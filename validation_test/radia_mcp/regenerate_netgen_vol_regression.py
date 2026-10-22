"""Propose prior-record updates with the common writer; never accept references."""
import argparse
import json
from pathlib import Path
import subprocess
import time

from tools.regression_lock import propose_json_regeneration
from validation_test.radia_mcp.netgen_vol_examples import CASES, ROOT, example_module
from validation_test.radia_mcp.netgen_vol_regression import approved_case


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reason", required=True)
    args = parser.parse_args(argv)
    if not args.reason.strip():
        parser.error("A reference-change reason is required")
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    if args.out.exists():
        raise FileExistsError(args.out)
    rows = []
    for name in CASES:
        approved, _, provenance = approved_case(name)
        started = time.perf_counter()
        with example_module(name) as module:
            data = json.loads(json.dumps(module.build_summary(), allow_nan=False))
        runtime = time.perf_counter() - started
        report = propose_json_regeneration(approved, actual=data, source_sha=head,
            conditions={"reference": provenance, "inputs": "Built-in synthetic inputs, no external VOL"},
            runtime_s=runtime, reason=args.reason, out=args.out / name)
        rows.append({"case": name, "proposal": f"{name}/proposal.json",
                     "changed": report["changed"], "byte_changed": report["byte_changed"]})
    # Index only; common writer owns provenance/hash/diff/exclusive file policy.
    (args.out / "index.json").write_text(json.dumps({"cases": rows,
        "approval": "Proposal only; user approval and review required"}, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
