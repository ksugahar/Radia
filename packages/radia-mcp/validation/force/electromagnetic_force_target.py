from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _repo_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def _write_markdown(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    lines = [
        "# Electromagnetic Force Target Pass",
        "",
        f"- pass: `{artifact['pass']}`",
        f"- force slots: `{summary['force_slot_count']}`",
        f"- public rows: `{summary['row_count']}`",
        f"- source-tool candidates: `{summary['source_tool_candidate_count']}`",
        f"- row gate: `{artifact['row_gate']['status']}`",
        f"- feedback gate: `{artifact['mcp_feedback']['feedback_gate']['status']}`",
        "",
        "## Target Counts",
        "",
        "| Target | Count |",
        "| --- | ---: |",
    ]
    for target, count in summary["target_counts"].items():
        lines.append(f"| {target} | {count} |")
    lines.extend(
        [
            "",
            "## Source Tool Counts",
            "",
            "| Tool | Count |",
            "| --- | ---: |",
        ]
    )
    for tool, count in summary["source_tool_counts"].items():
        lines.append(f"| {tool} | {count} |")
    lines.extend(
        [
            "",
            "This target pass verifies public analytic electromagnetic-force rows.",
            "It does not claim that live commercial/source-tool solvers have run.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract force_torque_motor slots and attach public EM-force analytic gates."
    )
    parser.add_argument("--source-json", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--artifact-id", default="em_force_target")
    parser.add_argument("--out-json", type=Path, default=None)
    parser.add_argument("--out-md", type=Path, default=None)
    args = parser.parse_args()

    repo_root = _repo_root_from_script()
    sys.path.insert(0, str(repo_root / "packages/radia-mcp/src"))

    from radia_mcp import __version__ as radia_mcp_version
    from radia_mcp.radia_ngsolve.em_force_target import (
        build_em_force_target_artifact,
        utc_now,
    )

    source_artifact = json.loads(args.source_json.read_text(encoding="utf-8"))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_json = args.out_json or args.out_dir / f"{args.artifact_id}.json"
    out_md = args.out_md or args.out_dir / f"{args.artifact_id}.md"
    command = " ".join(sys.argv)
    artifact = build_em_force_target_artifact(
        source_artifact,
        artifact_id=out_json.stem,
        source_artifact_id=args.source_json.stem,
        run_date_utc=utc_now(),
        radia_mcp_version=radia_mcp_version,
        command=f"python {command}",
    )
    out_json.write_text(json.dumps(artifact, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    _write_markdown(out_md, artifact)
    print(
        json.dumps(
            {
                "pass": artifact["pass"],
                "artifact": str(out_json),
                "force_slots": artifact["summary"]["force_slot_count"],
                "rows": artifact["summary"]["row_count"],
                "source_tool_candidates": artifact["summary"]["source_tool_candidate_count"],
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    if not artifact["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
