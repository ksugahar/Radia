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
        "# Autonomous Basic Learning Pass",
        "",
        f"- pass: `{artifact['pass']}`",
        f"- slots processed: `{summary['slot_count']}`",
        f"- basic rows: `{summary['row_count']}`",
        f"- source-tool candidates: `{summary['source_tool_candidate_count']}`",
        f"- row gate: `{artifact['row_gate']['status']}`",
        f"- feedback gate: `{artifact['mcp_feedback']['feedback_gate']['status']}`",
        "",
        "## Family Counts",
        "",
        "| Family | Count |",
        "| --- | ---: |",
    ]
    for family, count in summary["family_counts"].items():
        lines.append(f"| {family} | {count} |")
    lines.extend(
        [
            "",
            "## Lane Counts",
            "",
            f"- public: `{summary['public_lane_counts']}`",
            f"- source-tool: `{summary['source_tool_lane_counts']}`",
            "",
            "This is a basic-learning pass.  It verifies every queued seed has a",
            "public analogue row and a next action, but it does not claim every",
            "commercial/live source-tool solver has executed.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run an autonomous basic-learning pass over a source-native queue.")
    parser.add_argument("--queue-json", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--artifact-id", default="autonomous_basic_learning_160")
    parser.add_argument("--out-json", type=Path, default=None)
    parser.add_argument("--out-md", type=Path, default=None)
    parser.add_argument("--no-check-local-sources", action="store_true")
    args = parser.parse_args()

    repo_root = _repo_root_from_script()
    sys.path.insert(0, str(repo_root / "packages/radia-mcp/src"))

    from radia_mcp import __version__ as radia_mcp_version
    from radia_mcp.radia_ngsolve.loop_autolearn import (
        build_autonomous_basic_learning_artifact,
        utc_now,
    )

    queue = json.loads(args.queue_json.read_text(encoding="utf-8"))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_json = args.out_json or args.out_dir / f"{args.artifact_id}.json"
    out_md = args.out_md or args.out_dir / f"{args.artifact_id}.md"
    command = " ".join(sys.argv)
    artifact = build_autonomous_basic_learning_artifact(
        queue,
        artifact_id=out_json.stem,
        queue_id=args.queue_json.stem,
        run_date_utc=utc_now(),
        radia_mcp_version=radia_mcp_version,
        command=f"python {command}",
        check_local_sources=not args.no_check_local_sources,
    )
    out_json.write_text(json.dumps(artifact, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    _write_markdown(out_md, artifact)
    print(
        json.dumps(
            {
                "pass": artifact["pass"],
                "artifact": str(out_json),
                "slots": artifact["summary"]["slot_count"],
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
