"""Measure the cold-import and tools/list cost of radia-mcp servers.

This is validation evidence, not a fast regression test. Compare an untouched
baseline checkout with a candidate checkout and save the machine-readable JSON
next to this script.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
from typing import Any


DEFAULT_SERVERS = (
    "accelerator",
    "build123d",
    "cubit",
    "differential_forms",
    "fem",
    "force",
    "gmsh",
    "ih",
    "litz_transmission",
    "maglev",
    "magnetic_materials",
    "matlab",
    "motor",
    "radia_ngsolve",
    "topology_optimization",
)

# Public tools removed with the retired standalone-panel workflow. Keep this
# list explicit so the validation still fails for every unplanned API loss.
RETIRED_BASELINE_TOOLS = {
    "radia_ngsolve": frozenset(
        {
            "panel_add_param",
            "panel_describe_jp",
            "panel_gui_pitfalls",
            "panel_schema",
            "panel_widget_locations",
            "standalone_panels",
        }
    ),
}

_PROBE = r"""
import json
import os
import time

started = time.perf_counter()
module = __import__(
    f"radia_mcp.{os.environ['RADIA_MCP_BENCH_SERVER']}.server",
    fromlist=["mcp"],
)
import_seconds = time.perf_counter() - started
tools = module.mcp._tool_manager._tools
schemas = [
    {
        "name": name,
        "description": tool.description,
        "parameters": tool.parameters,
    }
    for name, tool in sorted(tools.items())
]
print(json.dumps({
    "import_seconds": import_seconds,
    "tool_count": len(tools),
    "schema_bytes": len(json.dumps(
        schemas, default=str, ensure_ascii=False, separators=(",", ":")
    ).encode("utf-8")),
    "tool_names": sorted(tools),
}))
"""


def _git_revision(source: Path) -> str | None:
    worktree = next(
        (parent for parent in (source, *source.parents) if (parent / ".git").exists()),
        source,
    )
    completed = subprocess.run(
        [
            "git",
            "-c",
            f"safe.directory={worktree.as_posix()}",
            "-C",
            str(source),
            "rev-parse",
            "HEAD",
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def _probe(source: Path, server: str, profile: str) -> dict[str, Any]:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(source)
    env["RADIA_MCP_TOOL_PROFILE"] = profile
    env["RADIA_MCP_BENCH_SERVER"] = server
    completed = subprocess.run(
        [sys.executable, "-c", _PROBE],
        capture_output=True,
        text=True,
        env=env,
        timeout=180,
    )
    if completed.returncode:
        raise RuntimeError(
            f"{server}/{profile} probe failed:\n{completed.stderr[-2000:]}"
        )
    return json.loads(completed.stdout.strip().splitlines()[-1])


def _measure(
    source: Path,
    server: str,
    profile: str,
    repeats: int,
) -> dict[str, Any]:
    rows = [_probe(source, server, profile) for _ in range(repeats)]
    first = rows[0]
    if any(row["tool_names"] != first["tool_names"] for row in rows[1:]):
        raise RuntimeError(f"non-deterministic tools/list for {server}/{profile}")
    return {
        "profile": profile,
        "import_seconds_median": statistics.median(
            row["import_seconds"] for row in rows
        ),
        "import_seconds_samples": [row["import_seconds"] for row in rows],
        "tool_count": first["tool_count"],
        "schema_bytes": first["schema_bytes"],
        "tool_names": first["tool_names"],
    }


def build_report(
    baseline_src: Path,
    candidate_src: Path,
    servers: tuple[str, ...],
    repeats: int,
) -> dict[str, Any]:
    measurements: dict[str, Any] = {}
    all_supported_compatible = True
    all_retired_absent = True
    all_not_larger = True
    any_smaller = False
    import_nonregression = True
    for server in servers:
        baseline = _measure(baseline_src, server, "core", repeats)
        core = _measure(candidate_src, server, "core", repeats)
        full = _measure(candidate_src, server, "full", repeats)
        baseline_names = set(baseline["tool_names"])
        full_names = set(full["tool_names"])
        retired_names = RETIRED_BASELINE_TOOLS.get(server, frozenset())
        supported_baseline_names = baseline_names - retired_names
        unexpected_missing = sorted(supported_baseline_names - full_names)
        retired_still_present = sorted(retired_names & full_names)
        compatible = not unexpected_missing
        retired_absent = not retired_still_present
        not_larger = core["schema_bytes"] <= baseline["schema_bytes"]
        smaller = core["schema_bytes"] < baseline["schema_bytes"]
        startup_ok = (
            core["import_seconds_median"]
            <= baseline["import_seconds_median"] * 1.25
        )
        all_supported_compatible &= compatible
        all_retired_absent &= retired_absent
        all_not_larger &= not_larger
        any_smaller |= smaller
        import_nonregression &= startup_ok
        measurements[server] = {
            "baseline": baseline,
            "candidate_core": core,
            "candidate_full": full,
            "tool_count_reduction_percent": round(
                100.0 * (baseline["tool_count"] - core["tool_count"])
                / baseline["tool_count"],
                2,
            ),
            "schema_size_reduction_percent": round(
                100.0 * (baseline["schema_bytes"] - core["schema_bytes"])
                / baseline["schema_bytes"],
                2,
            ),
            "cold_import_reduction_percent": round(
                100.0
                * (
                    baseline["import_seconds_median"]
                    - core["import_seconds_median"]
                )
                / baseline["import_seconds_median"],
                2,
            ),
            "retired_baseline_tools": sorted(retired_names),
            "unexpected_missing_full_tools": unexpected_missing,
            "retired_tools_still_present": retired_still_present,
            "full_profile_preserves_supported_baseline_tools": compatible,
            "retired_baseline_tools_absent": retired_absent,
            "core_schema_is_smaller": smaller,
            "core_schema_is_not_larger": not_larger,
            "cold_import_within_25_percent_of_baseline": startup_ok,
        }
    passed = (
        all_supported_compatible
        and all_retired_absent
        and all_not_larger
        and any_smaller
        and import_nonregression
    )
    return {
        "schema": "radia.validation.mcp-server-surface.v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "machine": platform.node(),
        "platform": platform.platform(),
        "python": sys.version,
        "repeats": repeats,
        "baseline": {
            "git_revision": _git_revision(baseline_src),
        },
        "candidate": {
            "git_revision": _git_revision(candidate_src),
        },
        "measurements": measurements,
        "acceptance": {
            "full_profile_preserves_all_supported_baseline_tools": (
                all_supported_compatible
            ),
            "retired_baseline_tools_absent": all_retired_absent,
            "core_schema_is_not_larger_for_every_server": all_not_larger,
            "core_schema_is_smaller_for_at_least_one_server": any_smaller,
            "cold_import_nonregression": import_nonregression,
            "passed": passed,
        },
    }


def main() -> int:
    repo_root = Path(__file__).parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-src", type=Path, required=True)
    parser.add_argument(
        "--candidate-src",
        type=Path,
        default=repo_root / "packages" / "radia-mcp" / "src",
    )
    parser.add_argument("--servers", nargs="+", default=DEFAULT_SERVERS)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("results_mcp_server_surface_20260902.json"),
    )
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    report = build_report(
        args.baseline_src.resolve(),
        args.candidate_src.resolve(),
        tuple(args.servers),
        args.repeats,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report["acceptance"], indent=2, sort_keys=True))
    print(args.output)
    return 0 if report["acceptance"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
