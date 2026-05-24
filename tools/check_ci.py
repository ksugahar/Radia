"""Check latest GitHub Actions CI runs without gh CLI.

Uses the public GitHub REST API (no auth required for public repos,
60 req/hour limit). Designed to replace `gh run list --limit N` in
the deploy skill after gh CLI was retired (2026-05-24).

Usage:
    python tools/check_ci.py                    # last 5 runs, any branch
    python tools/check_ci.py --branch main      # filter by branch
    python tools/check_ci.py --branch claude/foo --limit 3
    python tools/check_ci.py --watch            # poll until newest is done
    python tools/check_ci.py --sha 1234567      # check a specific sha

Exit code:
    0 — all checked runs are SUCCESS (or no runs found)
    1 — at least one FAILURE / CANCELLED in the checked set
    2 — at least one still IN_PROGRESS (when --watch off)
"""
from __future__ import annotations
import argparse
import json
import sys
import time
import urllib.parse
import urllib.request

REPO = "ksugahar/Radia"
API  = f"https://api.github.com/repos/{REPO}/actions/runs"


def fetch_runs(branch: str | None = None, sha: str | None = None,
                limit: int = 5) -> list[dict]:
    params: dict[str, str] = {"per_page": str(min(max(limit, 1), 100))}
    if branch:
        params["branch"] = branch
    if sha:
        params["head_sha"] = sha
    url = f"{API}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=30) as r:
        data = json.load(r)
    return data.get("workflow_runs", [])[:limit]


def fmt_run(r: dict) -> str:
    conclusion = r.get("conclusion") or "-"
    status     = r.get("status", "?")
    marker = {
        "success":   "OK ",
        "failure":   "FAIL",
        "cancelled": "CAN ",
        None:        "    ",
    }.get(r.get("conclusion"), "    ")
    if status == "in_progress":
        marker = "RUN "
    elif status == "queued":
        marker = "QUE "
    return (f"  {marker}  {r['created_at']}  "
            f"{r['name'][:38]:<38}  "
            f"sha={r['head_sha'][:7]}  "
            f"branch={r['head_branch']}  "
            f"{status}/{conclusion}")


def overall_exit(runs: list[dict]) -> int:
    if not runs:
        return 0
    failures   = [r for r in runs if r.get("conclusion") in ("failure", "cancelled", "timed_out")]
    inflight   = [r for r in runs if r.get("status") in ("in_progress", "queued")]
    if failures:
        return 1
    if inflight:
        return 2
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--branch", help="filter by head_branch")
    p.add_argument("--sha", help="filter by head_sha (7+ chars)")
    p.add_argument("--limit", type=int, default=5)
    p.add_argument("--watch", action="store_true",
                    help="poll every 30s until no in_progress / queued remains")
    args = p.parse_args()

    while True:
        runs = fetch_runs(args.branch, args.sha, args.limit)
        if not runs:
            print(f"(no runs found for branch={args.branch}, sha={args.sha})")
            return 0

        print(f"=== {len(runs)} run(s) [{REPO}] ===")
        for r in runs:
            print(fmt_run(r))

        rc = overall_exit(runs)
        if not args.watch:
            return rc
        if rc != 2:  # no in_progress / queued left
            return rc
        print("(still running — re-checking in 30s)\n")
        time.sleep(30)


if __name__ == "__main__":
    sys.exit(main())
