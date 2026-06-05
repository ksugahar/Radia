"""Check latest GitHub Actions CI runs without gh CLI.

Uses the GitHub REST API via tools/gh_api.py, which is TOKEN-AWARE: set
$GH_TOKEN (or write ~/.radia/gh_token) once to get the authenticated
5000 req/hr limit instead of the unauthenticated 60 req/hr that polling
exhausts.  Replaces `gh run list` after gh CLI was retired (2026-05-24).

Usage:
    python tools/check_ci.py                    # last 5 runs, any branch
    python tools/check_ci.py --branch main      # filter by branch
    python tools/check_ci.py --sha 1234567      # check a specific sha
    python tools/check_ci.py --watch            # poll until newest is done
    python tools/check_ci.py --rate             # just show rate-limit status

Exit code:
    0 -- all checked runs are SUCCESS (or no runs found)
    1 -- at least one FAILURE / CANCELLED in the checked set
    2 -- at least one still IN_PROGRESS (when --watch off)
"""
from __future__ import annotations

import argparse
import datetime
import os
import sys
import time
import urllib.parse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gh_api import authenticated, gh_get, rate_limit  # noqa: E402

REPO = "ksugahar/Radia"


def fetch_runs(branch=None, sha=None, limit=5):
    # When a sha is given, filter CLIENT-SIDE by prefix (the API's head_sha
    # param requires the FULL 40-char sha, so a short sha would match
    # nothing).  Fetch a wider page in that case.
    per = 40 if sha else min(max(limit, 1), 100)
    params = {"per_page": str(per)}
    if branch:
        params["branch"] = branch
    path = f"repos/{REPO}/actions/runs?{urllib.parse.urlencode(params)}"
    data, _ = gh_get(path)
    runs = data.get("workflow_runs", [])
    if sha:
        return [r for r in runs if r["head_sha"].startswith(sha)][:10]
    return runs[:limit]


def fmt_run(r: dict) -> str:
    conclusion = r.get("conclusion") or "-"
    status = r.get("status", "?")
    marker = {"success": "OK ", "failure": "FAIL", "cancelled": "CAN ",
              None: "    "}.get(r.get("conclusion"), "    ")
    if status == "in_progress":
        marker = "RUN "
    elif status == "queued":
        marker = "QUE "
    return (f"  {marker}  {r['created_at']}  {r['name'][:38]:<38}  "
            f"sha={r['head_sha'][:7]}  branch={r['head_branch']}  "
            f"{status}/{conclusion}")


def overall_exit(runs: list[dict]) -> int:
    if not runs:
        return 0
    if any(r.get("conclusion") in ("failure", "cancelled", "timed_out") for r in runs):
        return 1
    if any(r.get("status") in ("in_progress", "queued") for r in runs):
        return 2
    return 0


def _print_rate():
    core = rate_limit()
    reset_jst = datetime.datetime.fromtimestamp(
        core["reset"], datetime.timezone.utc) + datetime.timedelta(hours=9)
    auth = "authenticated (5000/hr)" if authenticated() else "UNAUTHENTICATED (60/hr)"
    print(f"GitHub API: {auth} | remaining {core['remaining']}/{core['limit']} "
          f"| resets {reset_jst:%H:%M:%S} JST")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--branch", help="filter by head_branch")
    p.add_argument("--sha", help="filter by head_sha (7+ chars)")
    p.add_argument("--limit", type=int, default=5)
    p.add_argument("--watch", action="store_true",
                   help="poll until no in_progress / queued remains")
    p.add_argument("--rate", action="store_true",
                   help="just print rate-limit status and exit")
    args = p.parse_args()

    if args.rate:
        _print_rate()
        return 0

    # Authenticated -> safe to poll at 30s; unauthenticated -> 90s to spare
    # the 60/hr budget.
    interval = 30 if authenticated() else 90
    while True:
        runs = fetch_runs(args.branch, args.sha, args.limit)
        if not runs:
            print(f"(no runs found for branch={args.branch}, sha={args.sha})")
            return 0
        print(f"=== {len(runs)} run(s) [{REPO}] ===")
        for r in runs:
            print(fmt_run(r))
        _print_rate()

        rc = overall_exit(runs)
        if not args.watch or rc != 2:
            return rc
        print(f"(still running -- re-checking in {interval}s)\n")
        time.sleep(interval)


if __name__ == "__main__":
    sys.exit(main())
