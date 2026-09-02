"""ci-monitor/monitor.py -- watch GitHub Actions CI runs to completion +
auto-fetch failing-step log tail.

See SKILL.md for usage.

CLI surface:

    python monitor.py [run_id ...]       # watch explicit IDs
    python monitor.py --auto N           # auto-discover N most recent
    python monitor.py --branch <ref> --auto N
    python monitor.py --poll <seconds>   # custom poll interval (default 30)
    python monitor.py --tail <N>         # log tail size on failure (default 80)

Exit codes:
    0 = every watched run ended in `success`
    1 = at least one watched run ended in `failure` / `cancelled` /
        `timed_out` / unknown conclusion
    2 = monitor itself errored (e.g. gh CLI not installed, JSON parse
        failure)

Output:
    - One line per state change per run, with timestamp.
    - On exit, one block per failed run with the failing-step log tail
      (via ``gh run view <id> --log-failed``).

Both stdout and stderr are piped through to the caller; stderr is used
only for monitor-level errors.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime
from typing import Iterable


def _gh(args: list[str], timeout: int = 60) -> tuple[int, str, str]:
    """Run gh CLI command, return (returncode, stdout, stderr)."""
    try:
        cp = subprocess.run(
            ["gh"] + args, capture_output=True, text=True,
            encoding="utf-8", errors="replace",
            timeout=timeout, check=False)
        return cp.returncode, cp.stdout, cp.stderr
    except FileNotFoundError:
        print("ERROR: 'gh' CLI not found in PATH", file=sys.stderr)
        sys.exit(2)
    except subprocess.TimeoutExpired:
        return 124, "", "timeout"


def discover_runs(branch: str | None, n: int) -> list[str]:
    """Find the N most recent CI runs (optionally on a given branch).

    Returns IDs as strings, newest first.  Filters to runs that are
    queued / in_progress (still running) -- if all matching runs are
    completed, returns the N most recent regardless of status.
    """
    args = ["run", "list", "--limit", str(max(n * 4, 10)),
            "--json", "databaseId,status,conclusion,name,headBranch"]
    if branch:
        args.extend(["--branch", branch])
    rc, out, err = _gh(args)
    if rc != 0:
        print(f"ERROR: gh run list failed: {err}", file=sys.stderr)
        sys.exit(2)
    try:
        runs = json.loads(out)
    except json.JSONDecodeError as exc:
        print(f"ERROR: cannot parse gh JSON: {exc}", file=sys.stderr)
        sys.exit(2)
    active = [r for r in runs if r["status"] != "completed"]
    pool = active if active else runs
    return [str(r["databaseId"]) for r in pool[:n]]


def get_run_state(run_id: str) -> dict:
    """Return {status, conclusion, name, headBranch} for one run, or
    raise on failure."""
    rc, out, err = _gh(["run", "view", run_id,
                          "--json", "status,conclusion,name,headBranch"])
    if rc != 0:
        return {"status": "error", "conclusion": "error",
                "name": "?", "headBranch": "?", "_err": err.strip()}
    try:
        return json.loads(out)
    except json.JSONDecodeError as exc:
        return {"status": "error", "conclusion": "error",
                "name": "?", "headBranch": "?", "_err": str(exc)}


def fmt_state(rid: str, st: dict) -> str:
    return (f"  CI run {rid} {st['name']:<10s} / "
            f"{st['headBranch']:<25s} status={st['status']:<12s} "
            f"conclusion={st['conclusion'] or '-'}")


def fetch_failing_log_tail(run_id: str, n_lines: int) -> str:
    """Run ``gh run view <id> --log-failed`` and return the last
    ``n_lines`` lines."""
    rc, out, err = _gh(["run", "view", run_id, "--log-failed"],
                        timeout=120)
    if rc != 0 and not out:
        return f"(could not fetch log: {err.strip()})"
    lines = (out or "").strip().splitlines()
    if not lines:
        return "(no failed-step log returned by gh run view --log-failed)"
    if len(lines) <= n_lines:
        return "\n".join(lines)
    return "...\n" + "\n".join(lines[-n_lines:])


def watch_runs(run_ids: Iterable[str], poll_seconds: int,
                tail_lines: int) -> int:
    """Poll ``run_ids`` until all are completed, emit per state change,
    print failing-step log tail at exit on failure.  Returns exit code
    (0 = all success, 1 = any failure)."""
    run_ids = list(run_ids)
    prev = {rid: None for rid in run_ids}
    final = {}

    while run_ids:
        for rid in list(run_ids):
            st = get_run_state(rid)
            key = (st["status"], st["conclusion"])
            if prev[rid] != key:
                ts = datetime.now().strftime("%H:%M:%S")
                print(f"[{ts}]" + fmt_state(rid, st), flush=True)
                prev[rid] = key
            if st["status"] == "completed":
                final[rid] = st
                run_ids.remove(rid)
            if st.get("status") == "error":
                # Treat persistent gh errors as terminal so we don't
                # infinite-loop on a permission / network issue.
                final[rid] = st
                run_ids.remove(rid)
        if run_ids:
            time.sleep(poll_seconds)

    print()
    print("=" * 70)
    bad = [rid for rid, st in final.items()
            if st["conclusion"] not in ("success", "skipped", "cancelled")
            or st.get("status") == "error"]
    if not bad:
        n_skipped = sum(1 for st in final.values()
                         if st["conclusion"] == "skipped")
        n_success = sum(1 for st in final.values()
                         if st["conclusion"] == "success")
        print(f"ALL GREEN ({n_success} success, {n_skipped} skipped, "
              f"0 failure)")
        return 0
    print(f"{len(bad)} run(s) FAILED:")
    for rid in bad:
        st = final[rid]
        print()
        print(f"=== FAILURE: {st.get('name','?')} / "
              f"{st.get('headBranch','?')} (id {rid}) ===")
        if st.get("status") == "error":
            print(f"  monitor error: {st.get('_err','?')}")
            continue
        print(fetch_failing_log_tail(rid, tail_lines))
    return 1


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("ids", nargs="*",
                    help="CI run IDs (database IDs from gh run list).  "
                         "If omitted, the script auto-discovers via "
                         "--auto / --branch.")
    p.add_argument("--auto", type=int, default=0,
                    help="auto-discover the N most recent runs "
                         "(prefers in_progress/queued).  Ignored if "
                         "explicit IDs are given.")
    p.add_argument("--branch", default=None,
                    help="restrict --auto discovery to runs of this "
                         "branch / tag ref.")
    p.add_argument("--poll", type=int, default=30,
                    help="poll interval in seconds (default 30)")
    p.add_argument("--tail", type=int, default=80,
                    help="failing-step log tail size in lines "
                         "(default 80)")
    args = p.parse_args()

    if args.ids:
        run_ids = args.ids
    elif args.auto > 0:
        run_ids = discover_runs(args.branch, args.auto)
        if not run_ids:
            print("No runs found.", file=sys.stderr)
            return 2
    else:
        # Default: auto-discover 3 most recent.
        run_ids = discover_runs(None, 3)
        if not run_ids:
            print("No runs found.", file=sys.stderr)
            return 2

    print(f"Watching {len(run_ids)} CI run(s): {' '.join(run_ids)}")
    print(f"Poll interval: {args.poll}s, log tail on failure: "
          f"{args.tail} lines")
    print()
    return watch_runs(run_ids, args.poll, args.tail)


if __name__ == "__main__":
    sys.exit(main())
