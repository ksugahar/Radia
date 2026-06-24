---
name: ci-monitor
description: Watch one or more GitHub Actions CI runs to completion, then on failure auto-fetch the failing job's log tail so the AI agent can diagnose without manually clicking through `gh run view`. Use this skill (a) immediately after `git push` of a release tag (the release-qud flow's Phase 7), (b) after pushing a fix-forward commit on main, (c) any time you want to know when N CI runs finish AND get the failing-step output in one go. Default: poll every 30s, emit a state-change line per run, exit when all runs are completed (success or failure), then for each failure, fetch the last ~80 lines of the failed step's log via `gh run view --log-failed`. CI run IDs can be passed explicitly or auto-discovered from the most recent N runs of the current repo.
---

# ci-monitor

Single skill, two phases:

1. **Watch**: poll `gh run view <id>` for one or more CI run IDs and
   emit a state-change line per run when status / conclusion changes.
   Exits when every watched run is `completed`.
2. **Triage on failure**: for any run that ended in `failure`,
   automatically pull the failing-step log tail with
   `gh run view <id> --log-failed` and print to stdout, so the
   AI agent's chat thread has the diagnosis context immediately.

## When to use

- **Right after `git push` of a release** (Phase 7 of `release-qud`):
  three CI runs trigger (main + tag1 + tag2), and we need to know
  which finish, which fail, and why before proceeding to deploy.
- **After a fix-forward commit** on a still-open release: confirm
  the fix actually fixed the CI, not just compiled locally.
- **When the user says "is CI green yet?"** for any push.
- **When a test you just edited landed on main**: poll once instead
  of asking the user to `gh run watch` manually.

## When NOT to use

- For a single fast (~30s) CI run you can just `gh run watch` once.
- For local test runs (use pytest directly).
- If you need the FULL log of a failure (not just the last failing
  step) — fall back to `gh run view <id> --log` for the complete
  job log; this skill only fetches the failing-step tail.

## Quick reference

Default mode: auto-discover the 3 most recent in_progress / queued
runs and watch them.

```bash
python .Codex/skills/ci-monitor/monitor.py
```

Explicit run IDs (release-qud Phase 7 idiom):

```bash
python .Codex/skills/ci-monitor/monitor.py 25275163545 25275166586 25275166607
```

Custom poll interval (default 30 s):

```bash
python .Codex/skills/ci-monitor/monitor.py --poll 60 25275166586
```

Watch the most recent N runs of a specific branch / tag:

```bash
python .Codex/skills/ci-monitor/monitor.py --branch v4.27.0 --auto 1
python .Codex/skills/ci-monitor/monitor.py --branch main   --auto 3
```

## Output format

Per state change (per run):

```
[16:14:23]  CI / v4.27.0       status=in_progress   conclusion=-
[16:18:41]  CI / v4.27.0       status=completed     conclusion=success
```

When ALL watched runs are completed, exits with:
- 0 if every run was `success`
- 1 if any run was `failure` / `cancelled` / `timed_out`

After exit, on any failure run, the script prints:

```
=== FAILURE: CI / v4.27.0 (id 25275166586) ===
[last ~80 lines of the failing step from gh run view --log-failed]
```

so the calling AI agent can read the failure cause without a
follow-up tool call.

## Implementation: monitor.py

The runner is `monitor.py` in this skill directory.  See its
docstring for the full CLI surface.

## Integration with release-qud

`release-qud` Phase 7 is the canonical caller.  The release-qud flow
ends with the main/tag CI runs; call
this skill with their IDs once they're queued, and the next thing
you read in the chat is either "all green, proceed to phase 8" or
"failed: here's the log tail, fix and re-push".
