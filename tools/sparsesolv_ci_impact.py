"""Skip native AMS work only for verified test/dev dependency-only changes."""

import json
import os
from pathlib import Path
import subprocess
import tomllib


def test_dependencies_only(before: str, after: str) -> bool:
    """All configuration outside the two development extras must be identical."""
    old, new = tomllib.loads(before), tomllib.loads(after)
    if old == new:
        return False
    for document in (old, new):
        extras = document.get("project", {}).get("optional-dependencies", {})
        for name in ("test", "dev"):
            values = extras.pop(name, [])
            if not isinstance(values, list) or not all(isinstance(x, str) for x in values):
                return False
    return old == new


def native_required(event_name, event, git):
    if event_name not in ("push", "pull_request"):
        return True, "manual or unknown event"
    # A PR checkout is GitHub's merge commit; its first parent is the base.
    # Pushes must inspect the entire event range, not just the final commit.
    base = "HEAD^1" if event_name == "pull_request" else event["before"]
    if not base or set(base) == {"0"}:
        return True, "no comparison base"
    changed = git("diff", "--name-only", "-z", base, "HEAD", "--").split("\0")
    if set(filter(None, changed)) != {"pyproject.toml"}:
        return True, "changes beyond pyproject.toml"
    if test_dependencies_only(git("show", f"{base}:pyproject.toml"),
                              git("show", "HEAD:pyproject.toml")):
        return False, "only test/dev optional dependencies changed"
    return True, "runtime, build, packaging, or other configuration changed"


def main():
    root = Path(__file__).resolve().parents[1]

    def git(*args):
        return subprocess.check_output(
            ["git", "-c", f"safe.directory={root.as_posix()}", "-C", str(root), *args],
            encoding="utf-8", stderr=subprocess.PIPE,
        )

    try:
        event = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text(encoding="utf-8"))
        required, reason = native_required(os.environ["GITHUB_EVENT_NAME"], event, git)
    except (OSError, ValueError, KeyError, TypeError, AttributeError, subprocess.SubprocessError) as exc:
        required, reason = True, f"impact unknown ({type(exc).__name__}); fail closed"
    print(f"Native AMS required={required}: {reason}")
    with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as output:
        output.write(f"required={str(required).lower()}\n")


if __name__ == "__main__":
    main()
