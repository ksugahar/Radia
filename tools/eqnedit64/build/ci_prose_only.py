"""Decide whether an Eqnedit64 CI event changes only prose.

Prose-only means every changed path is under docs/equation/ or is the
release-skill prose (SKILL.md, agents/openai.yaml). Those changes get a
parse/link check and no native build. Release-skill scripts are exercised by
the native tests, so they keep the native lanes. Tags, manual runs, missing
comparison bases and any error keep the native lanes (fail closed).

Writes ``native=true|false`` to $GITHUB_OUTPUT. With ``--check-links`` it
instead verifies that relative Markdown links in the prose paths resolve.
"""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PROSE_PREFIXES = ("docs/equation/",)
PROSE_FILES = {
    ".agents/skills/release-eqnedit64/SKILL.md",
    ".agents/skills/release-eqnedit64/agents/openai.yaml",
}
LINK = re.compile(r"\]\(([^)\s]+)\)")


def is_prose(path):
    return path.startswith(PROSE_PREFIXES) or path in PROSE_FILES


def git(*args):
    return subprocess.check_output(
        ["git", "-c", f"safe.directory={ROOT.as_posix()}", "-C", str(ROOT), *args],
        encoding="utf-8")


def changed_paths():
    name = os.environ.get("GITHUB_EVENT_NAME", "")
    event = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text(encoding="utf-8"))
    if name not in ("push", "pull_request") or event.get("ref", "").startswith("refs/tags/"):
        return None
    base = event.get("before") if name == "push" else event["pull_request"]["base"]["sha"]
    if not base or set(base) == {"0"}:
        return None
    if name == "pull_request":
        base = git("merge-base", base, "HEAD").strip()
    raw = git("diff", "--name-only", "-z", "--find-renames", base, "HEAD", "--")
    return [p for p in raw.split("\0") if p]


def native_required():
    try:
        changed = changed_paths()
    except Exception as exc:  # unknown impact keeps the native lanes
        return True, f"impact unknown ({type(exc).__name__}); native kept"
    if changed is None:
        return True, "release/manual/unknown event; native kept"
    if changed and all(is_prose(p) for p in changed):
        return False, "prose-only change: " + ", ".join(sorted(changed))
    return True, "native inputs changed"


def check_links():
    paths = [p for p in git("ls-files", "docs/equation", *PROSE_FILES).splitlines()
             if p.endswith(".md")]
    broken = []
    for relative in paths:
        source = ROOT / relative
        for target in LINK.findall(source.read_text(encoding="utf-8")):
            if re.match(r"[a-z]+:", target) or target.startswith("#"):
                continue
            if not (source.parent / target.split("#", 1)[0]).exists():
                broken.append(f"{relative}: {target}")
    if broken:
        print("Broken relative links:\n  " + "\n  ".join(broken))
        return 1
    print(f"OK: {len(paths)} Eqnedit64 prose files, relative links resolve")
    return 0


def main(argv):
    if "--check-links" in argv:
        return check_links()
    native, reason = native_required()
    print(f"native={str(native).lower()} ({reason})")
    with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as output:
        output.write(f"native={str(native).lower()}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
