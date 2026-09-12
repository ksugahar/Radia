"""Skip native AMS work only for verified non-native configuration changes."""

from fnmatch import fnmatchcase
import json
import os
from pathlib import Path
import subprocess
import tomllib


# Checked against both workflow trigger lists by the lightweight contract suite.
NATIVE_PATHS = (
    "src/ext/sparsesolv/**", "src/matlab/radia_mex.cpp",
    "matlab/+radia/+internal/callMex.m", "matlab/+radia/setup.m",
    "tests/matlab/test_mex_runtime_setup.m", "matlab/+radia/+sparsesolv/**",
    "matlab/+radia/+python/sparsesolv.m", "matlab/+radia/+python/electromagnetValidation.m",
    "matlab/+radia/+python/esrfExamples.m", "matlab/+radia/+python/staticElectromagnet.m",
    "matlab/python_api_parity_manifest.json", "tests/matlab/*sparsesolv*",
    "validation_test/ngsolve_matlab_parity/run_sparsesolv_parity.py",
    "src/radia/__init__.py", "CMakeLists.txt", "Build.ps1", "tools/run_test_tier.py",
    "tools/native_build_provenance.ps1",
    "tools/sparsesolv_ci_impact.py", "tests/test_sparsesolv_ci_contract.py",
    ".github/workflows/sparsesolv.yml",
)


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


def ams_manifest_unchanged(before: str, after: str) -> bool:
    def relevant(text):
        data = json.loads(text)
        if data["schema"] != "radia.test-tier-manifest.v1":
            raise ValueError("unknown manifest schema")
        profiles = data.pop("profiles")
        data.pop("impact_rules", None)  # AMS invokes the profile without --since.
        selected = {}
        name = "sparsesolv"
        while name:
            if name in selected:
                raise ValueError("cyclic AMS profile inheritance")
            profile = profiles[name]
            paths = profile.get("paths", [])
            if not isinstance(paths, list) or not all(isinstance(p, str) for p in paths):
                raise ValueError("invalid AMS paths")
            budget = profile.get("max_elapsed_seconds")
            if budget is not None and (not isinstance(budget, (int, float)) or budget <= 0):
                raise ValueError("invalid AMS time budget")
            selected[name] = profile
            name = profile.get("extends")
            if name is not None and not isinstance(name, str):
                raise ValueError("invalid AMS parent")
        data["profiles"] = selected
        paths = [path for profile in selected.values() for path in profile.get("paths", [])]
        if not paths or len(paths) != len(set(paths)):
            raise ValueError("empty or duplicate AMS test paths")
        return data
    return relevant(before) == relevant(after)


def native_required(event_name, event, git):
    if event_name not in ("push", "pull_request"):
        return True, "manual or unknown event"
    # A PR checkout is GitHub's merge commit; its first parent is the base.
    # Pushes must inspect the entire event range, not just the final commit.
    base = "HEAD^1" if event_name == "pull_request" else event["before"]
    if not base or set(base) == {"0"}:
        return True, "no comparison base"
    changed = git("diff", "--name-only", "-z", base, "HEAD", "--").split("\0")
    changed = set(filter(None, changed))
    if any(fnmatchcase(path, pattern) for path in changed for pattern in NATIVE_PATHS):
        return True, "native source, runtime boundary, or CI contract changed"
    checks = {"pyproject.toml": test_dependencies_only,
              "tests/test_tier_manifest.json": ams_manifest_unchanged}
    configs = changed.intersection(checks)
    if not configs:
        return True, "no verifiable configuration-only impact"
    for path in sorted(configs):
        if not checks[path](git("show", f"{base}:{path}"), git("show", f"HEAD:{path}")):
            return True, f"AMS-relevant configuration changed: {path}"
    return False, "configuration changes do not affect native AMS; other files are outside its trigger paths"


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
