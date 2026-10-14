from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

ALLOWED_MARKDOWN = {
    "AGENTS.md",
    "BUILD.md",
    "CHANGELOG.md",
    "CLAUDE.md",
    "CONTRIBUTING.md",
    "MEMORY.md",
    "README.md",
    "ROADMAP.md",
    "SECURITY.md",
}

BANNED_SUFFIXES = {
    ".exp",
    ".jou",
    ".lib",
    ".log",
    ".obj",
    ".pdb",
    ".slxc",
    ".tex",
    ".vtu",
}

BANNED_DIRECTORIES = {"logs", "slprj", "temp", "work"}


def test_repository_root_contains_no_scratch_artifacts() -> None:
    violations: list[str] = []
    for path in ROOT.iterdir():
        name = path.name
        if path.is_dir() and name.lower() in BANNED_DIRECTORIES:
            violations.append(name + "/")
        elif path.is_file() and path.suffix.lower() in BANNED_SUFFIXES:
            violations.append(name)
        elif path.is_file() and path.suffix.lower() == ".md" and name not in ALLOWED_MARKDOWN:
            violations.append(name)

    assert not violations, (
        "Repository root is source space, not scratch storage. Move generated "
        f"artifacts to C:\\temp and documents under docs/: {sorted(violations)}"
    )
