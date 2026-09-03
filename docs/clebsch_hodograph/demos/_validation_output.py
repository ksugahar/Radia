"""Shared output location for Clebsch/Hodograph validation records."""

from pathlib import Path


def validation_output(filename: str) -> Path:
    repo_root = Path(__file__).resolve().parents[3]
    output = (
        repo_root
        / "validation_test"
        / "clebsch_hodograph"
        / "demos"
        / filename
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    return output
