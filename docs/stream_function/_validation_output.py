"""Route durable stream-function evidence outside the public docs tree."""

from pathlib import Path


DOCS_DIR = Path(__file__).resolve().parent
VALIDATION_DIR = DOCS_DIR.parents[1] / "validation_test" / "stream_function" / "demos"


def validation_output(filename, output_dir=None):
    """Return the validation JSON path, preserving explicit test output dirs."""
    directory = DOCS_DIR if output_dir is None else Path(output_dir)
    if directory.resolve() == DOCS_DIR:
        directory = VALIDATION_DIR
    directory.mkdir(parents=True, exist_ok=True)
    return directory / filename


def validation_json_for_basename(output_basename, canonical_filename):
    """Route a docs-default basename to validation and preserve custom basenames."""
    base = Path(output_basename)
    if base.parent.resolve() == DOCS_DIR:
        path = VALIDATION_DIR / canonical_filename
    else:
        path = base.with_suffix(".json")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path
