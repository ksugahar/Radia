"""Route durable MagLev evidence outside the public docs tree."""

from pathlib import Path


DOCS_DIR = Path(__file__).resolve().parent
VALIDATION_DIR = DOCS_DIR.parents[2] / "validation_test" / "maglev" / "demos"


def validation_output(filename, output_dir=None):
    """Route docs-owned JSON to validation while preserving copied test runs."""
    directory = DOCS_DIR if output_dir is None else Path(output_dir).resolve()
    try:
        relative = directory.relative_to(DOCS_DIR)
    except ValueError:
        target_dir = directory
    else:
        target_dir = VALIDATION_DIR / relative
    target_dir.mkdir(parents=True, exist_ok=True)
    return target_dir / filename
