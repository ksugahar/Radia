"""Contract: shipped knowledge and prompts do not send agents to APIs that do not exist."""

from pathlib import Path

import pytest

SOURCE = Path(__file__).resolve().parents[1] / "src" / "radia_mcp"

# Each entry was copied into recipes but never existed in Radia or SparseSolv.
STALE = (
    "ssn.BiCGStab",
    "_compute_panel_local_radii",
    "radia.panels.calc_carstensen_loss",
)


def _texts():
    for path in sorted(SOURCE.rglob("*")):
        if path.suffix in {".py", ".md"} and path.is_file():
            yield path, path.read_text(encoding="utf-8")


@pytest.mark.parametrize("token", STALE)
def test_knowledge_does_not_reference_missing_api(token):
    hits = [str(path.relative_to(SOURCE)) for path, text in _texts() if token in text]
    assert not hits, f"{token} referenced in {hits}"


def test_sparsesolv_exports_no_bicgstab():
    ssn = pytest.importorskip("radia.sparsesolv_ngsolve")
    assert not hasattr(ssn, "BiCGStab")
