"""Source contract for portable Mathematica derivation artifacts, no CAS run."""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DERIVATIONS = ROOT / "validation_test/maglev/research_cln/axifem"


def test_axifem_reference_links_resolve_inside_repository():
    document = ROOT / "docs/axifem/AXIFEM.md"
    source = document.read_text(encoding="utf-8")
    assert "W%3A" not in source and "W:/30_CauerLadderNetwork" not in source
    for filename in ("bem_disk_axisym_cauer.wls", "disk_bem_cauer.py"):
        links = re.findall(r"\[`" + re.escape(filename) + r"`\]\(([^)]+)\)", source)
        assert len(links) == 1
        target = (document.parent / links[0]).resolve()
        assert target.is_relative_to(ROOT.resolve())
        assert target.is_file()


@pytest.mark.parametrize("order,expected", [(1, 1), (2, 2)])
def test_henrotte_json_exports_are_relative_to_the_script(order, expected):
    source = (DERIVATIONS / f"derive_quad_q{order}_henrotte.wls").read_text(
        encoding="utf-8"
    )
    exports = re.findall(r"^Export\[.*$", source, re.MULTILINE)
    assert len(exports) == expected
    assert "DirectoryName[$InputFileName]" in source
    assert all(line.startswith("Export[FileNameJoin[{") for line in exports)
    assert all('"JSON"]' in line for line in exports)
    assert not re.search(r"[A-Za-z]:[/\\]", "\n".join(exports))
    if order == 2:
        assert "exportDir = DirectoryName[$InputFileName];" in source
    for line in exports:
        root = "DirectoryName[$InputFileName]" if order == 1 else "exportDir"
        assert line.startswith(f"Export[FileNameJoin[{{{root}, ")
