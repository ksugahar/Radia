"""Source contract for portable Mathematica derivation artifacts, no CAS run."""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DERIVATIONS = ROOT / "validation_test/maglev/research_cln/axifem"


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
