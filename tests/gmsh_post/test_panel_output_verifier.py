"""Guards: the panel helper delegates strict format and field checks.

Assumes: the committed analytic post-only fixture is independently verified.
"""
from pathlib import Path

from validation_test.panels.verify_gmsh_output import verify_msh


def test_panel_verifier_uses_strict_current_format(tmp_path):
    fixture = Path(__file__).resolve().parents[2] / "validation_test/gmsh_post/fixtures/unit_cube_post.msh"
    assert verify_msh(fixture, verbose=False), fixture
    for version in ((4, 1), (2, 2)):
        rejected = tmp_path / "unsupported.msh"
        header = ".".join(map(str, version))
        rejected.write_text(fixture.read_text(encoding="utf-8").replace("4.2 0 8", f"{header} 0 8", 1), encoding="utf-8")
        assert not verify_msh(rejected, verbose=False), header
