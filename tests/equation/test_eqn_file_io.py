r"""Reading and writing a .eqn file, including one whose name is Japanese.

This came from pointing the converter at a real lecture file for the first
time.  read_eqn took a narrow std::ifstream, which on Windows reads the path
in the ANSI code page -- so an ASCII name worked and every file in this lab's
own tree did not.  The failure was a clean "cannot open", but the reason was
invisible from Python, and the whole retirement goal for Equation Editor runs
through this function.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

B = chr(92)

# The header Equation Editor 3 writes at the front of a raw MTEF stream:
# version, platform, product, product version, product subversion.
MTEF3 = bytes([3, 1, 1, 3, 10, 10, 4])


def test_a_written_file_reads_back(tmp_path):
    data = equation.tex_to_mtef(B + "frac{a}{b}")
    p = tmp_path / "plain.eqn"
    equation.write_eqn(str(p), data)
    assert equation.read_eqn(str(p)) == data


@pytest.mark.parametrize("name", [
    "第8回.eqn",
    "透磁率.eqn",
    "Kelvin変換.eqn",
    "with space.eqn",
])
def test_a_japanese_name_is_a_name_like_any_other(tmp_path, name):
    data = equation.tex_to_mtef(B + "vec{" + B + "bm{B}}")
    p = tmp_path / name
    equation.write_eqn(str(p), data)
    assert equation.read_eqn(str(p)) == data
    assert p.exists()


def test_a_missing_file_says_so(tmp_path):
    with pytest.raises(RuntimeError):
        equation.read_eqn(str(tmp_path / "第9回.eqn"))


def test_the_mtef_route_writes_a_vector_as_bm():
    r"""MTEF's vector typeface is a vector, so it takes \bm like every other
    path.  It wrote \mathbf, which is now the upright face -- the one route
    that had not heard the rule."""
    mtef = equation.tex_to_mtef(B + "vec{" + B + "bm{B}}")
    out = equation.mtef_to_tex(mtef)
    assert B + "bm{" in out
    assert B + "mathbf{" not in out
