"""Every decoration template must draw as itself.

The renderer chose an accent glyph with an if-chain that fell through to the
hat for anything it did not list, and it listed five of nineteen
embellishments.  So the Decoration palette's prime, double prime, triple
prime, double dot, triple dot, cancel, frown and smile all drew a hat, while
the LaTeX and the Office MathML they copied were correct.  A wrong picture
beside a right paste is the worst shape for this defect: nothing downstream
disagrees, so only a person looking at the screen can catch it -- and one did,
after the release.

These tests compare renderings against each other instead of against stored
pixels, so they stay meaningful when metrics or fonts change.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "build"))

eqnedit_core = pytest.importorskip("eqnedit_core")

# The Decoration palette, in the order palettes.cpp declares it.
DECORATIONS = [
    "hat", "tilde", "bar", "vec", "dot", "ddot", "dddot",
    "prime", "dprime", "tprime", "strike", "frown", "smile",
]


def render(kind: str) -> str:
    """Put one decoration around the same base and return its SVG."""
    equation = eqnedit_core.Equation()
    equation.insert_template(kind)
    equation.insert_text("a")
    return equation.svg()


@pytest.fixture(scope="module")
def rendered():
    return {kind: render(kind) for kind in DECORATIONS}


def test_every_decoration_renders(rendered):
    for kind, svg in rendered.items():
        assert svg.strip(), f"{kind} rendered nothing"


def test_no_two_decorations_render_alike(rendered):
    """The fall-through made eight of these identical to the hat."""
    seen: dict[str, str] = {}
    collisions = []
    for kind, svg in rendered.items():
        if svg in seen:
            collisions.append(f"{kind} renders exactly like {seen[svg]}")
        seen[svg] = kind
    assert not collisions, "; ".join(collisions)


@pytest.mark.parametrize("kind", ["prime", "dprime", "tprime"])
def test_prime_is_a_suffix_not_an_accent(kind, rendered):
    """TeX sets a' to the right of the letter, not centred above it.

    A suffix widens the line; an accent does not.  Comparing against the hat,
    which is a true accent over the same base, keeps this independent of the
    font's actual advance widths.
    """
    assert rendered[kind] != rendered["hat"]
    plain = render_plain()
    assert width_of(rendered[kind]) > width_of(plain), (
        f"{kind} did not advance past the base, so it was drawn as an accent")


def render_plain() -> str:
    equation = eqnedit_core.Equation()
    equation.insert_text("a")
    return equation.svg()


def width_of(svg: str) -> float:
    """Read the width attribute the SVG header carries."""
    import re
    match = re.search(r'width="([0-9.]+)', svg)
    assert match, "the SVG carries no width"
    return float(match.group(1))


def test_the_palette_and_the_renderer_agree_on_what_exists():
    """A template the palette offers must be one insert_template accepts."""
    known = set(eqnedit_core.Equation.templates())
    missing = [kind for kind in DECORATIONS if kind not in known]
    assert not missing, f"palette offers templates the model rejects: {missing}"
