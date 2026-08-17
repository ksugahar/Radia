"""Which spans of a Markdown file are equations, and which are not.

The scan has to be conservative in one specific direction.  Missing an equation
costs the user one manual edit; *inventing* one rewrites their prose or their
code when the file is saved, silently.  Most of the cases below are therefore
about what must NOT be treated as math.

The other guarantee is exactness: loading and saving an untouched file has to
reproduce it byte for byte, or the editor cannot be trusted with a file it did
not create.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
MarkdownDoc = equation.MarkdownDoc
MdSegment = equation.MdSegment


DETECTION = [
    ("no math here at all", []),
    ("inline $x^{2}$ math", ["x^{2}"]),
    ("display $$a+b$$ math", ["a+b"]),
    ("both $x$ and $$y$$", ["x", "y"]),
    (r"TeX delimiters \(x\) and \[y\]", ["x", "y"]),
    # Prose: a space follows the opening dollar.
    ("costs $ 5 and $ 6", []),
    # Prose: the closing dollar would be followed by a digit.
    ("from $5 to $6", []),
    # A code span cannot be inside inline math -- without this rule the second
    # dollar opens a span that closes on the one inside the backticks.
    ("prices $5 and $6, where `$HOME` is set", []),
    (r"a \$5 note and \$6 more", []),
    ("```\necho $PATH $HOME\n```\n", []),
    ("$a$\n```\n$PATH\n```\n$b$", ["a", "b"]),
    ("$$\n\\frac{a}{b}\n$$", ["\n\\frac{a}{b}\n"]),
    ("let $u$ be $v$", ["u", "v"]),
    ("one $ alone", []),
]

ROUND_TRIP = [md for md, _ in DETECTION] + [
    "# heading\n\ntext with $x_{i}$ and a list\n\n- $a$\n- `code`\n\n```py\nx = 1\n```\n",
    "trailing newline\n",
    "no trailing newline",
    "",
    "\n\n\n",
    r"backslash \\ and \$ and $y$",
]


@pytest.mark.parametrize("markdown,expect", DETECTION)
def test_detection(markdown, expect):
    doc = MarkdownDoc()
    doc.load(markdown)
    assert [doc.math_latex(i) for i in range(doc.math_count())] == expect


@pytest.mark.parametrize("markdown", ROUND_TRIP)
def test_round_trip_is_byte_exact(markdown):
    doc = MarkdownDoc()
    doc.load(markdown)
    assert doc.text() == markdown


def test_editing_one_equation_leaves_every_other_byte_alone():
    doc = MarkdownDoc()
    doc.load("before $x$ middle $$y$$ after")
    assert doc.set_math_latex(0, "z^{2}")
    assert doc.text() == "before $z^{2}$ middle $$y$$ after"


def test_inline_and_display_are_distinguished():
    doc = MarkdownDoc()
    doc.load("$a$ and $$b$$")
    assert not doc.math_is_display(0)
    assert doc.math_is_display(1)


def test_delimiters_are_preserved_not_normalised():
    """A file written with \\( \\) must not come back as $ $."""
    doc = MarkdownDoc()
    doc.load(r"\(x\)")
    doc.set_math_latex(0, "y")
    assert doc.text() == r"\(y\)"


def test_a_code_span_keeps_its_backticks_out_of_the_body():
    """open/body/close means what the header says it means.

    A viewer renders `body`, so leaving the fences in it prints them on the
    page; every caller having to strip them is the same defect repeated.
    """
    doc = MarkdownDoc()
    doc.load("set `x = 1` here")
    span = next(s for s in doc.segments() if s.kind == MdSegment.CodeSpan)
    assert span.open == "`"
    assert span.body == "x = 1"
    assert span.close == "`"


@pytest.mark.parametrize("markdown", [
    "a `b` c",
    "``a ` b``",
    "`` ` ``",
    "`unclosed",
    "`a`\n`b`",
])
def test_splitting_a_code_span_still_rebuilds_the_file(markdown):
    doc = MarkdownDoc()
    doc.load(markdown)
    assert doc.text() == markdown
