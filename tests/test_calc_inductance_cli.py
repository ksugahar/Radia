"""Fast contracts for the canonical IH inductance command line."""

import re

from radia.panels import calc_inductance


def test_help_renders_without_argparse_percent_interpolation_failure():
    help_text = calc_inductance.build_argparser().format_help()

    # argparse rewraps help at the terminal width, so the phrase may be
    # split after the space or after the hyphen.  Matching the literal
    # substring made this test depend on how many options precede it.
    assert re.search(r"\+25-30%\s+over-\s*estimate", help_text)
    assert "--wp-loop-dof" in help_text
