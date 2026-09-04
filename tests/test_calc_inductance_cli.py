"""Fast contracts for the canonical IH inductance command line."""

from radia.panels import calc_inductance


def test_help_renders_without_argparse_percent_interpolation_failure():
    help_text = calc_inductance.build_argparser().format_help()

    assert "+25-30% over-estimate" in help_text
    assert "--wp-loop-dof" in help_text
