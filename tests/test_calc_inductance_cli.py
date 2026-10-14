"""Fast contracts for the canonical IH inductance command line."""

import re
import pytest

from radia.panels import calc_inductance


def test_help_renders_without_argparse_percent_interpolation_failure():
    help_text = calc_inductance.build_argparser().format_help()

    # argparse rewraps help at the terminal width, so the phrase may be
    # split after the space or after the hyphen.  Matching the literal
    # substring made this test depend on how many options precede it.
    assert re.search(r"\+25-30%\s+over-\s*estimate", help_text)
    assert "--wp-loop-dof" in help_text


def test_strong_esim_rejected_before_coil_solve(monkeypatch):
    args = calc_inductance.build_argparser().parse_args([
        '--coil-solver','peec','--frequency','1000','--sigma','5.8e7',
        '--coil-step','not-read.step','--vol','not-read.vol',
        '--coupling-mode','strong','--impedance-model','esim'])
    monkeypatch.setattr(calc_inductance, '_solve_coil_peec',
                        lambda *_: (_ for _ in ()).throw(AssertionError('coil accessed')))
    result = calc_inductance.run_inductance(args)
    assert result['status'] == 'error'
    assert 'ESIM is unsupported' in result['error']


def test_sibc_esim_per_panel_rejected_before_coil_solve(monkeypatch):
    args = calc_inductance.build_argparser().parse_args([
        '--coil-solver','peec','--frequency','1000','--sigma','5.8e7',
        '--coil-step','not-read.step','--vol','not-read.vol',
        '--impedance-model','sibc','--esim-per-panel'])
    monkeypatch.setattr(calc_inductance, '_solve_coil_peec',
                        lambda *_: (_ for _ in ()).throw(AssertionError('coil accessed')))
    result = calc_inductance.run_inductance(args)
    assert result['status'] == 'error'
    assert '--esim-per-panel requires --impedance-model esim' in result['error']


@pytest.mark.parametrize('extra', [[], ['--esim-per-panel', '--h1-order', '2'],
                                  ['--esim-per-panel', '--coupling-mode', 'strong']])
def test_panel_acceleration_rejected_on_unsupported_path_before_coil(monkeypatch, extra):
    args = calc_inductance.build_argparser().parse_args([
        '--coil-solver', 'peec', '--frequency', '1000', '--sigma', '5.8e7',
        '--coil-step', 'not-read.step', '--vol', 'not-read.vol',
        '--impedance-model', 'esim', '--esim-panel-evaluator', 'table', *extra])
    monkeypatch.setattr(calc_inductance, '_solve_coil_peec',
                        lambda *_: (_ for _ in ()).throw(AssertionError('coil accessed')))
    result = calc_inductance.run_inductance(args)
    assert result['status'] == 'error'
    expected = 'ESIM is unsupported' if '--coupling-mode' in extra else 'requires weak P1 per-panel ESIM'
    assert expected in result['error']
