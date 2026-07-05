"""Smoke test: every ``panels/calc_*.py`` has to at least ``--help`` cleanly.

Catches import regressions fast (missing deps, syntax errors, bad paths).
Does NOT validate numerical correctness — that's what the per-panel
golden tests (``test_peec_bem_golden.py``, ``test_fem_coilmesh_golden.py``,
``test_peec_inductance_golden.py``) are for.

Runs in < 1 s per script because ``argparse --help`` parses the CLI then
exits before any heavy ``import radia`` work.  Use this as a CI-cheap
regression net; full E2E goldens stay marked ``@pytest.mark.slow`` and
run on LAB + 100号機 only.

Future work: add full numerical E2E goldens for
  * calc_accel_magnet.py  (EM Picard, iron yoke)
  * calc_accel_hdiv.py    (EM HDiv-VIM, iron yoke)
  * calc_pcb_peec.py      (PCB, FastHenry .inp)
  * calc_fem_kelvin.py with --impedance-model esim  (ESIM + Karl iter.)
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parent.parent.parent
PANELS_DIR = ROOT / "src" / "radia" / "panels"


def _calc_scripts() -> list[Path]:
    # calc_common.py is a helper module (calc_main, progress, etc.) not
    # a CLI; it has no argparse and its --help exits 0 with no output.
    return [p for p in sorted(PANELS_DIR.glob("calc_*.py"))
            if p.stem != "calc_common"]


@pytest.mark.parametrize("script", _calc_scripts(), ids=lambda p: p.stem)
def test_calc_help(script: Path):
    """``python calc_*.py --help`` must exit 0 within 30 s."""
    proc = subprocess.run(
        [sys.executable, str(script), "--help"],
        capture_output=True, text=True, timeout=30,
    )
    assert proc.returncode == 0, (
        f"{script.name} --help exited {proc.returncode}\n"
        f"STDOUT:\n{proc.stdout[-500:]}\n"
        f"STDERR:\n{proc.stderr[-500:]}"
    )
    # Sanity: argparse must print at least "usage:" to stdout
    assert "usage" in proc.stdout.lower(), (
        f"{script.name} --help produced no argparse usage line:\n"
        f"{proc.stdout[:500]}"
    )
