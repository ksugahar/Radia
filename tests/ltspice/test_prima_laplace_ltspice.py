"""LTspice solves the generated PRIMA Laplace subcircuits to the model's response.

Skipped unless LTspice.exe and an interactive Windows desktop are available.
"""
from __future__ import annotations

import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from radia import prima_export as pe
from radia.ltspice.parser.asc_parser import _find_ltspice_exe, _has_interactive_windows_desktop

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "matlab"))
from prima_port_model_python_reference import (  # noqa: E402
    coupled_rl_admittance, rc_two_port_impedance)

_EXE = _find_ltspice_exe()
pytestmark = pytest.mark.skipif(
    _EXE is None or not _has_interactive_windows_desktop(),
    reason="LTspice.exe with an interactive Windows desktop is required")

FREQUENCIES = (10.0, 300.0, 5e3, 1e5)


def _solve(tmp_path: Path, model, probes: list[str]) -> dict[tuple[int, int], complex]:
    """Drive port 1 with a 1 A AC current (other ports open) and read the probes."""
    library = pe.write_ltspice_subckt(model, tmp_path / "dut.lib", name="DUT")
    pins = " ".join(f"x{k} 0" for k in range(1, model.port_count + 1))
    measures = "\n".join(f".meas AC m{p}_{i} FIND {probe} AT={f}"
                         for p, probe in enumerate(probes) for i, f in enumerate(FREQUENCIES))
    loads = "\n".join(f"Ropen{k} x{k} 0 1e15" for k in range(2, model.port_count + 1))
    netlist = tmp_path / "prima.cir"
    netlist.write_text(
        f"* PRIMA Laplace check\n.include {library.name}\nI1 0 x1 AC 1\nX1 {pins} DUT\n{loads}\n"
        f".ac list {' '.join(map(str, FREQUENCIES))}\n{measures}\n.end\n", encoding="ascii")
    subprocess.run([_EXE, "-b", str(netlist)], cwd=tmp_path, timeout=120, check=True)
    log = netlist.with_suffix(".log").read_text(encoding="utf-8", errors="replace")
    values = {}
    for p, i, db, phase in re.findall(
            r"^m(\d+)_(\d+):.*?=\(([-+0-9.eE]+)dB,([-+0-9.eE]+)\S*\)", log, re.MULTILINE):
        values[(int(p), int(i))] = 10.0 ** (float(db) / 20.0) * np.exp(1j * math.radians(float(phase)))
    assert len(values) == len(probes) * len(FREQUENCIES), log
    return values


def test_coupled_rl_port_impedance(tmp_path):
    model = coupled_rl_admittance()           # improper impedance: realized as Y
    values = _solve(tmp_path, model, ["V(x1)"])
    for i, frequency in enumerate(FREQUENCIES):
        expected = model.inverse().response(2j * math.pi * frequency)[0, 0]
        assert values[(0, i)] == pytest.approx(expected, rel=1e-8)


def test_two_port_rc_impedance_column(tmp_path):
    model = rc_two_port_impedance()
    values = _solve(tmp_path, model, ["V(x1)", "V(x2)"])
    for i, frequency in enumerate(FREQUENCIES):
        Z = model.response(2j * math.pi * frequency)
        assert values[(0, i)] == pytest.approx(Z[0, 0], rel=1e-8)
        assert values[(1, i)] == pytest.approx(Z[1, 0], rel=1e-8)
